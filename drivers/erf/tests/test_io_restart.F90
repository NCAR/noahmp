! ===========================================================================
! test_io_restart -- checkpoint / restart round-trip (NoahmpWriteRestartMod /
! NoahmpReadRestartMod).
!
! The restart writer/reader serialize the FULL prognostic state to a per-level
! NetCDF-4 file using collective parallel I/O (NF90_MPIIO, comm=NoahmpIO%comm),
! at working precision (NF90_DOUBLE when kind_noahmp==8, else NF90_REAL) so a
! restart reproduces a cold run bitwise (see dev/spec-io-restart.md).
!
! The program creates its own MPI world (tio_mpi_init -> MPI_COMM_WORLD) and
! runs single-rank, so one block's hyperslab covers the whole global domain.
!
! Selector argv(1):
!   roundtrip (default) -- fill every allocated prognostic array with a distinct
!       per-variable pattern, WriteRestart, verify the file's precision + layer
!       metadata, zero the arrays, ReadRestart, and assert BIT-EXACT recovery.
!       Registered as a normal (must-pass) CTest case.
!   mismatch            -- write a checkpoint, then corrupt NSOIL and ReadRestart;
!       the layer-geometry guard must abort. Detected by OUTPUT (PASS on the
!       reader's "layer mismatch" diagnostic, FAIL on the "did NOT occur"
!       fall-through below) rather than WILL_FAIL: WILL_FAIL would also go green
!       if WriteRestart aborted first (masking a broken writer without ever
!       reaching the guard). Links test_abort_handler.cpp for a clean exit.
! ===========================================================================
program test_io_restart

  use netcdf
  use Machine,               only : kind_noahmp
  use NoahmpIOVarType,       only : NoahmpIO_type
  use NoahmpWriteRestartMod, only : NoahmpWriteRestart
  use NoahmpReadRestartMod,  only : NoahmpReadRestart
  use NoahmpTestIOSupport

  implicit none

  ! visit modes for the single-source-of-truth state list
  integer, parameter :: FILL = 1, ZERO = 2, CHECK = 3

  type(NoahmpIO_type), pointer :: blk
  integer :: comm
  integer, parameter :: NX = 4, NY = 3, NSOIL = 4, NSNOW = 3
  character(len=32)  :: which
  character(len=64)  :: dir

  call get_command_argument(1, which)
  if (len_trim(which) == 0) which = "roundtrip"

  call tio_reset()
  call tio_mpi_init(comm)
  call io_setup_block(blk, level=0, blkid=0, nblocks=1, nx=NX, ny=NY, &
                      nsoil=NSOIL, nsnow=NSNOW, comm=comm)

  select case (trim(which))
  case ("mismatch")
     dir = "rst_mismatch"
     call visit_state(blk, FILL)
     call NoahmpWriteRestart(blk, trim(dir), 1)
     ! Corrupt the run's layer count: the reader must reject the file and abort.
     blk%NSOIL = NSOIL + 1
     call NoahmpReadRestart(blk, trim(dir), 1)
     ! Unreachable if the guard fired.
     write(0,'(A)') "test_io_restart: expected layer-mismatch abort did NOT occur"
     call exit(0)

  case ("roundtrip")
     dir = "rst_roundtrip"
     call visit_state(blk, FILL)
     call NoahmpWriteRestart(blk, trim(dir), 1)

     ! Checkpoint file contract: precision matches kind, layer metadata present.
     call verify_file_meta(trim(dir))

     ! Wipe the in-memory state, restore it, and require bit-exact recovery.
     call visit_state(blk, ZERO)
     call NoahmpReadRestart(blk, trim(dir), 1)
     call visit_state(blk, CHECK)

     call tio_mpi_finalize()
     call tio_finish("test_io_restart")

  case default
     write(0,'(A)') "test_io_restart: unknown selector '"//trim(which)//"'"
     call exit(2)
  end select

contains

  ! Deterministic per-(variable,index) reference values. Any value works for the
  ! bit-exact round trip (no type conversion on write/read); distinct tags guard
  ! against the reader/writer transposing two variables.
  pure function expr_r(tag, i, k, j) result(v)
    integer, intent(in) :: tag, i, k, j
    real(kind=kind_noahmp) :: v
    v = real(tag*100000 + i*1000 + k*100 + j, kind_noahmp)
  end function expr_r

  pure function expr_i(tag, i, k, j) result(v)
    integer, intent(in) :: tag, i, k, j
    integer :: v
    v = tag*100000 + i*1000 + k*100 + j
  end function expr_i

  ! One list, three behaviors (fill / zero / check), so the array set can never
  ! drift between filling and verifying. Optional fields are handled uniformly:
  ! skipped in every mode when unallocated (matching the writer's allocated()
  ! guard), so they neither get written nor asserted.
  subroutine visit_state(b, mode)
    type(NoahmpIO_type), intent(inout) :: b
    integer,             intent(in)    :: mode

    ! --- soil (3D) ---
    call op3d(b%TSLB,    1, mode, "TSLB")
    call op3d(b%SMOIS,   2, mode, "SMOIS")
    call op3d(b%SH2O,    3, mode, "SH2O")
    if (allocated(b%SMOISEQ)) call op3d(b%SMOISEQ, 4, mode, "SMOISEQ")
    ! --- snow layers (3D, negative-indexed) ---
    call op3d(b%TSNOXY,  5, mode, "TSNOXY")
    call op3d(b%SNICEXY, 6, mode, "SNICEXY")
    call op3d(b%SNLIQXY, 7, mode, "SNLIQXY")
    call op3d(b%ZSNSOXY, 8, mode, "ZSNSOXY")
    ! --- snowpack scalars (2D) + ISNOWXY (int) ---
    call op2d(b%SNOW,    10, mode, "SNOW")
    call op2d(b%SNOWH,   11, mode, "SNOWH")
    call op2d(b%SNOWC,   12, mode, "SNOWC")
    call op2d(b%CANWAT,  13, mode, "CANWAT")
    call op2d(b%ACSNOM,  14, mode, "ACSNOM")
    call op2d(b%ACSNOW,  15, mode, "ACSNOW")
    call op2di(b%ISNOWXY, 70, mode, "ISNOWXY")
    ! --- canopy / surface (2D) ---
    call op2d(b%TVXY,    16, mode, "TVXY")
    call op2d(b%TGXY,    17, mode, "TGXY")
    call op2d(b%CANICEXY,18, mode, "CANICEXY")
    call op2d(b%CANLIQXY,19, mode, "CANLIQXY")
    call op2d(b%EAHXY,   20, mode, "EAHXY")
    call op2d(b%TAHXY,   21, mode, "TAHXY")
    call op2d(b%CMXY,    22, mode, "CMXY")
    call op2d(b%CHXY,    23, mode, "CHXY")
    call op2d(b%FWETXY,  24, mode, "FWETXY")
    call op2d(b%QSFC,    25, mode, "QSFC")
    call op2d(b%TSK,     60, mode, "TSK")        ! c_kind_noahmp boundary field
    call op2d(b%QSNOWXY, 26, mode, "QSNOWXY")
    call op2d(b%QRAINXY, 27, mode, "QRAINXY")
    ! --- albedo history (2D) ---
    call op2d(b%SNEQVOXY,28, mode, "SNEQVOXY")
    call op2d(b%ALBOLDXY,29, mode, "ALBOLDXY")
    call op2d(b%TAUSSXY, 30, mode, "TAUSSXY")
    call op2d(b%ALBEDO,  31, mode, "ALBEDO")
    ! --- aquifer / groundwater (2D) ---
    call op2d(b%ZWTXY,     32, mode, "ZWTXY")
    call op2d(b%WAXY,      33, mode, "WAXY")
    call op2d(b%WTXY,      34, mode, "WTXY")
    call op2d(b%SMCWTDXY,  35, mode, "SMCWTDXY")
    call op2d(b%DEEPRECHXY,36, mode, "DEEPRECHXY")
    call op2d(b%RECHXY,    37, mode, "RECHXY")
    ! --- phenology (2D) ---
    call op2d(b%LAI,     38, mode, "LAI")
    call op2d(b%XSAIXY,  39, mode, "XSAIXY")
    ! --- accumulators / carried state (2D) ---
    call op2d(b%SFCRUNOFF,40, mode, "SFCRUNOFF")
    call op2d(b%UDRUNOFF, 41, mode, "UDRUNOFF")
    call op2d(b%SMSTAV,   42, mode, "SMSTAV")
    call op2d(b%SMSTOT,   43, mode, "SMSTOT")
    call op2d(b%EMISS,    61, mode, "EMISS")     ! c_kind_noahmp boundary field
    call op2d(b%GRDFLX,   44, mode, "GRDFLX")
    ! --- soil-cycle accumulators (SOIL_UPDATE_STEPS>1 carry; ETRANI is 3D over soil) ---
    call op2d(b%ACC_SSOILXY,  71, mode, "ACC_SSOILXY")
    call op2d(b%ACC_QINSURXY, 72, mode, "ACC_QINSURXY")
    call op2d(b%ACC_QSEVAXY,  73, mode, "ACC_QSEVAXY")
    call op2d(b%ACC_DWATERXY, 74, mode, "ACC_DWATERXY")
    call op2d(b%ACC_PRCPXY,   75, mode, "ACC_PRCPXY")
    call op2d(b%ACC_ECANXY,   76, mode, "ACC_ECANXY")
    call op2d(b%ACC_ETRANXY,  77, mode, "ACC_ETRANXY")
    call op2d(b%ACC_EDIRXY,   78, mode, "ACC_EDIRXY")
    call op3d(b%ACC_ETRANIXY, 79, mode, "ACC_ETRANIXY")
    call op2d(b%ACC_GLAFLWXY, 80, mode, "ACC_GLAFLWXY")
    ! --- optional carbon / dveg / lake (only if allocated) ---
    if (allocated(b%LFMASSXY)) call op2d(b%LFMASSXY, 45, mode, "LFMASSXY")
    if (allocated(b%RTMASSXY)) call op2d(b%RTMASSXY, 46, mode, "RTMASSXY")
    if (allocated(b%STMASSXY)) call op2d(b%STMASSXY, 47, mode, "STMASSXY")
    if (allocated(b%WOODXY))   call op2d(b%WOODXY,   48, mode, "WOODXY")
    if (allocated(b%GRAINXY))  call op2d(b%GRAINXY,  49, mode, "GRAINXY")
    if (allocated(b%GDDXY))    call op2d(b%GDDXY,    50, mode, "GDDXY")
    if (allocated(b%WSLAKEXY)) call op2d(b%WSLAKEXY, 62, mode, "WSLAKEXY")
  end subroutine visit_state

  ! Assumed-shape dummies index 1..size, so negative snow lower bounds are handled
  ! transparently and fill/check use the identical convention.
  subroutine op3d(arr, tag, mode, name)
    real(kind=kind_noahmp), intent(inout) :: arr(:,:,:)
    integer,      intent(in) :: tag, mode
    character(*), intent(in) :: name
    integer :: i, k, j
    real(kind=kind_noahmp) :: want
    do j = 1, size(arr,3)
       do k = 1, size(arr,2)
          do i = 1, size(arr,1)
             select case (mode)
             case (FILL);  arr(i,k,j) = expr_r(tag,i,k,j)
             case (ZERO);  arr(i,k,j) = 0.0_kind_noahmp
             case (CHECK)
                want = expr_r(tag,i,k,j)
                if (arr(i,k,j) /= want) call fail3(name, i, k, j)
             end select
          end do
       end do
    end do
  end subroutine op3d

  subroutine op2d(arr, tag, mode, name)
    real(kind=kind_noahmp), intent(inout) :: arr(:,:)
    integer,      intent(in) :: tag, mode
    character(*), intent(in) :: name
    integer :: i, j
    real(kind=kind_noahmp) :: want
    do j = 1, size(arr,2)
       do i = 1, size(arr,1)
          select case (mode)
          case (FILL);  arr(i,j) = expr_r(tag,i,0,j)
          case (ZERO);  arr(i,j) = 0.0_kind_noahmp
          case (CHECK)
             want = expr_r(tag,i,0,j)
             if (arr(i,j) /= want) call fail3(name, i, 0, j)
          end select
       end do
    end do
  end subroutine op2d

  subroutine op2di(arr, tag, mode, name)
    integer,      intent(inout) :: arr(:,:)
    integer,      intent(in) :: tag, mode
    character(*), intent(in) :: name
    integer :: i, j
    do j = 1, size(arr,2)
       do i = 1, size(arr,1)
          select case (mode)
          case (FILL);  arr(i,j) = expr_i(tag,i,0,j)
          case (ZERO);  arr(i,j) = 0
          case (CHECK)
             if (arr(i,j) /= expr_i(tag,i,0,j)) call fail3(name, i, 0, j)
          end select
       end do
    end do
  end subroutine op2di

  subroutine fail3(name, i, k, j)
    character(*), intent(in) :: name
    integer,      intent(in) :: i, k, j
    call tio_expect(.false., "restart roundtrip mismatch: "//trim(name))
    write(0,'(A,A,A,I0,A,I0,A,I0,A)') "    (", trim(name), " at [", i, ",", k, ",", j, "])"
  end subroutine fail3

  ! Reopen the checkpoint serially and assert the precision + layer-geometry
  ! contract that guarantees a bit-exact restart.
  subroutine verify_file_meta(d)
    character(*), intent(in) :: d
    integer :: ncid, vid, fn_soil, fn_snow, xt, expect_rtype
    call ncchk(nf90_open(trim(d)//"/Level_0.nc", NF90_NOWRITE, ncid), "reopen checkpoint")
    call ncchk(nf90_get_att(ncid, NF90_GLOBAL, "NSOIL", fn_soil), "get NSOIL attr")
    call ncchk(nf90_get_att(ncid, NF90_GLOBAL, "NSNOW", fn_snow), "get NSNOW attr")
    call tio_expect_eq_i(fn_soil, NSOIL, "checkpoint NSOIL attribute")
    call tio_expect_eq_i(fn_snow, NSNOW, "checkpoint NSNOW attribute")

    expect_rtype = NF90_REAL
    if (kind_noahmp == 8) expect_rtype = NF90_DOUBLE
    call ncchk(nf90_inq_varid(ncid, "TSLB", vid), "inq TSLB")
    call ncchk(nf90_inquire_variable(ncid, vid, xtype=xt), "xtype TSLB")
    call tio_expect_eq_i(xt, expect_rtype, "TSLB stored at working precision")

    call ncchk(nf90_inq_varid(ncid, "ISNOWXY", vid), "inq ISNOWXY")
    call ncchk(nf90_inquire_variable(ncid, vid, xtype=xt), "xtype ISNOWXY")
    call tio_expect_eq_i(xt, NF90_INT, "ISNOWXY stored as NF90_INT")
    call ncchk(nf90_close(ncid), "close checkpoint")
  end subroutine verify_file_meta

end program test_io_restart
