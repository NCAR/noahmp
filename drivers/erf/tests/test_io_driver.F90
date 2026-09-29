! ===========================================================================
! test_io_driver -- full cold-init chain + NoahmpDriverMain over two steps.
!
! Ties the NetCDF read path to the physics driver: it reproduces ERF's exact
! per-block init sequence (ERF_NOAHMP_Init.cpp)
!
!   ReadNamelist -> ReadLandHeader -> VarInitDefault -> ReadTable
!                -> ReadLandMain -> InitMain
!
! from a wrfinput/WPS file, stages benign atmospheric forcing, then advances
! NoahmpDriverMain for TWO iterations (itimestep = 1, 2) -- the itimestep==1
! initial-guess branch followed by a normal step -- and asserts the run reaches
! the end and leaves physically plausible, finite surface state on every land
! cell.
!
! Like the read test it is generic: it synthesizes its own small wrfinput by
! default (all cells land, non-glacier veg/soil categories), or reads any real
! local file given via NOAHMP_TEST_WRFINPUT. The program owns the MPI world.
!
! namelist.erf and NoahmpTable.TBL are staged into the working directory by
! CMake (shared with config_smoke). Exit code 0 = pass.
! ===========================================================================
program test_io_driver

  use iso_c_binding,        only : C_INT
  use netcdf
  use Machine,              only : kind_noahmp, c_kind_noahmp
  use NoahmpIO_fi,          only : NoahmpIO_vect, NoahmpIOTypeVectInit_fi
  use NoahmpIOVarType,      only : NoahmpIO_type
  use NoahmpIOVarInitMod,   only : NoahmpIOVarInitDefault
  use NoahmpReadNamelistMod,only : NoahmpReadNamelist
  use NoahmpReadTableMod,   only : NoahmpReadTable
  use NoahmpReadLandMod,    only : NoahmpReadLandHeader, NoahmpReadLandMain
  use NoahmpInitMainMod,    only : NoahmpInitMain
  use NoahmpDriverMainMod,  only : NoahmpDriverMain
  use NoahmpTestIOSupport

  implicit none

  type(NoahmpIO_type), pointer :: blk
  integer :: comm
  integer :: nx, ny, nsoil
  integer :: it
  logical :: synthetic
  character(len=1024) :: ext_path, path
  integer :: elen, estat
  real(kind=c_kind_noahmp), allocatable :: tsk0(:,:)   ! TSK snapshot before stepping

  call tio_reset()
  call tio_mpi_init(comm)

  ! Fixture selection (synthetic by default; any real WPS file via the env var).
  call get_environment_variable("NOAHMP_TEST_WRFINPUT", ext_path, elen, estat)
  synthetic = (estat /= 0 .or. elen == 0)
  if (synthetic) then
     nx = 4; ny = 4; nsoil = 4
     path = "test_wrfinput_driver.nc"
     call make_wrfinput(trim(path), nx, ny, nsoil)
  else
     path = trim(ext_path)
     call query_dims(trim(path), nx, ny, nsoil)
     write(*,'(A,A,A,I0,A,I0,A,I0)') "Using external wrfinput: ", trim(path), &
          "  nx=", nx, " ny=", ny, " nsoil=", nsoil
  end if

  ! ---- size + wire one block spanning the whole domain (mirrors ERF init) ----
  call NoahmpIOTypeVectInit_fi(0_C_INT, 1_C_INT)
  blk => NoahmpIO_vect(0)%NoahmpIO(0)

  ! Allocate the coupled pointer scalars (C++-owned in production). ReadNamelist
  ! writes xstart/xend/ystart/yend from the namelist (defaults 0/-1), and
  ! ReadLandHeader reads the header at (xstart,ystart); the REAL tile bounds are
  ! set afterwards, exactly as ERF's per-block init does (ERF_NOAHMP_Init.cpp).
  call seti(blk%XSTART, 0);      call seti(blk%XEND, -1)
  call seti(blk%YSTART, 0);      call seti(blk%YEND, -1)
  call seti(blk%IDS, 0);         call seti(blk%IDE, nx-1)
  call seti(blk%JDS, 0);         call seti(blk%JDE, ny-1)
  call seti(blk%KDS, 1);         call seti(blk%KDE, 2)
  call seti(blk%ITS, 0);         call seti(blk%ITE, nx-1)
  call seti(blk%JTS, 0);         call seti(blk%JTE, ny-1)
  call seti(blk%KTS, 1);         call seti(blk%KTE, 2)
  call seti(blk%IMS, 0);         call seti(blk%IME, nx-1)
  call seti(blk%JMS, 0);         call seti(blk%JME, ny-1)
  call seti(blk%KMS, 1);         call seti(blk%KME, 2)
  call seti(blk%NSOIL, nsoil);   call seti(blk%NSNOW, 4)
  call seti(blk%NUMRAD, 2);      call seti(blk%RANK, 0)
  call seti(blk%BLKID, 0);       call seti(blk%LEVEL, 0)
  call seti(blk%COMM, comm)
  allocate(blk%ITIMESTEP)
  call seti(blk%NTIME, 0)                       ! coupled pointer; InitMain sets it
  call setr(blk%DTBL, 3600.0_c_kind_noahmp)
  call setr(blk%ZLVL, 10.0_c_kind_noahmp)

  blk%YR     = 2023
  blk%JULIAN = 229.0_kind_noahmp   ! ~Aug 17, matching namelist start date

  ! ---- ERF cold-init sequence ----
  call NoahmpReadNamelist(blk)                 ! options, nsoil/nsnow, DTBL, soil layers
  blk%erf_setup_file_lev = trim(path)          ! namelist path is a placeholder; override
  blk%erf_setup_file_01  = trim(path)
  call NoahmpReadLandHeader(blk)               ! grid globals + offsets from the file

  ! Real tile/domain/memory bounds span the whole block (ERF sets these AFTER the
  ! header read; ReadNamelist's 0/-1 placeholders would otherwise size arrays to
  ! zero). xsglobal/xoffset were just set from the file by ReadLandHeader.
  blk%XSTART = 0; blk%XEND = nx-1; blk%YSTART = 0; blk%YEND = ny-1
  blk%IDS = 0; blk%IDE = nx-1; blk%JDS = 0; blk%JDE = ny-1
  blk%ITS = 0; blk%ITE = nx-1; blk%JTS = 0; blk%JTE = ny-1
  blk%IMS = 0; blk%IME = nx-1; blk%JMS = 0; blk%JME = ny-1

  call NoahmpIOVarInitDefault(blk)             ! allocate all storage
  call NoahmpReadTable(blk)                    ! parameter tables
  call NoahmpReadLandMain(blk)                 ! land/soil state from the file
  call NoahmpInitMain(blk)                     ! cold-init the prognostic state

  ! Guard against a vacuous pass: the state arrays must actually span the domain
  ! (a zero-size allocation would make every all() assertion below trivially true).
  call tio_expect(size(blk%TSK,1) == nx .and. size(blk%TSK,2) == ny, &
                  "state arrays span the full domain (non-empty)")

  ! ---- benign atmospheric forcing (lowest level; driver mirrors to level 2) ----
  blk%T_PHY(:,1,:)   = 290.0_c_kind_noahmp
  blk%QV_CURR(:,1,:) = 0.006_c_kind_noahmp
  blk%U_PHY(:,1,:)   = 3.0_c_kind_noahmp
  blk%V_PHY(:,1,:)   = 1.0_c_kind_noahmp
  blk%P8W(:,1,:)     = 1.0e5_c_kind_noahmp
  blk%SWDOWN         = 400.0_c_kind_noahmp
  blk%GLW            = 340.0_c_kind_noahmp
  blk%COSZEN         = 0.6_c_kind_noahmp
  blk%RAINBL         = 0.0_c_kind_noahmp
  blk%SR             = 0.0_c_kind_noahmp
  blk%MP_RAINNC      = 0.0_c_kind_noahmp
  blk%MP_SNOW        = 0.0_c_kind_noahmp
  blk%MP_GRAUP       = 0.0_c_kind_noahmp
  blk%MP_HAIL        = 0.0_c_kind_noahmp

  ! Snapshot TSK so we can prove the driver actually produced output. Every
  ! asserted "plausible" value below is already true at cold-init (TSK is read
  ! from the file at ~290 K, inside [150,400]; HFX/TSLB start finite), so a driver
  ! that silently did nothing would pass all of them -- this delta catches that.
  tsk0 = blk%TSK

  ! ---- advance two steps: itimestep==1 (initial guess) then a normal step ----
  do it = 1, 2
     blk%ITIMESTEP = it
     call NoahmpDriverMain(blk)
     ! State must stay finite after each step. abs(x) < huge is false for both NaN
     ! and +/-Inf, so it is a single finiteness test. (Exact magnitudes are not
     ! asserted: with the crude constant forcing here they are not meaningful --
     ! the point is that the read -> init -> driver pipeline runs without blowing
     ! up. TSK additionally must stay in a broad physical range.)
     call tio_expect(all(abs(blk%HFX)  < huge(1.0_c_kind_noahmp)), "HFX finite")
     call tio_expect(all(abs(blk%TSLB) < huge(1.0_c_kind_noahmp)), "TSLB finite")
     call tio_expect(all(abs(blk%TSK)  < huge(1.0_c_kind_noahmp)), "TSK finite")
     call tio_expect(all(blk%TSK > 150.0_c_kind_noahmp .and. blk%TSK < 400.0_c_kind_noahmp), &
                     "TSK physically plausible [K]")
     write(*,'(A,I0,A)') "  NoahmpDriverMain iteration ", it, " completed"
  end do

  ! Not a no-op: the driver must have written surface temperature (a driver that
  ! returned without updating TSK would leave it bit-identical to the cold-init
  ! read and still satisfy every finiteness/range check above).
  call tio_expect(any(blk%TSK /= tsk0), "driver updated TSK (not a no-op)")

  call tio_mpi_finalize()
  call tio_finish("test_io_driver")

contains

  ! allocate + set an integer(C_INT) coupled pointer scalar
  subroutine seti(p, v)
    integer(C_INT), pointer, intent(out) :: p
    integer,                 intent(in)  :: v
    allocate(p)
    p = v
  end subroutine seti

  ! allocate + set a real coupled pointer scalar
  subroutine setr(p, v)
    real(kind=c_kind_noahmp), pointer, intent(out) :: p
    real(kind=c_kind_noahmp),          intent(in)  :: v
    allocate(p)
    p = v
  end subroutine setr

  ! Read nx/ny/nsoil from a wrfinput/WPS file's standard dimensions.
  subroutine query_dims(fpath, qnx, qny, qnsoil)
    character(*), intent(in)  :: fpath
    integer,      intent(out) :: qnx, qny, qnsoil
    integer :: ncid, did
    call ncchk(nf90_open(fpath, NF90_NOWRITE, ncid), "open external wrfinput")
    call ncchk(nf90_inq_dimid(ncid, "west_east", did), "dimid west_east")
    call ncchk(nf90_inquire_dimension(ncid, did, len=qnx), "len west_east")
    call ncchk(nf90_inq_dimid(ncid, "south_north", did), "dimid south_north")
    call ncchk(nf90_inquire_dimension(ncid, did, len=qny), "len south_north")
    call ncchk(nf90_inq_dimid(ncid, "soil_layers_stag", did), "dimid soil_layers_stag")
    call ncchk(nf90_inquire_dimension(ncid, did, len=qnsoil), "len soil_layers_stag")
    call ncchk(nf90_close(ncid), "close external wrfinput")
  end subroutine query_dims

end program test_io_driver
