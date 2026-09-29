! ===========================================================================
! test_fortran_alloc -- pure-Fortran storage + allocation (Tier 1).
!
! Exercises the Fortran side of the boundary directly (no C++): size the block
! array via the bind(C) entry point NoahmpIOTypeVectInit_fi, populate a block's
! dimension/option scalars, run the generated NoahmpIOVarInitDefault allocator,
! and verify the resulting array shapes/bounds. Also checks the precision probe
! (NoahmpRealSize_fi) and the idempotent re-init guard. Exit code 0 = pass.
! ===========================================================================
program test_fortran_alloc

  use iso_c_binding,     only : C_INT, C_SIZE_T
  use Machine,           only : c_kind_noahmp
  use NoahmpIO_fi,       only : NoahmpIO_vect, NoahmpIOTypeVectInit_fi, &
                                NoahmpRealSize_fi
  use NoahmpIOVarType,   only : NoahmpIO_type
  use NoahmpIOVarInitMod, only : NoahmpIOVarInitDefault

  implicit none

  integer(C_INT) :: level, nblocks
  type(NoahmpIO_type), pointer :: blk
  integer :: nfail
  integer, parameter :: XS=1, XE=4, YS=1, YE=3, KMS=1, KME=3, NSOIL=4, NSNOW=3, NUMRAD=2
  real(kind=c_kind_noahmp) :: probe
  integer(C_SIZE_T) :: rsize

  nfail = 0
  level = 0
  nblocks = 1

  ! ---- size the module-global block array once ----
  call NoahmpIOTypeVectInit_fi(level, nblocks)
  call expect(allocated(NoahmpIO_vect(level)%NoahmpIO), "level allocated")
  call expect_eq_i(size(NoahmpIO_vect(level)%NoahmpIO), nblocks, "nblocks")

  blk => NoahmpIO_vect(level)%NoahmpIO(0)

  ! ---- dimension scalars ----
  ! The coupled dimension scalars are C++-OWNED pointer components (normally
  ! wired to C++ memory by ScalarInitDefault). In this pure-Fortran path there is
  ! no C++ owner, so we give each pointer its own Fortran target before use. The
  ! IOPT_* options are ordinary (host-only) integers, assigned directly.
  allocate(blk%XSTART); blk%XSTART = XS
  allocate(blk%XEND);   blk%XEND   = XE
  allocate(blk%YSTART); blk%YSTART = YS
  allocate(blk%YEND);   blk%YEND   = YE
  allocate(blk%KMS);    blk%KMS    = KMS
  allocate(blk%KME);    blk%KME    = KME
  allocate(blk%NSOIL);  blk%NSOIL  = NSOIL
  allocate(blk%NSNOW);  blk%NSNOW  = NSNOW
  allocate(blk%NUMRAD); blk%NUMRAD = NUMRAD
  allocate(blk%ITIMESTEP)   ! NoahmpIOVarInitDefault initializes this coupled scalar
  blk%IOPT_SOIL = 1; blk%IOPT_ALB = 1; blk%IOPT_WETLAND = 0; blk%SF_URBAN_PHYSICS = 0

  ! ---- allocate all storage ----
  call NoahmpIOVarInitDefault(blk)

  ! ---- 3-D soil arrays: [XSTART:XEND, 1:NSOIL, YSTART:YEND] ----
  call check_shape3(blk%SMOIS, XE-XS+1, NSOIL, YE-YS+1, "SMOIS shape")
  call expect_eq_i(lbound(blk%SMOIS,2), 1,     "SMOIS lbound layer")
  call expect_eq_i(ubound(blk%SMOIS,2), NSOIL, "SMOIS ubound layer")
  call check_shape3(blk%TSLB,  XE-XS+1, NSOIL, YE-YS+1, "TSLB shape")

  ! ---- 1-D soil arrays: DZS(1:NSOIL) ----
  call expect_eq_i(size(blk%DZS),   NSOIL, "DZS size")
  call expect_eq_i(lbound(blk%DZS,1), 1,   "DZS lbound")

  ! ---- 2-D field: XLAT[XSTART:XEND, YSTART:YEND] ----
  call expect_eq_i(size(blk%XLAT,1), XE-XS+1, "XLAT dim1")
  call expect_eq_i(size(blk%XLAT,2), YE-YS+1, "XLAT dim2")

  ! ---- precision probe agrees with the Fortran coupling real kind ----
  rsize = NoahmpRealSize_fi()
  call expect(rsize == 4_C_SIZE_T .or. rsize == 8_C_SIZE_T, "RealSize 4 or 8")
  call expect_eq_i(int(rsize), storage_size(probe)/8, "RealSize == sizeof(real)")

  ! ---- idempotent re-init with the SAME size is a no-op, not an abort ----
  call NoahmpIOTypeVectInit_fi(level, nblocks)
  call expect_eq_i(size(NoahmpIO_vect(level)%NoahmpIO), nblocks, "re-init same size")

  if (nfail == 0) then
     write(*,'(A)') "PASS: test_fortran_alloc"
     call exit(0)
  else
     write(*,'(A,I0,A)') "FAILED: test_fortran_alloc (", nfail, " check(s))"
     call exit(1)
  end if

contains

  subroutine expect(cond, name)
    logical, intent(in) :: cond
    character(*), intent(in) :: name
    if (.not. cond) then
       write(0,'(A,A)') "  FAIL: ", name
       nfail = nfail + 1
    end if
  end subroutine expect

  subroutine expect_eq_i(got, want, name)
    integer, intent(in) :: got, want
    character(*), intent(in) :: name
    if (got /= want) then
       write(0,'(A,A,A,I0,A,I0)') "  FAIL: ", name, " got=", got, " want=", want
       nfail = nfail + 1
    end if
  end subroutine expect_eq_i

  subroutine check_shape3(a, n1, n2, n3, name)
    real(kind=c_kind_noahmp), allocatable, intent(in) :: a(:,:,:)
    integer, intent(in) :: n1, n2, n3
    character(*), intent(in) :: name
    call expect(allocated(a), name//" allocated")
    if (allocated(a)) then
       call expect_eq_i(size(a,1), n1, name//" d1")
       call expect_eq_i(size(a,2), n2, name//" d2")
       call expect_eq_i(size(a,3), n3, name//" d3")
    end if
  end subroutine check_shape3

end program test_fortran_alloc
