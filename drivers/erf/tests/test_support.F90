! ===========================================================================
! Fortran test-support module for the Noah-MP ERF-driver regression tests.
!
! These bind(C) helpers are compiled INTO the test executables (never into
! libnoahmp) and reach into the same module-global block array the C++ side
! drives through NoahmpIO_type -- i.e. NoahmpIO_vect(level)%NoahmpIO(blkid).
! They let a C++ test:
!   * set the branch scalars (IOPT_*) that NoahmpIOVarInitDefault reads but that
!     are normally populated by ReadNamelist -- so VarInitDefault is safe to call
!     without a namelist file (noahmp_test_prep_block),
!   * read/write Fortran-owned arrays and table scalars to prove zero-copy memory
!     sharing across the C<->Fortran boundary and that ReadTable populated state.
!
! Everything uses real(c_kind_noahmp)/noahmp_real so it stays precision-correct
! whether the library was built single or double precision.
! ===========================================================================
module NoahmpTestSupport

  use iso_c_binding, only : C_INT
  use Machine,       only : c_kind_noahmp
  use NoahmpIO_fi,   only : NoahmpIO_vect, NLEVEL_MIN, NLEVEL_MAX
  use NoahmpIOVarType, only : NoahmpIO_type

  implicit none
  private

  public :: noahmp_test_prep_block
  public :: noahmp_test_read_xlat, noahmp_test_read_smois
  public :: noahmp_test_write_hfx, noahmp_test_read_table

contains

  ! Resolve a pointer to the module-global block, aborting loudly on a bad index
  ! (a test bug, not a library condition).
  subroutine get_block(level, blkid, blk)
    integer(C_INT),               intent(in)  :: level, blkid
    type(NoahmpIO_type), pointer, intent(out) :: blk
    blk => null()
    if (level < NLEVEL_MIN .or. level > NLEVEL_MAX) then
       write(0,*) "test_support: level out of range:", level
       error stop 2
    end if
    if (.not. allocated(NoahmpIO_vect(level)%NoahmpIO)) then
       write(0,*) "test_support: level not initialized:", level
       error stop 2
    end if
    if (blkid < lbound(NoahmpIO_vect(level)%NoahmpIO,1) .or. &
        blkid > ubound(NoahmpIO_vect(level)%NoahmpIO,1)) then
       write(0,*) "test_support: blkid out of range:", blkid
       error stop 2
    end if
    blk => NoahmpIO_vect(level)%NoahmpIO(blkid)
  end subroutine get_block

  ! Set the config-option scalars that gate allocations in NoahmpIOVarInitDefault
  ! to benign values, standing in for a ReadNamelist call. IOPT_SOIL=1/IOPT_ALB=1
  ! skip the soil-composition, SNICAR and urban/wetland allocation branches.
  subroutine noahmp_test_prep_block(level, blkid) bind(C, name="noahmp_test_prep_block")
    integer(C_INT), value, intent(in) :: level, blkid
    type(NoahmpIO_type), pointer :: blk
    call get_block(level, blkid, blk)
    blk%IOPT_SOIL        = 1
    blk%IOPT_ALB         = 1
    blk%IOPT_WETLAND     = 0
    blk%SF_URBAN_PHYSICS = 0
  end subroutine noahmp_test_prep_block

  ! Read a 2-D coupled array element (latitude) straight from Fortran storage.
  function noahmp_test_read_xlat(level, blkid, i, j) result(val) &
       bind(C, name="noahmp_test_read_xlat")
    integer(C_INT), value, intent(in) :: level, blkid, i, j
    real(kind=c_kind_noahmp) :: val
    type(NoahmpIO_type), pointer :: blk
    call get_block(level, blkid, blk)
    val = blk%XLAT(i, j)
  end function noahmp_test_read_xlat

  ! Read a 3-D coupled array element (soil moisture) -- verifies column-major
  ! (i, layer, j) index agreement across the boundary.
  function noahmp_test_read_smois(level, blkid, i, k, j) result(val) &
       bind(C, name="noahmp_test_read_smois")
    integer(C_INT), value, intent(in) :: level, blkid, i, k, j
    real(kind=c_kind_noahmp) :: val
    type(NoahmpIO_type), pointer :: blk
    call get_block(level, blkid, blk)
    val = blk%SMOIS(i, k, j)
  end function noahmp_test_read_smois

  ! Write a 2-D coupled array element (sensible heat flux) from Fortran so the
  ! C++ view can read it back -- proves the Fortran->C++ direction of sharing.
  subroutine noahmp_test_write_hfx(level, blkid, i, j, val) &
       bind(C, name="noahmp_test_write_hfx")
    integer(C_INT), value, intent(in) :: level, blkid, i, j
    real(kind=c_kind_noahmp), value, intent(in) :: val
    type(NoahmpIO_type), pointer :: blk
    call get_block(level, blkid, blk)
    blk%HFX(i, j) = val
  end subroutine noahmp_test_write_hfx

  ! Return a Fortran-owned table scalar populated by ReadTable so a C++ test can
  ! confirm the parse landed. which: 0 -> ZBOT_TABLE, 1 -> CSOIL_TABLE.
  function noahmp_test_read_table(level, blkid, which) result(val) &
       bind(C, name="noahmp_test_read_table")
    integer(C_INT), value, intent(in) :: level, blkid, which
    real(kind=c_kind_noahmp) :: val
    type(NoahmpIO_type), pointer :: blk
    call get_block(level, blkid, blk)
    select case (which)
    case (0)
       val = blk%ZBOT_TABLE
    case default
       val = blk%CSOIL_TABLE
    end select
  end function noahmp_test_read_table

end module NoahmpTestSupport
