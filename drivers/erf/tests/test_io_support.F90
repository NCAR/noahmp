! ===========================================================================
! test_io_support -- shared fixtures/helpers for the NetCDF I/O regression
! tests (read / write / restart / checkpoint) of the Noah-MP ERF driver.
!
! Compiled INTO the I/O test executables (never into libnoahmp), like
! test_support.F90, but with a private module dir to avoid .mod clashes.
!
! The read path (NoahmpReadLandMod) opens a real wrfinput/WPS-style NetCDF
! file. The local ChisholmView fixture is >200 MB and not committed, so instead
! of hard-coding it these helpers SYNTHESIZE a tiny, self-contained wrfinput
! (make_wrfinput) carrying every global attribute and variable the reader
! needs. The generated file is therefore generic -- the read test works with no
! external data -- while still exercising the real reader end to end. Set
! NOAHMP_TEST_WRFINPUT to additionally point the header reader at any real local
! WPS/wrfinput file (see test_io_readland).
!
! The write/restart paths (NoahmpWriteLandMod / Noahmp*RestartMod) use collective
! parallel NetCDF-4 (NF90_MPIIO, comm=NoahmpIO%comm). The test PROGRAM creates
! its own MPI world (tio_mpi_init -> MPI_COMM_WORLD) and runs single-rank; the
! block's hyperslab then covers the whole global domain.
!
! Everything uses real(c_kind_noahmp)/kind_noahmp so it stays precision-correct
! whether the library was built single or double precision.
! ===========================================================================
module NoahmpTestIOSupport

  use iso_c_binding,      only : C_INT
  use netcdf
  use Machine,            only : kind_noahmp, c_kind_noahmp
  use NoahmpIO_fi,        only : NoahmpIO_vect, NoahmpIOTypeVectInit_fi
  use NoahmpIOVarType,    only : NoahmpIO_type
  use NoahmpIOVarInitMod, only : NoahmpIOVarInitDefault

  implicit none
  private

  ! block/MPI setup + fixture generation
  public :: io_setup_block, make_wrfinput, tio_mpi_init, tio_mpi_finalize, ncchk
  ! deterministic reference field values (shared by generator and assertions)
  public :: wrf_xlat, wrf_xlong, wrf_terrain, wrf_tmn, wrf_tsk, wrf_canwat, &
            wrf_tslb, wrf_smois, wrf_ivgtyp, wrf_isltyp
  ! minimal assertion harness (module-global failure counter)
  public :: tio_reset, tio_expect, tio_expect_eq_i, tio_expect_close, &
            tio_nfail, tio_finish

  integer, save :: nfail = 0

  ! Grid dimensions built by WEST-EAST/SOUTH-NORTH_GRID_DIMENSION-1 convention.
  real(kind=kind_noahmp), parameter :: FILL_DX = 1000.0_kind_noahmp

contains

  ! ---- deterministic reference values (i,j are 0-based grid indices) --------
  ! Layer-independent soil values so the reader's vertical interpolation (file
  ! layers == model layers here) returns them unchanged, keeping assertions crisp.

  pure function wrf_xlat(i, j) result(v)
    integer, intent(in) :: i, j
    real(kind=c_kind_noahmp) :: v
    v = 35.0_c_kind_noahmp + 0.10_c_kind_noahmp*i + 0.01_c_kind_noahmp*j
  end function wrf_xlat

  pure function wrf_xlong(i, j) result(v)
    integer, intent(in) :: i, j
    real(kind=c_kind_noahmp) :: v
    v = -98.0_c_kind_noahmp + 0.10_c_kind_noahmp*i + 0.01_c_kind_noahmp*j
  end function wrf_xlong

  pure function wrf_terrain(i, j) result(v)
    integer, intent(in) :: i, j
    real(kind=c_kind_noahmp) :: v
    v = 100.0_c_kind_noahmp + real(i, c_kind_noahmp) + 10.0_c_kind_noahmp*j
  end function wrf_terrain

  pure function wrf_tmn(i, j) result(v)
    integer, intent(in) :: i, j
    real(kind=c_kind_noahmp) :: v
    v = 285.0_c_kind_noahmp + 0.5_c_kind_noahmp*i + 0.1_c_kind_noahmp*j
  end function wrf_tmn

  pure function wrf_tsk(i, j) result(v)
    integer, intent(in) :: i, j
    real(kind=c_kind_noahmp) :: v
    v = 290.0_c_kind_noahmp + real(i, c_kind_noahmp) + real(j, c_kind_noahmp)
  end function wrf_tsk

  pure function wrf_canwat(i, j) result(v)
    integer, intent(in) :: i, j
    real(kind=c_kind_noahmp) :: v
    v = 0.10_c_kind_noahmp * i
  end function wrf_canwat

  pure function wrf_tslb(i, j) result(v)
    integer, intent(in) :: i, j
    real(kind=c_kind_noahmp) :: v
    v = 280.0_c_kind_noahmp + real(i, c_kind_noahmp) + 10.0_c_kind_noahmp*j
  end function wrf_tslb

  pure function wrf_smois(i, j) result(v)
    integer, intent(in) :: i, j
    real(kind=c_kind_noahmp) :: v
    v = 0.20_c_kind_noahmp + 0.01_c_kind_noahmp*i
  end function wrf_smois

  pure function wrf_ivgtyp(i, j) result(v)
    integer, intent(in) :: i, j
    integer :: v
    v = 1 + mod(i + j, 10)
  end function wrf_ivgtyp

  pure function wrf_isltyp(i, j) result(v)
    integer, intent(in) :: i, j
    integer :: v
    v = 1 + mod(i + 2*j, 8)
  end function wrf_isltyp

  ! ---- assertion harness ---------------------------------------------------

  subroutine tio_reset()
    nfail = 0
  end subroutine tio_reset

  integer function tio_nfail()
    tio_nfail = nfail
  end function tio_nfail

  subroutine tio_expect(cond, name)
    logical,      intent(in) :: cond
    character(*), intent(in) :: name
    if (.not. cond) then
       write(0,'(A,A)') "  FAIL: ", name
       nfail = nfail + 1
    end if
  end subroutine tio_expect

  subroutine tio_expect_eq_i(got, want, name)
    integer,      intent(in) :: got, want
    character(*), intent(in) :: name
    if (got /= want) then
       write(0,'(A,A,A,I0,A,I0)') "  FAIL: ", name, " got=", got, " want=", want
       nfail = nfail + 1
    end if
  end subroutine tio_expect_eq_i

  subroutine tio_expect_close(got, want, name)
    real(kind=c_kind_noahmp), intent(in) :: got, want
    character(*),             intent(in) :: name
    real(kind=c_kind_noahmp) :: tol
    tol = 1.0e-3_c_kind_noahmp * max(1.0_c_kind_noahmp, abs(want))
    if (abs(got - want) > tol) then
       write(0,'(A,A,A,ES14.6,A,ES14.6)') "  FAIL: ", name, " got=", got, " want=", want
       nfail = nfail + 1
    end if
  end subroutine tio_expect_close

  ! Print a one-line summary and terminate with the CTest exit convention.
  subroutine tio_finish(name)
    character(*), intent(in) :: name
    if (nfail == 0) then
       write(*,'(A,A)') "PASS: ", name
       call exit(0)
    else
       write(*,'(A,A,A,I0,A)') "FAILED: ", name, " (", nfail, " check(s))"
       call exit(1)
    end if
  end subroutine tio_finish

  ! ---- MPI world (the test program owns it) --------------------------------

  ! Initialize MPI if not already up and return the Fortran MPI_COMM_WORLD handle
  ! to hand to NoahmpIO%comm for the collective parallel-NetCDF writers.
  subroutine tio_mpi_init(comm)
    use mpi
    integer, intent(out) :: comm
    integer :: ierr
    logical :: inited
    call MPI_Initialized(inited, ierr)
    if (.not. inited) call MPI_Init(ierr)
    comm = MPI_COMM_WORLD
  end subroutine tio_mpi_init

  subroutine tio_mpi_finalize()
    use mpi
    integer :: ierr
    logical :: finalized
    call MPI_Finalized(finalized, ierr)
    if (.not. finalized) call MPI_Finalize(ierr)
  end subroutine tio_mpi_finalize

  ! ---- NetCDF status check for the fixtures (abort loudly on a test bug) ----
  subroutine ncchk(status, context)
    integer,      intent(in) :: status
    character(*), intent(in) :: context
    if (status /= nf90_noerr) then
       write(0,'(A,A,A,A)') "test_io_support NetCDF error [", trim(context), "]: ", &
                            trim(nf90_strerror(status))
       error stop 3
    end if
  end subroutine ncchk

  ! ---- block setup ---------------------------------------------------------

  ! Size the module-global block array and populate one block's coupled scalars
  ! (normally wired from C++ by ScalarInitDefault) plus the ERF/grid geometry the
  ! I/O routines read, then allocate all storage via NoahmpIOVarInitDefault.
  ! A single block covers the whole nx*ny domain (xoffset=yoffset=0).
  subroutine io_setup_block(blk, level, blkid, nblocks, nx, ny, nsoil, nsnow, comm)
    type(NoahmpIO_type), pointer, intent(out) :: blk
    integer, intent(in) :: level, blkid, nblocks, nx, ny, nsoil, nsnow, comm

    integer(C_INT) :: lvl, nb
    integer :: k

    lvl = level
    nb  = nblocks
    call NoahmpIOTypeVectInit_fi(lvl, nb)
    blk => NoahmpIO_vect(level)%NoahmpIO(blkid)

    ! Coupled dimension/index scalars are C++-owned pointer components; in this
    ! pure-Fortran path give each its own target before use (cf. test_fortran_alloc).
    allocate(blk%XSTART);    blk%XSTART    = 0
    allocate(blk%XEND);      blk%XEND      = nx - 1
    allocate(blk%YSTART);    blk%YSTART    = 0
    allocate(blk%YEND);      blk%YEND      = ny - 1
    allocate(blk%KMS);       blk%KMS       = 1
    allocate(blk%KME);       blk%KME       = 3
    allocate(blk%NSOIL);     blk%NSOIL     = nsoil
    allocate(blk%NSNOW);     blk%NSNOW     = nsnow
    allocate(blk%NUMRAD);    blk%NUMRAD    = 2
    allocate(blk%ITIMESTEP)                       ! set by NoahmpIOVarInitDefault
    allocate(blk%RANK);      blk%RANK      = 0
    allocate(blk%BLKID);     blk%BLKID     = blkid
    allocate(blk%LEVEL);     blk%LEVEL     = level
    allocate(blk%COMM);      blk%COMM      = comm
    allocate(blk%DTBL);      blk%DTBL      = 3600.0_c_kind_noahmp

    ! Option branch scalars that gate allocations (stand-in for ReadNamelist).
    blk%IOPT_SOIL        = 1
    blk%IOPT_ALB         = 1
    blk%IOPT_WETLAND     = 0
    blk%SF_URBAN_PHYSICS = 0

    ! Global grid geometry (normally from ReadLandHeader); a single block spans it.
    blk%xsglobal = nx
    blk%ysglobal = ny
    blk%xoffset  = 0
    blk%yoffset  = 0

    ! Soil-layer thicknesses used to build DZS and the restart NSOIL geometry.
    if (allocated(blk%soil_thick_input)) deallocate(blk%soil_thick_input)
    allocate(blk%soil_thick_input(nsoil))
    if (nsoil == 4) then
       blk%soil_thick_input = [0.10_kind_noahmp, 0.30_kind_noahmp, &
                               0.60_kind_noahmp, 1.00_kind_noahmp]
    else
       do k = 1, nsoil
          blk%soil_thick_input(k) = 1.0_kind_noahmp / real(nsoil, kind_noahmp)
       end do
    end if

    call NoahmpIOVarInitDefault(blk)
  end subroutine io_setup_block

  ! ---- synthetic wrfinput / WPS fixture ------------------------------------

  ! Write a minimal but complete wrfinput-style NetCDF file for an nx*ny domain
  ! with nsoil soil layers, carrying every global attribute and variable that
  ! NoahmpReadLandHeader / NoahmpReadLandMain read. Serial NetCDF (no MPI). The
  ! WEST-EAST/SOUTH-NORTH_GRID_DIMENSION are nx+1/ny+1 (staggered) so the reader
  ! derives xsglobal=nx, ysglobal=ny. Variables carry the deterministic wrf_*
  ! reference values so the read test can assert exact recovery.
  subroutine make_wrfinput(path, nx, ny, nsoil)
    character(*), intent(in) :: path
    integer,      intent(in) :: nx, ny, nsoil

    integer :: ncid, d_time, d_we, d_sn, d_soil
    integer :: i, j, k
    integer, dimension(3) :: dims2d       ! (we, sn, Time)
    integer, dimension(4) :: dims3d       ! (we, sn, soil, Time) -- wrfinput soil order
    integer :: v_xlat, v_xlong, v_xland, v_seaice, v_hgt, v_tmn, v_mfx, v_mfy, &
               v_ivg, v_isl, v_canwat, v_tsk, v_snow, v_snowc, v_snowh, &
               v_dzs, v_tslb, v_smois, v_vegfra, v_lai, v_shdmin, v_shdmax
    real(kind=kind_noahmp), allocatable :: a2(:,:,:)      ! (we, sn, Time)
    real(kind=kind_noahmp), allocatable :: a3(:,:,:,:)    ! (we, sn, soil, Time)
    integer,                allocatable :: i2(:,:,:)      ! (we, sn, Time)
    real(kind=kind_noahmp), allocatable :: dzs(:,:)       ! (soil, Time)

    call ncchk(nf90_create(trim(path), IOR(NF90_CLOBBER, NF90_NETCDF4), ncid), "create wrfinput")

    call ncchk(nf90_def_dim(ncid, "Time",             1,     d_time), "def Time")
    call ncchk(nf90_def_dim(ncid, "west_east",        nx,    d_we),   "def west_east")
    call ncchk(nf90_def_dim(ncid, "south_north",      ny,    d_sn),   "def south_north")
    call ncchk(nf90_def_dim(ncid, "soil_layers_stag", nsoil, d_soil), "def soil")

    ! Global attributes read by NoahmpReadLandHeader (grid-dimension attrs are the
    ! staggered counts; the reader subtracts one).
    call ncchk(nf90_put_att(ncid, NF90_GLOBAL, "WEST-EAST_GRID_DIMENSION",   nx+1), "att WE")
    call ncchk(nf90_put_att(ncid, NF90_GLOBAL, "SOUTH-NORTH_GRID_DIMENSION", ny+1), "att SN")
    call ncchk(nf90_put_att(ncid, NF90_GLOBAL, "DX", FILL_DX), "att DX")
    call ncchk(nf90_put_att(ncid, NF90_GLOBAL, "DY", FILL_DX), "att DY")
    call ncchk(nf90_put_att(ncid, NF90_GLOBAL, "TRUELAT1", 30.0_kind_noahmp), "att TRUELAT1")
    call ncchk(nf90_put_att(ncid, NF90_GLOBAL, "TRUELAT2", 60.0_kind_noahmp), "att TRUELAT2")
    call ncchk(nf90_put_att(ncid, NF90_GLOBAL, "STAND_LON", -98.0_kind_noahmp), "att STAND_LON")
    call ncchk(nf90_put_att(ncid, NF90_GLOBAL, "MAP_PROJ", 1), "att MAP_PROJ")
    call ncchk(nf90_put_att(ncid, NF90_GLOBAL, "GRID_ID", 1), "att GRID_ID")
    call ncchk(nf90_put_att(ncid, NF90_GLOBAL, "ISWATER", 17), "att ISWATER")
    call ncchk(nf90_put_att(ncid, NF90_GLOBAL, "ISLAKE", -1), "att ISLAKE")
    call ncchk(nf90_put_att(ncid, NF90_GLOBAL, "ISURBAN", 13), "att ISURBAN")
    call ncchk(nf90_put_att(ncid, NF90_GLOBAL, "ISICE", 15), "att ISICE")
    call ncchk(nf90_put_att(ncid, NF90_GLOBAL, "MMINLU", "MODIFIED_IGBP_MODIS_NOAH"), "att MMINLU")
    ! Read from the level-0 setup file (same file here) to build x/y offsets.
    call ncchk(nf90_put_att(ncid, NF90_GLOBAL, "I_PARENT_START", 1), "att IPS")
    call ncchk(nf90_put_att(ncid, NF90_GLOBAL, "J_PARENT_START", 1), "att JPS")
    call ncchk(nf90_put_att(ncid, NF90_GLOBAL, "PARENT_GRID_RATIO", 1), "att PGR")

    dims2d = [d_we, d_sn, d_time]
    dims3d = [d_we, d_sn, d_soil, d_time]

    call def2d("XLAT",   v_xlat)
    call def2d("XLONG",  v_xlong)
    call def2d("XLAND",  v_xland)
    call def2d("SEAICE", v_seaice)
    call def2d("HGT",    v_hgt)
    call def2d("TMN",    v_tmn)
    call def2d("MAPFAC_MX", v_mfx)
    call def2d("MAPFAC_MY", v_mfy)
    call ncchk(nf90_def_var(ncid, "IVGTYP", NF90_INT, dims2d, v_ivg), "def IVGTYP")
    call ncchk(nf90_def_var(ncid, "ISLTYP", NF90_INT, dims2d, v_isl), "def ISLTYP")
    call def2d("CANWAT", v_canwat)
    call def2d("TSK",    v_tsk)
    call def2d("SNOW",   v_snow)
    call def2d("SNOWC",  v_snowc)
    call def2d("SNOWH",  v_snowh)
    call ncchk(nf90_def_var(ncid, "DZS", NF90_FLOAT, [d_soil, d_time], v_dzs), "def DZS")
    call ncchk(nf90_def_var(ncid, "TSLB",  NF90_FLOAT, dims3d, v_tslb),  "def TSLB")
    call ncchk(nf90_def_var(ncid, "SMOIS", NF90_FLOAT, dims3d, v_smois), "def SMOIS")
    call def2d("VEGFRA", v_vegfra)
    call def2d("LAI",    v_lai)
    call def2d("SHDMIN", v_shdmin)
    call def2d("SHDMAX", v_shdmax)
    call ncchk(nf90_enddef(ncid), "enddef")

    ! ---- fill 2-D fields ----
    allocate(a2(nx, ny, 1), i2(nx, ny, 1))
    call put2d(v_xlat,   fill2d(wrf_xlat))
    call put2d(v_xlong,  fill2d(wrf_xlong))
    call put2d(v_hgt,    fill2d(wrf_terrain))
    call put2d(v_tmn,    fill2d(wrf_tmn))
    call put2d(v_tsk,    fill2d(wrf_tsk))
    call put2d(v_canwat, fill2d(wrf_canwat))
    ! land everywhere, no sea-ice, unit map factors, benign veg/snow.
    call put2d(v_xland,  const2d(1.0_kind_noahmp))
    call put2d(v_seaice, const2d(0.0_kind_noahmp))
    call put2d(v_mfx,    const2d(1.0_kind_noahmp))
    call put2d(v_mfy,    const2d(1.0_kind_noahmp))
    call put2d(v_snow,   const2d(0.0_kind_noahmp))
    call put2d(v_snowc,  const2d(0.0_kind_noahmp))
    call put2d(v_snowh,  const2d(0.0_kind_noahmp))
    call put2d(v_vegfra, const2d(50.0_kind_noahmp))
    call put2d(v_lai,    const2d(2.0_kind_noahmp))
    call put2d(v_shdmin, const2d(10.0_kind_noahmp))
    call put2d(v_shdmax, const2d(80.0_kind_noahmp))

    do j = 1, ny
       do i = 1, nx
          i2(i, j, 1) = wrf_ivgtyp(i-1, j-1)
       end do
    end do
    call ncchk(nf90_put_var(ncid, v_ivg, i2), "put IVGTYP")
    do j = 1, ny
       do i = 1, nx
          i2(i, j, 1) = wrf_isltyp(i-1, j-1)
       end do
    end do
    call ncchk(nf90_put_var(ncid, v_isl, i2), "put ISLTYP")

    ! ---- soil layers (layer-independent so interpolation is a no-op) ----
    allocate(dzs(nsoil, 1))
    dzs(:, 1) = real(soil_thick(nsoil), kind_noahmp)
    call ncchk(nf90_put_var(ncid, v_dzs, dzs), "put DZS")

    ! On-disk order is (west_east, south_north, soil, Time), matching how
    ! get_netcdf_soillevel reads the hyperslab; values are layer-independent.
    allocate(a3(nx, ny, nsoil, 1))
    do k = 1, nsoil
       do j = 1, ny
          do i = 1, nx
             a3(i, j, k, 1) = real(wrf_tslb(i-1, j-1), kind_noahmp)
          end do
       end do
    end do
    call ncchk(nf90_put_var(ncid, v_tslb, a3), "put TSLB")
    do k = 1, nsoil
       do j = 1, ny
          do i = 1, nx
             a3(i, j, k, 1) = real(wrf_smois(i-1, j-1), kind_noahmp)
          end do
       end do
    end do
    call ncchk(nf90_put_var(ncid, v_smois, a3), "put SMOIS")

    call ncchk(nf90_close(ncid), "close wrfinput")

  contains

    subroutine def2d(name, vid)
      character(*), intent(in)  :: name
      integer,      intent(out) :: vid
      call ncchk(nf90_def_var(ncid, name, NF90_FLOAT, dims2d, vid), "def "//name)
    end subroutine def2d

    subroutine put2d(vid, arr)
      integer,                intent(in) :: vid
      real(kind=kind_noahmp), intent(in) :: arr(:,:,:)
      call ncchk(nf90_put_var(ncid, vid, arr), "put var")
    end subroutine put2d

    function fill2d(f) result(arr)
      interface
         pure function f(ii, jj) result(vv)
           import :: c_kind_noahmp
           integer, intent(in) :: ii, jj
           real(kind=c_kind_noahmp) :: vv
         end function f
      end interface
      real(kind=kind_noahmp) :: arr(nx, ny, 1)
      integer :: ii, jj
      do jj = 1, ny
         do ii = 1, nx
            arr(ii, jj, 1) = real(f(ii-1, jj-1), kind_noahmp)
         end do
      end do
    end function fill2d

    function const2d(c) result(arr)
      real(kind=kind_noahmp), intent(in) :: c
      real(kind=kind_noahmp) :: arr(nx, ny, 1)
      arr = c
    end function const2d

  end subroutine make_wrfinput

  ! Soil-layer thicknesses matching io_setup_block (kept in sync).
  pure function soil_thick(nsoil) result(t)
    integer, intent(in) :: nsoil
    real(kind=kind_noahmp) :: t(nsoil)
    integer :: k
    if (nsoil == 4) then
       t = [0.10_kind_noahmp, 0.30_kind_noahmp, 0.60_kind_noahmp, 1.00_kind_noahmp]
    else
       do k = 1, nsoil
          t(k) = 1.0_kind_noahmp / real(nsoil, kind_noahmp)
       end do
    end if
  end function soil_thick

  ! =========================================================================
  ! bind(C) shims for the C++ driver test (test_io_driver_cpp.cpp). They let the
  ! C++ side own the MPI world and reuse the Fortran fixture generator without
  ! pulling in MPI/NetCDF headers, and set the few Fortran-only NoahmpIO fields
  ! (erf_setup_file_*, YR/JULIAN) that the C++ mirror does not expose. Accessible
  ! from C via their bind(C) names regardless of module privacy.
  ! =========================================================================

  ! Bring up MPI (idempotent) and return the Fortran MPI_COMM_WORLD handle to
  ! store directly in NoahmpIO%comm (already a Fortran handle -- no c2f needed).
  function noahmp_test_mpi_init_c() result(comm) bind(C, name="noahmp_test_mpi_init_c")
    use iso_c_binding, only : C_INT
    integer(C_INT) :: comm
    integer :: fcomm
    call tio_mpi_init(fcomm)
    comm = fcomm
  end function noahmp_test_mpi_init_c

  subroutine noahmp_test_mpi_finalize_c() bind(C, name="noahmp_test_mpi_finalize_c")
    call tio_mpi_finalize()
  end subroutine noahmp_test_mpi_finalize_c

  ! Generate the synthetic wrfinput fixture at the given (C) path.
  subroutine noahmp_test_make_wrfinput_c(path, plen, nx, ny, nsoil) &
       bind(C, name="noahmp_test_make_wrfinput_c")
    use iso_c_binding, only : C_INT, C_CHAR
    character(kind=C_CHAR), intent(in) :: path(*)
    integer(C_INT), value, intent(in) :: plen, nx, ny, nsoil
    call make_wrfinput(cstr(path, plen), nx, ny, nsoil)
  end subroutine noahmp_test_make_wrfinput_c

  ! Point a block's setup-file fields at a path (Fortran-only NoahmpIO members).
  subroutine noahmp_test_set_setup_file_c(level, blkid, path, plen) &
       bind(C, name="noahmp_test_set_setup_file_c")
    use iso_c_binding, only : C_INT, C_CHAR
    integer(C_INT), value, intent(in) :: level, blkid, plen
    character(kind=C_CHAR), intent(in) :: path(*)
    character(len=plen) :: f
    f = cstr(path, plen)
    NoahmpIO_vect(level)%NoahmpIO(blkid)%erf_setup_file_lev = f
    NoahmpIO_vect(level)%NoahmpIO(blkid)%erf_setup_file_01  = f
  end subroutine noahmp_test_set_setup_file_c

  ! Set the Fortran-only time scalars (ERF will wire these from its clock; the
  ! test sets sane values so phenology/lengths are well-defined).
  subroutine noahmp_test_set_time_c(level, blkid, yr, julian) &
       bind(C, name="noahmp_test_set_time_c")
    use iso_c_binding, only : C_INT
    integer(C_INT), value, intent(in) :: level, blkid, yr
    real(kind=c_kind_noahmp), value, intent(in) :: julian
    NoahmpIO_vect(level)%NoahmpIO(blkid)%YR     = yr
    NoahmpIO_vect(level)%NoahmpIO(blkid)%JULIAN = julian
  end subroutine noahmp_test_set_time_c

  ! Read west_east/south_north/soil_layers_stag dims from a wrfinput/WPS file.
  subroutine noahmp_test_query_dims_c(path, plen, nx, ny, nsoil) &
       bind(C, name="noahmp_test_query_dims_c")
    use iso_c_binding, only : C_INT, C_CHAR
    integer(C_INT), value, intent(in)  :: plen
    integer(C_INT),        intent(out) :: nx, ny, nsoil
    character(kind=C_CHAR), intent(in) :: path(*)
    character(len=plen) :: f
    integer :: ncid, did
    f = cstr(path, plen)
    call ncchk(nf90_open(trim(f), NF90_NOWRITE, ncid), "open ext")
    call ncchk(nf90_inq_dimid(ncid, "west_east", did), "we");   call ncchk(nf90_inquire_dimension(ncid, did, len=nx), "we")
    call ncchk(nf90_inq_dimid(ncid, "south_north", did), "sn"); call ncchk(nf90_inquire_dimension(ncid, did, len=ny), "sn")
    call ncchk(nf90_inq_dimid(ncid, "soil_layers_stag", did), "soil"); call ncchk(nf90_inquire_dimension(ncid, did, len=nsoil), "soil")
    call ncchk(nf90_close(ncid), "close ext")
  end subroutine noahmp_test_query_dims_c

  ! Read a field element straight from the Fortran module-global block, so the C++
  ! driver test can prove its C++ view and the Fortran storage alias the SAME
  ! element at a given (i,j)/(i,k,j). The C++ side drives everything through its
  ! own views, so a transposed or offset index map would be self-consistent and
  ! invisible without a cross-language read like this. (2-D: HFX; 3-D column-major
  ! (i,layer,j): TSLB.)
  function noahmp_test_read_hfx_c(level, blkid, i, j) result(val) &
       bind(C, name="noahmp_test_read_hfx_c")
    use iso_c_binding, only : C_INT
    integer(C_INT), value, intent(in) :: level, blkid, i, j
    real(kind=c_kind_noahmp) :: val
    val = NoahmpIO_vect(level)%NoahmpIO(blkid)%HFX(i, j)
  end function noahmp_test_read_hfx_c

  function noahmp_test_read_tslb_c(level, blkid, i, k, j) result(val) &
       bind(C, name="noahmp_test_read_tslb_c")
    use iso_c_binding, only : C_INT
    integer(C_INT), value, intent(in) :: level, blkid, i, k, j
    real(kind=c_kind_noahmp) :: val
    val = NoahmpIO_vect(level)%NoahmpIO(blkid)%TSLB(i, k, j)
  end function noahmp_test_read_tslb_c

  ! Convert a C character buffer of known length to a Fortran string.
  pure function cstr(path, plen) result(f)
    use iso_c_binding, only : C_CHAR
    character(kind=C_CHAR), intent(in) :: path(*)
    integer,               intent(in) :: plen
    character(len=plen) :: f
    integer :: k
    do k = 1, plen
       f(k:k) = path(k)
    end do
  end function cstr

end module NoahmpTestIOSupport
