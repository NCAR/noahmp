! ===========================================================================
! test_io_readland -- NetCDF READ path (NoahmpReadLandMod).
!
! Drives the real NoahmpReadLandHeader + NoahmpReadLandMain against a
! wrfinput/WPS-style NetCDF file and asserts the header globals and per-cell
! field values land in the NoahmpIO arrays.
!
! Generic by construction: instead of the >200 MB local ChisholmView fixture,
! the test SYNTHESIZES a tiny self-contained wrfinput (make_wrfinput) with the
! deterministic wrf_* reference fields, so it runs anywhere with no external
! data and asserts EXACT recovery. To additionally exercise a real local file,
! set NOAHMP_TEST_WRFINPUT=/path/to/wrfinput (any WPS-type file with the
! standard west_east / south_north / soil_layers_stag dims); that mode sizes the
! block from the file's own dimensions and asserts sanity (extents + finite,
! in-range values) since the true values are unknown.
!
! The program owns the MPI world (tio_mpi_init); the read path itself is serial
! NetCDF, but keeping MPI up matches the coupled runtime and the other I/O tests.
! Exit code 0 = pass.
! ===========================================================================
program test_io_readland

  use netcdf
  use Machine,             only : c_kind_noahmp
  use NoahmpIOVarType,     only : NoahmpIO_type
  use NoahmpReadLandMod,   only : NoahmpReadLandHeader, NoahmpReadLandMain
  use NoahmpTestIOSupport

  implicit none

  type(NoahmpIO_type), pointer :: blk
  integer :: comm
  integer :: nx, ny, nsoil
  integer :: i, j, k
  logical :: synthetic
  character(len=1024) :: ext_path
  character(len=1024) :: path
  integer :: elen, estat

  call tio_reset()
  call tio_mpi_init(comm)

  ! Choose the fixture: real file from the environment, or a synthetic one.
  call get_environment_variable("NOAHMP_TEST_WRFINPUT", ext_path, elen, estat)
  synthetic = (estat /= 0 .or. elen == 0)

  if (synthetic) then
     nx = 4; ny = 3; nsoil = 4
     path = "test_wrfinput_synth.nc"
     call make_wrfinput(trim(path), nx, ny, nsoil)
  else
     path = trim(ext_path)
     call query_dims(trim(path), nx, ny, nsoil)
     write(*,'(A,A,A,I0,A,I0,A,I0)') "Using external wrfinput: ", trim(path), &
          "  nx=", nx, " ny=", ny, " nsoil=", nsoil
  end if

  ! Build a single block that spans the whole nx*ny domain, then point the reader
  ! at the chosen file (both the current level and the level-0 setup file).
  call io_setup_block(blk, level=0, blkid=0, nblocks=1, nx=nx, ny=ny, &
                      nsoil=nsoil, nsnow=3, comm=comm)
  blk%erf_setup_file_lev = trim(path)
  blk%erf_setup_file_01  = trim(path)

  call NoahmpReadLandHeader(blk)

  ! ---- header globals ----
  call tio_expect_eq_i(blk%xsglobal, nx, "xsglobal = WE_GRID_DIM-1")
  call tio_expect_eq_i(blk%ysglobal, ny, "ysglobal = SN_GRID_DIM-1")
  call tio_expect_eq_i(blk%xoffset,  0,  "xoffset (top level) = 0")
  call tio_expect_eq_i(blk%yoffset,  0,  "yoffset (top level) = 0")

  call NoahmpReadLandMain(blk)

  ! ReadLandMain sets these regardless of the file.
  call tio_expect_eq_i(blk%itimestep, 1, "itimestep reset to 1")
  call tio_expect(.not. blk%restart_flag, "restart_flag cleared")
  do k = 1, nsoil
     call tio_expect_close(real(blk%dzs(k), c_kind_noahmp), &
          real(blk%soil_thick_input(k), c_kind_noahmp), "dzs = soil_thick_input")
  end do

  if (synthetic) then
     ! ---- exact per-cell recovery (0-based grid indices) ----
     do j = 0, ny-1
        do i = 0, nx-1
           call tio_expect_close(blk%xlat(i,j),    wrf_xlat(i,j),    "xlat")
           call tio_expect_close(blk%xlong(i,j),   wrf_xlong(i,j),   "xlong")
           call tio_expect_close(blk%terrain(i,j), wrf_terrain(i,j), "terrain")
           call tio_expect_close(blk%tmn(i,j),     wrf_tmn(i,j),     "tmn")
           call tio_expect_close(blk%tsk(i,j),     wrf_tsk(i,j),     "tsk")
           call tio_expect_close(blk%canwat(i,j),  wrf_canwat(i,j),  "canwat")
           call tio_expect_eq_i(blk%ivgtyp(i,j),   wrf_ivgtyp(i,j),  "ivgtyp (nint)")
           call tio_expect_eq_i(blk%isltyp(i,j),   wrf_isltyp(i,j),  "isltyp (nint)")
           ! soil fields are layer-independent, so interpolation is a no-op.
           do k = 1, nsoil
              call tio_expect_close(blk%tslb(i,k,j),  wrf_tslb(i,j),  "tslb")
              call tio_expect_close(blk%smois(i,k,j), wrf_smois(i,j), "smois")
           end do
        end do
     end do
  else
     ! ---- sanity for an unknown real file ----
     call tio_expect(all(abs(blk%xlat)  <= 90.0_c_kind_noahmp),  "xlat in [-90,90]")
     call tio_expect(all(abs(blk%xlong) <= 360.0_c_kind_noahmp), "xlong in [-360,360]")
     call tio_expect(all(blk%ivgtyp >= 0), "ivgtyp non-negative")
     call tio_expect(all(blk%isltyp >= 0), "isltyp non-negative")
     call tio_expect(all(blk%tslb > 150.0_c_kind_noahmp .and. blk%tslb < 400.0_c_kind_noahmp), &
                     "tslb physically plausible [K]")
  end if

  call tio_mpi_finalize()
  call tio_finish("test_io_readland")

contains

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

end program test_io_readland
