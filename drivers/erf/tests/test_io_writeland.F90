! ===========================================================================
! test_io_writeland -- NetCDF land-output WRITE path (NoahmpWriteLandMod).
!
! NoahmpWriteLand serializes the diagnostic/coupling output subset to
! lnd<filenum>/Level_<L>.nc using the same collective parallel NetCDF-4 pattern
! as the restart writer (NF90_MPIIO, comm=NoahmpIO%comm), but at NF90_FLOAT
! (lossy output precision). The program creates its own MPI world and runs
! single-rank, so one block's hyperslab covers the whole domain.
!
! The test fills a few representative output fields, calls WriteLand, then
! REOPENS the produced file serially and asserts the dimensions and the written
! values round-trip (within float tolerance). Exit code 0 = pass.
! ===========================================================================
program test_io_writeland

  use netcdf
  use Machine,             only : kind_noahmp, c_kind_noahmp
  use NoahmpIOVarType,     only : NoahmpIO_type
  use NoahmpWriteLandMod,  only : NoahmpWriteLand
  use NoahmpTestIOSupport

  implicit none

  type(NoahmpIO_type), pointer :: blk
  integer :: comm
  integer, parameter :: NX = 4, NY = 3, NSOIL = 4, NSNOW = 3
  integer, parameter :: FILENUM = 1
  integer :: i, j, k
  character(len=64) :: fname

  call tio_reset()
  call tio_mpi_init(comm)
  call io_setup_block(blk, level=0, blkid=0, nblocks=1, nx=NX, ny=NY, &
                      nsoil=NSOIL, nsnow=NSNOW, comm=comm)

  ! Fill representative output fields with deterministic values (0-based grid).
  do j = 0, NY-1
     do i = 0, NX-1
        blk%TERRAIN(i,j) = wrf_terrain(i,j)
        blk%HFX(i,j)     = 10.0_c_kind_noahmp + i + j
        blk%TSK(i,j)     = wrf_tsk(i,j)
        do k = 1, NSOIL
           blk%TSLB(i,k,j) = 280.0_c_kind_noahmp + i + k + 10.0_c_kind_noahmp*j
        end do
     end do
  end do

  call NoahmpWriteLand(blk, FILENUM, 1)

  ! WriteLand builds "lnd" + zero-padded(>=5) filenum + "/Level_<L>.nc".
  fname = "lnd00001/Level_0.nc"
  call check_output(trim(fname))

  call tio_mpi_finalize()
  call tio_finish("test_io_writeland")

contains

  subroutine check_output(fpath)
    character(*), intent(in) :: fpath
    integer :: ncid, did, vid, dnx, dny, dns
    real(kind=kind_noahmp) :: r2(NX,NY), r3(NX,NSOIL,NY)

    call ncchk(nf90_open(fpath, NF90_NOWRITE, ncid), "reopen land output")

    ! ---- dimensions ----
    call ncchk(nf90_inq_dimid(ncid, "NX", did), "dimid NX")
    call ncchk(nf90_inquire_dimension(ncid, did, len=dnx), "len NX")
    call ncchk(nf90_inq_dimid(ncid, "NY", did), "dimid NY")
    call ncchk(nf90_inquire_dimension(ncid, did, len=dny), "len NY")
    call ncchk(nf90_inq_dimid(ncid, "NSOIL", did), "dimid NSOIL")
    call ncchk(nf90_inquire_dimension(ncid, did, len=dns), "len NSOIL")
    call tio_expect_eq_i(dnx, blk%xsglobal, "output NX = xsglobal")
    call tio_expect_eq_i(dny, blk%ysglobal, "output NY = ysglobal")
    call tio_expect_eq_i(dns, NSOIL,        "output NSOIL")

    ! ---- 2D field values (file element (i+1,j+1) == blk field at (i,j)) ----
    call ncchk(nf90_inq_varid(ncid, "TERRAIN", vid), "inq TERRAIN")
    call ncchk(nf90_get_var(ncid, vid, r2), "get TERRAIN")
    do j = 0, NY-1
       do i = 0, NX-1
          call tio_expect_close(real(r2(i+1,j+1), c_kind_noahmp), wrf_terrain(i,j), "out TERRAIN")
       end do
    end do

    call ncchk(nf90_inq_varid(ncid, "TSK", vid), "inq TSK")
    call ncchk(nf90_get_var(ncid, vid, r2), "get TSK")
    do j = 0, NY-1
       do i = 0, NX-1
          call tio_expect_close(real(r2(i+1,j+1), c_kind_noahmp), wrf_tsk(i,j), "out TSK")
       end do
    end do

    call ncchk(nf90_inq_varid(ncid, "HFX", vid), "inq HFX")
    call ncchk(nf90_get_var(ncid, vid, r2), "get HFX")
    do j = 0, NY-1
       do i = 0, NX-1
          call tio_expect_close(real(r2(i+1,j+1), c_kind_noahmp), &
               10.0_c_kind_noahmp + i + j, "out HFX")
       end do
    end do

    ! ---- 3D soil field (NX, NSOIL, NY) ----
    call ncchk(nf90_inq_varid(ncid, "TSLB", vid), "inq TSLB")
    call ncchk(nf90_get_var(ncid, vid, r3), "get TSLB")
    do j = 0, NY-1
       do k = 1, NSOIL
          do i = 0, NX-1
             call tio_expect_close(real(r3(i+1,k,j+1), c_kind_noahmp), &
                  280.0_c_kind_noahmp + i + k + 10.0_c_kind_noahmp*j, "out TSLB")
          end do
       end do
    end do

    call ncchk(nf90_close(ncid), "close land output")
  end subroutine check_output

end program test_io_writeland
