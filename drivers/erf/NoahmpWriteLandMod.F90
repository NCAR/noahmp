module NoahmpWriteLandMod

!!!  Module to initialize Noah-MP 2-D variables

   use mpi
   use netcdf
   use NoahmpIOVarType

   implicit none

contains

   subroutine NoahmpWriteLand(NoahmpIO, filenum)

      implicit none

      type(NoahmpIO_type), intent(inout) :: NoahmpIO
      integer, intent(in) :: filenum

      ! local variables
      integer :: ierr, ncid, nx, ny, terrain, snowh, start(2), count(2)
      integer :: shbxy, evbxy
      character(len=5) :: ts_str
      character(len=1) :: lev_str
      character(len=100) :: dir, filename
      logical :: ex

      write (ts_str, '(I5.5)') filenum
      write (lev_str, '(I1.1)') NoahmpIO%LEVEL

      dir = "lnd"//trim(ts_str)
      inquire (file=trim(dir), exist=ex)
      if (.not. ex) then
         call execute_command_line("mkdir -p "//trim(dir), exitstat=ierr)
         if (ierr /= 0) then
            print *, "Failed to create directory: ", trim(dir)
            stop
         end if
      end if

      filename = trim(dir)//"/Level_"//trim(lev_str)//".nc"
      ierr = nf90_create(trim(filename), IOR(NF90_CLOBBER, IOR(NF90_NETCDF4, NF90_MPIIO)), &
                         ncid, comm=NoahmpIO%comm, info=MPI_INFO_NULL)

      if (ierr /= nf90_noerr) then
         print *, "NetCDF create failed: ", trim(nf90_strerror(ierr))
         stop
      end if

      ! Define dimensions and variable
      ierr = nf90_def_dim(ncid, "NX", NoahmpIO%xsglobal, nx)
      ierr = nf90_def_dim(ncid, "NY", NoahmpIO%ysglobal, ny)
      ierr = nf90_def_var(ncid, "TERRAIN", NF90_FLOAT, (/nx, ny/), terrain)
      ierr = nf90_def_var(ncid, "SNOWH", NF90_FLOAT, (/nx, ny/), snowh)
      ierr = nf90_def_var(ncid, "SHBXY", NF90_FLOAT, (/nx, ny/), shbxy)
      ierr = nf90_def_var(ncid, "EVBXY", NF90_FLOAT, (/nx, ny/), evbxy)

      ! End definition mode
      ierr = nf90_enddef(ncid)

      ! Select portion to write
      start = (/NoahmpIO%xstart-NoahmpIO%xoffset+1, NoahmpIO%ystart-NoahmpIO%yoffset+1/)
      count = (/NoahmpIO%xend-NoahmpIO%xstart+1, NoahmpIO%yend-NoahmpIO%ystart+1/)

      ! Write data
      ierr = nf90_put_var(ncid, terrain, NoahmpIO%TERRAIN, start=start, count=count)
      ierr = nf90_put_var(ncid, snowh, NoahmpIO%SNOWH, start=start, count=count)
      ierr = nf90_put_var(ncid, shbxy, NoahmpIO%SHBXY, start=start, count=count)
      ierr = nf90_put_var(ncid, evbxy, NoahmpIO%EVBXY, start=start, count=count)

      ! Close file
      ierr = nf90_close(ncid)

   end subroutine NoahmpWriteLand

end module NoahmpWriteLandMod
