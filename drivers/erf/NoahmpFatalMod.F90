module NoahmpFatalMod

! Centralized fatal-exit and NetCDF error checking for the ERF Noah-MP driver.
! Leaf module (depends only on netcdf) so every module can `use` it without a
! cycle. NoahmpIO_abort routes through the C shim NoahmpIO_fatal_c (NoahmpIO.cpp),
! whose host-installed handler (e.g. amrex::Abort) propagates the failure across
! all ranks -- keeping Noah-MP free of a direct MPI/AMReX dependency (see
! NoahmpFatal.H).

   use netcdf, only : nf90_noerr, nf90_strerror

   implicit none
   private

   public :: NoahmpIO_abort, check_nc

contains

   ! Terminate after a fatal error. Routes through the host-installed handler for
   ! cross-rank propagation; the trailing `error stop` is the serial fallback.
   subroutine NoahmpIO_abort()

      implicit none

      interface
         subroutine NoahmpIO_fatal_c() bind(C, name="NoahmpIO_fatal_c")
         end subroutine NoahmpIO_fatal_c
      end interface

      flush(0)
      call NoahmpIO_fatal_c()
      error stop 1

   end subroutine NoahmpIO_abort

   ! Abort with context if a NetCDF call did not succeed.
   subroutine check_nc(status, context)

      implicit none

      integer,          intent(in) :: status   ! return code from an nf90_* call
      character(len=*), intent(in) :: context  ! variable/operation name for diagnostics

      if (status /= nf90_noerr) then
         print *, "NetCDF error [", trim(context), "]: ", trim(nf90_strerror(status))
         call NoahmpIO_abort()
      end if

   end subroutine check_nc

end module NoahmpFatalMod
