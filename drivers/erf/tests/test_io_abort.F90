! ===========================================================================
! test_io_abort -- NoahmpFatalMod fatal-exit + NetCDF error checking.
!
! NoahmpFatalMod is the centralized abort path shared by every NetCDF I/O
! module: check_nc(status, ctx) aborts via NoahmpIO_abort() on any non-success
! status, and NoahmpIO_abort() routes through the host handler (C shim
! NoahmpIO_fatal_c) so a failure propagates across ranks, with `error stop` as
! the serial fallback.
!
! Selector argv(1):
!   ok       -- check_nc(NF90_NOERR, ...) must RETURN normally (no abort); the
!               program then exits 0. Registered as a normal must-pass case.
!   checkbad -- check_nc(<a NetCDF error>, ...) must abort. WILL_FAIL.
!   abort    -- NoahmpIO_abort() called directly must not return. WILL_FAIL.
!
! Built with test_abort_handler.cpp linked in, which installs a fatal handler
! that terminates via std::_Exit(7) instead of SIGABRT, so CTest's WILL_FAIL
! inverts the abort cases cleanly (a signal death is not inverted). The "ok"
! case never triggers the handler.
! ===========================================================================
program test_io_abort

  use netcdf,         only : nf90_noerr, nf90_ebadid
  use NoahmpFatalMod, only : check_nc, NoahmpIO_abort

  implicit none

  character(len=32) :: which

  call get_command_argument(1, which)

  select case (trim(which))
  case ("ok")
     ! A success status must fall through without aborting.
     call check_nc(nf90_noerr, "test_io_abort: success path")
     write(*,'(A)') "PASS: test_io_abort (check_nc returned on NF90_NOERR)"
     call exit(0)

  case ("checkbad")
     ! A genuine NetCDF error status must abort inside check_nc.
     call check_nc(nf90_ebadid, "test_io_abort: forced NetCDF error")
     ! Unreachable if check_nc aborted as required.
     write(0,'(A)') "test_io_abort: check_nc did NOT abort on an error status"
     call exit(0)

  case ("abort")
     ! NoahmpIO_abort must terminate; it is [[noreturn]] on the C side.
     call NoahmpIO_abort()
     write(0,'(A)') "test_io_abort: NoahmpIO_abort returned unexpectedly"
     call exit(0)

  case default
     write(0,'(A)') "test_io_abort: unknown selector '"//trim(which)//"'"
     call exit(2)
  end select

end program test_io_abort
