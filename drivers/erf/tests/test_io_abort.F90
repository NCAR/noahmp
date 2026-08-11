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
!   checkbad -- check_nc(<a NetCDF error>, ...) must abort. Detected by output:
!               PASS on the forced-error context, FAIL if it reaches "did NOT
!               abort" (check_nc prints its error line BEFORE it decides to
!               abort, so an exit-code test alone could not tell the two apart).
!   abort    -- NoahmpIO_abort() called directly must not return. Detected by
!               output: PASS on the pre-call marker, FAIL on "returned
!               unexpectedly" (the call itself emits no diagnostic).
!
! The abort cases are pinned by PASS/FAIL_REGULAR_EXPRESSION in tests/CMakeLists.txt
! rather than WILL_FAIL -- WILL_FAIL accepts ANY non-zero exit (a wrong-reason
! abort, a mistyped selector) as a spurious pass. test_abort_handler.cpp is still
! linked so the abort terminates cleanly with its output flushed; the "ok" case
! never triggers the handler.
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
     ! NoahmpIO_abort must terminate; it is [[noreturn]] on the C side. It emits
     ! no diagnostic of its own, so print a marker BEFORE the call: the test is
     ! detected by output (PASS on this marker, FAIL on the "returned
     ! unexpectedly" line below), not by exit code -- see tests/CMakeLists.txt.
     write(*,'(A)') "test_io_abort: reached direct-abort call"
     flush(6)
     call NoahmpIO_abort()
     write(0,'(A)') "test_io_abort: NoahmpIO_abort returned unexpectedly"
     call exit(0)

  case default
     write(0,'(A)') "test_io_abort: unknown selector '"//trim(which)//"'"
     call exit(2)
  end select

end program test_io_abort
