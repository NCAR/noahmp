! ===========================================================================
! test_fortran_guards -- input-validation / abort paths (Tier 1).
!
! NoahmpIOTypeVectInit_fi validates its arguments and aborts (write diagnostic +
! NoahmpIO_abort -> std::abort with no handler installed) rather than corrupting
! module state. Each abort is a separate CTest case registered WILL_FAIL, so a
! non-zero termination = the guard fired = the test passes. If a call instead
! returns normally the program exits 0, which WILL_FAIL flags as a failure --
! catching a guard that silently stopped firing.
!
! Which guard to trigger is selected by argv(1): level | nblocks | reinit.
! ===========================================================================
program test_fortran_guards

  use iso_c_binding, only : C_INT
  use NoahmpIO_fi,   only : NoahmpIOTypeVectInit_fi, NLEVEL_MAX

  implicit none
  character(len=32) :: which
  integer(C_INT) :: level, nb

  call get_command_argument(1, which)

  select case (trim(which))
  case ("level")
     ! level above NLEVEL_MAX -> out-of-range abort
     level = NLEVEL_MAX + 1
     nb    = 1
     call NoahmpIOTypeVectInit_fi(level, nb)
  case ("nblocks")
     ! NBlocks < 1 -> abort
     level = 0
     nb    = 0
     call NoahmpIOTypeVectInit_fi(level, nb)
  case ("reinit")
     ! second init with a DIFFERENT size -> refuse/abort
     level = 0
     call NoahmpIOTypeVectInit_fi(level, 1_C_INT)
     call NoahmpIOTypeVectInit_fi(level, 2_C_INT)
  case default
     write(0,'(A)') "test_fortran_guards: unknown selector '"//trim(which)//"'"
     call exit(2)
  end select

  ! Should be unreachable: the guard above was expected to abort.
  write(0,'(A)') "test_fortran_guards: expected abort did NOT occur for '"// &
                 trim(which)//"'"
  call exit(0)

end program test_fortran_guards
