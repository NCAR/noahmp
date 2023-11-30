module Machine

!!! define machine-related constants and parameters
!!! To define real data type precision, use "-DOUBLE_PREC" in CPPFLAG in user_build_options file
!!! By default, Noah-MP uses single precision

! ------------------------ Code history -----------------------------------
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! -------------------------------------------------------------------------

  implicit none
  save
  private

#ifdef DOUBLE_PREC
  integer, public, parameter :: kind_noahmp = 8 ! double precision
#else
  integer, public, parameter :: kind_noahmp = 4 ! single precision
#endif

  integer,                public, parameter :: undefined_int  = -99999       ! undefined integer for variable initialization
  real(kind=kind_noahmp), public, parameter :: undefined_real = -99999.0     ! undefined real for variable initializatin
  integer,                public, parameter :: undefined_int_neg  = -99999   ! undefined integer negative for variable initialization
  real(kind=kind_noahmp), public, parameter :: undefined_real_neg = -99999.0 ! undefined real negative for variable initializatin

end module Machine
