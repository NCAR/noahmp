module Machine

!!! define machine-related constants and parameters
!!! To define real data type precision, use "-DOUBLE_PREC" in CPPFLAG in user_build_options file
!!! By default, Noah-MP uses single precision

! ------------------------ Code history -----------------------------------
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! -------------------------------------------------------------------------

  use, intrinsic :: iso_c_binding, only : C_DOUBLE, C_FLOAT

  implicit none
  save
  private

#ifdef DOUBLE_PREC
  integer, public, parameter :: kind_noahmp   = 8        ! double precision
  integer, public, parameter :: c_kind_noahmp = C_DOUBLE ! matching C kind
#else
  integer, public, parameter :: kind_noahmp   = 4        ! single precision
  integer, public, parameter :: c_kind_noahmp = C_FLOAT  ! matching C kind
#endif

  integer,                public, parameter :: undefined_int  = -9999       ! undefined integer for variable initialization
  real(kind=kind_noahmp), public, parameter :: undefined_real = -9999.0     ! undefined real for variable initializatin
  integer,                public, parameter :: undefined_int_neg  = -9999   ! undefined integer negative for variable initialization
  real(kind=kind_noahmp), public, parameter :: undefined_real_neg = -9999.0 ! undefined real negative for variable initializatin

end module Machine
