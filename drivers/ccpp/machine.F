      module machine

!! \section arg_table_machine
!! \htmlinclude machine.html
!!

      implicit none

      integer, parameter :: kind_sngl_prec = 4                          &
     &,                     kind_dbl_prec = 8                           &
# ifdef __PGI
     &,                     kind_qdt_prec = 8                           &
# else
     &,                     kind_qdt_prec = 16                          &
# endif
     &,                     kind_integer = 4                            &
     &,                     kind_logical = 4
     &,                     kind_io4  = kind_sngl_prec                  &
     &,                     kind_ior = kind_dbl_prec                    &
     &,                     kind_grid = kind_dbl_prec

! Physics single precision flag
#ifndef SINGLE_PREC
      integer, parameter :: kind_phys = kind_dbl_prec
#else
      integer, parameter :: kind_phys = kind_sngl_prec
#endif

      integer, parameter :: kind_io8 = kind_phys                        &! Note kind_io8 is not always 8 bytes

! Dynamics single precision flag
#ifdef OVERLOAD_R4
      integer, parameter :: kind_dyn  = kind_sngl_prec
#else
      integer, parameter :: kind_dyn  = kind_dbl_prec
#endif

!
      real(kind=kind_phys), parameter :: mprec = 1.e-12           ! machine precision to restrict dep
      real(kind=kind_phys), parameter :: grib_undef = 9.99e20     ! grib undefine value
!
      end module machine
