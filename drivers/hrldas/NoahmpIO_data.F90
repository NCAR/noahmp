module NoahmpIO_data

  use NoahmpIOVarType, ONLY: NoahmpIO_type

  implicit none

  ! ---------------------------------------------------------------------------
  ! Public variable of NoahmpIO_type
  ! --------------------------------------------------------------------------- 
  type(NoahmpIO_type), save, target, public :: NoahmpIO

end module NoahmpIO_data
