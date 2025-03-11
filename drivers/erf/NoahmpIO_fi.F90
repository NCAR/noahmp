module NoahmpIO_fi

  use NoahmpIOVarType, ONLY: NoahmpIO_type
  use NoahmpIOVarInitMod, ONLY: NoahmpIOVarInitDefault
  use NoahmpInitMainMod, ONLY: NoahmpInitMain
  use NoahmpReadNamelistMod, ONLY: NoahmpReadNamelist
  use NoahmpReadTableMod, ONLY: NoahmpReadTable

  use, intrinsic :: iso_c_binding

  implicit none

  ! ---------------------------------------------------------------------------
  ! Public variable of NoahmpIO_type
  ! --------------------------------------------------------------------------- 
  type(NoahmpIO_type), save, target, public :: NoahmpIO

  ! ---------------------------------------------------------------------------
  ! Mirror of extern C struct
  ! ---------------------------------------------------------------------------
  type, bind(c), public :: NoahmpIO_type_fi
    type(C_PTR)                                         ::  ids,ide, &          ! d -> domain
                                                            jds,jde, &          ! d -> domain
                                                            kds,kde, &          ! d -> domain
                                                            ims,ime, &          ! m -> memory
                                                            jms,jme, &          ! m -> memory
                                                            kms,kme, &          ! m -> memory
                                                            its,ite, &          ! t -> tile
                                                            jts,jte, &          ! t -> tile
                                                            kts,kte             ! t -> tile

    type(C_PTR) :: xstart, xend, ystart, yend
    type(C_PTR) :: nsoil, nsnow
    type(C_PTR) :: XLAT, WSLAKEXY
  end type NoahmpIO_type_fi

contains

  subroutine NoahmpIOTypeInit_fi(NoahmpIO_cptr) bind(C, name="NoahmpIOTypeInit_fi")
    use, intrinsic :: iso_c_binding
    implicit none
    type(NoahmpIO_type_fi), intent(inout) :: NoahmpIO_cptr

    call C_F_POINTER(NoahmpIO_cptr%XSTART, NoahmpIO%XSTART)
    call C_F_POINTER(NoahmpIO_cptr%XEND,   NoahmpIO%XEND)
    call C_F_POINTER(NoahmpIO_cptr%YSTART, NoahmpIO%YSTART)
    call C_F_POINTER(NoahmpIO_cptr%YEND,   NoahmpIO%YEND)

    call C_F_POINTER(NoahmpIO_cptr%NSOIL,  NoahmpIO%NSOIL)
    call C_F_POINTER(NoahmpIO_cptr%NSNOW,  NoahmpIO%NSNOW)

    call C_F_POINTER(NoahmpIO_cptr%IDS,    NoahmpIO%IDS)
    call C_F_POINTER(NoahmpIO_cptr%IDE,    NoahmpIO%IDE)
    call C_F_POINTER(NoahmpIO_cptr%JDS,    NoahmpIO%JDS)
    call C_F_POINTER(NoahmpIO_cptr%JDE,    NoahmpIO%JDE)
    call C_F_POINTER(NoahmpIO_cptr%KDS,    NoahmpIO%KDS)
    call C_F_POINTER(NoahmpIO_cptr%KDE,    NoahmpIO%KDE)

    call C_F_POINTER(NoahmpIO_cptr%IMS,    NoahmpIO%IMS)
    call C_F_POINTER(NoahmpIO_cptr%IME,    NoahmpIO%IME)
    call C_F_POINTER(NoahmpIO_cptr%JMS,    NoahmpIO%JMS)
    call C_F_POINTER(NoahmpIO_cptr%JME,    NoahmpIO%JME)
    call C_F_POINTER(NoahmpIO_cptr%KMS,    NoahmpIO%KMS)
    call C_F_POINTER(NoahmpIO_cptr%KME,    NoahmpIO%KME)

    call C_F_POINTER(NoahmpIO_cptr%ITS,    NoahmpIO%ITS)
    call C_F_POINTER(NoahmpIO_cptr%ITE,    NoahmpIO%ITE)
    call C_F_POINTER(NoahmpIO_cptr%JTS,    NoahmpIO%JTS)
    call C_F_POINTER(NoahmpIO_cptr%JTE,    NoahmpIO%JTE)
    call C_F_POINTER(NoahmpIO_cptr%KTS,    NoahmpIO%KTS)
    call C_F_POINTER(NoahmpIO_cptr%KTE,    NoahmpIO%KTE)
  end subroutine NoahmpIOTypeInit_fi

  subroutine NoahmpIOVarInitDefault_fi(NoahmpIO_cptr) bind(C, name="NoahmpIOVarInitDefault_fi")
    use, intrinsic :: iso_c_binding
    implicit none
    type(NoahmpIO_type_fi), intent(inout) :: NoahmpIO_cptr

    call NoahmpIOVarInitDefault(NoahmpIO)

    NoahmpIO_cptr%XLAT = C_LOC(NoahmpIO%XLAT)
    NoahmpIO_cptr%WSLAKEXY = C_LOC(NoahmpIO%WSLAKEXY)
  end subroutine NoahmpIOVarInitDefault_fi

  subroutine NoahmpInitMain_fi(NoahmpIO_cptr) bind(C, name="NoahmpInitMain_fi")
    use, intrinsic :: iso_c_binding 
    implicit none 
    type(NoahmpIO_type_fi), intent(inout) :: NoahmpIO_cptr
    call NoahmpInitMain(NoahmpIO)
  end subroutine NoahmpInitMain_fi

  subroutine NoahmpReadTable_fi(NoahmpIO_cptr) bind(C, name="NoahmpReadTable_fi")
    use, intrinsic :: iso_c_binding 
    implicit none 
    type(NoahmpIO_type_fi), intent(inout) :: NoahmpIO_cptr
    call NoahmpReadTable(NoahmpIO)
  end subroutine NoahmpReadTable_fi

  subroutine NoahmpReadNamelist_fi(NoahmpIO_cptr) bind(C, name="NoahmpReadNamelist_fi")
    use, intrinsic :: iso_c_binding 
    implicit none 
    type(NoahmpIO_type_fi), intent(inout) :: NoahmpIO_cptr
    call NoahmpReadNamelist(NoahmpIO)
  end subroutine NoahmpReadNamelist_fi

end module NoahmpIO_fi
