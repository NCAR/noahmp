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

  subroutine NoahmpIOVarInitDefault_fi(NoahmpIO_cptr) bind(C, name="NoahmpIOVarInitDefault_fi")
    use, intrinsic :: iso_c_binding
    implicit none
    type(NoahmpIO_type_fi), intent(inout) :: NoahmpIO_cptr

    call C_F_POINTER(NoahmpIO_cptr%XSTART, NoahmpIO%XSTART)
    call C_F_POINTER(NoahmpIO_cptr%XEND,   NoahmpIO%XEND)
    call C_F_POINTER(NoahmpIO_cptr%YSTART, NoahmpIO%YSTART)
    call C_F_POINTER(NoahmpIO_cptr%YEND,   NoahmpIO%YEND)
    call C_F_POINTER(NoahmpIO_cptr%KDS,    NoahmpIO%KDS)
    call C_F_POINTER(NoahmpIO_cptr%KDE,    NoahmpIO%KDE)
    call C_F_POINTER(NoahmpIO_cptr%NSOIL,  NoahmpIO%NSOIL)
    call C_F_POINTER(NoahmpIO_cptr%NSNOW,  NoahmpIO%NSNOW)

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
