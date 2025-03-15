module NoahmpIO_fi

  use NoahmpIOVarType, ONLY: NoahmpIO_type
  use NoahmpIOVarInitMod, ONLY: NoahmpIOVarInitDefault
  use NoahmpInitMainMod, ONLY: NoahmpInitMain
  use NoahmpReadNamelistMod, ONLY: NoahmpReadNamelist
  use NoahmpReadTableMod, ONLY: NoahmpReadTable
  use NoahmpReadLandMod, ONLY: NoahmpReadLandHeader, NoahmpReadLandMain
  use NoahmpDriverMainMod, ONLY: NoahmpDriverMain

  use, intrinsic :: iso_c_binding

  implicit none

  ! ---------------------------------------------------------------------------
  ! Public variable of NoahmpIO_type
  ! --------------------------------------------------------------------------- 
  type(NoahmpIO_type), save, target, public :: NoahmpIO

  ! ---------------------------------------------------------------------------
  ! Mirror of extern C struct
  !
  ! DEVNOTE :: The order variables between C struct and the binded Fortran type
  !            should be consistent for memory managment. 
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
    type(C_PTR) :: itimestep, ntime
    type(C_PTR) :: llanduse
    type(C_PTR) :: XLAT, WSLAKEXY
    type(C_PTR) :: U_PHY, T_PHY, V_PHY, QV_CURR
    type(C_PTR) :: SHBXY, EVBXY
  end type NoahmpIO_type_fi

contains

  subroutine NoahmpIOTypeInit_fi(NoahmpIO_cptr) bind(C, name="NoahmpIOTypeInit_fi")
    use, intrinsic :: iso_c_binding
    implicit none
    type(NoahmpIO_type_fi), intent(inout) :: NoahmpIO_cptr

    call C_F_POINTER(NoahmpIO_cptr%XSTART, NoahmpIO%XSTART)
    call C_F_POINTER(NoahmpIO_cptr%XEND, NoahmpIO%XEND)
    call C_F_POINTER(NoahmpIO_cptr%YSTART, NoahmpIO%YSTART)
    call C_F_POINTER(NoahmpIO_cptr%YEND, NoahmpIO%YEND)

    call C_F_POINTER(NoahmpIO_cptr%NSOIL, NoahmpIO%NSOIL)
    call C_F_POINTER(NoahmpIO_cptr%NSNOW, NoahmpIO%NSNOW)

    call C_F_POINTER(NoahmpIO_cptr%IDS, NoahmpIO%IDS)
    call C_F_POINTER(NoahmpIO_cptr%IDE, NoahmpIO%IDE)
    call C_F_POINTER(NoahmpIO_cptr%JDS, NoahmpIO%JDS)
    call C_F_POINTER(NoahmpIO_cptr%JDE, NoahmpIO%JDE)
    call C_F_POINTER(NoahmpIO_cptr%KDS, NoahmpIO%KDS)
    call C_F_POINTER(NoahmpIO_cptr%KDE, NoahmpIO%KDE)

    call C_F_POINTER(NoahmpIO_cptr%IMS, NoahmpIO%IMS)
    call C_F_POINTER(NoahmpIO_cptr%IME, NoahmpIO%IME)
    call C_F_POINTER(NoahmpIO_cptr%JMS, NoahmpIO%JMS)
    call C_F_POINTER(NoahmpIO_cptr%JME, NoahmpIO%JME)
    call C_F_POINTER(NoahmpIO_cptr%KMS, NoahmpIO%KMS)
    call C_F_POINTER(NoahmpIO_cptr%KME, NoahmpIO%KME)

    call C_F_POINTER(NoahmpIO_cptr%ITS, NoahmpIO%ITS)
    call C_F_POINTER(NoahmpIO_cptr%ITE, NoahmpIO%ITE)
    call C_F_POINTER(NoahmpIO_cptr%JTS, NoahmpIO%JTS)
    call C_F_POINTER(NoahmpIO_cptr%JTE, NoahmpIO%JTE)
    call C_F_POINTER(NoahmpIO_cptr%KTS, NoahmpIO%KTS)
    call C_F_POINTER(NoahmpIO_cptr%KTE, NoahmpIO%KTE)

    call C_F_POINTER(NoahmpIO_cptr%ITIMESTEP, NoahmpIO%ITIMESTEP)
    call C_F_POINTER(NoahmpIO_cptr%NTIME, NoahmpIO%NTIME)

    call C_F_POINTER(NoahmpIO_cptr%LLANDUSE, NoahmpIO%LLANDUSE)
    NoahmpIO%LLANDUSE => NoahmpIO%LLANDUSE(1:256)
  end subroutine NoahmpIOTypeInit_fi

  subroutine NoahmpIOVarInitDefault_fi(NoahmpIO_cptr) bind(C, name="NoahmpIOVarInitDefault_fi")
    use, intrinsic :: iso_c_binding
    implicit none
    type(NoahmpIO_type_fi), intent(inout) :: NoahmpIO_cptr

    call NoahmpIOVarInitDefault(NoahmpIO)

    NoahmpIO_cptr%XLAT = C_LOC(NoahmpIO%XLAT)
    NoahmpIO_cptr%WSLAKEXY = C_LOC(NoahmpIO%WSLAKEXY)
    
    NoahmpIO_cptr%T_PHY = C_LOC(NoahmpIO%T_PHY)
    NoahmpIO_cptr%U_PHY = C_LOC(NoahmpIO%U_PHY)
    NoahmpIO_cptr%V_PHY = C_LOC(NoahmpIO%V_PHY)
    NoahmpIO_cptr%QV_CURR = C_LOC(NoahmpIO%QV_CURR)
    NoahmpIO_cptr%SHBXY = C_LOC(NoahmpIO%SHBXY)
    NoahmpIO_cptr%EVBXY = C_LOC(NoahmpIO%EVBXY)
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

  subroutine NoahmpReadLandHeader_fi(NoahmpIO_cptr) bind(C, name="NoahmpReadLandHeader_fi")
    use, intrinsic :: iso_c_binding 
    implicit none 
    type(NoahmpIO_type_fi), intent(inout) :: NoahmpIO_cptr
    call NoahmpReadLandHeader(NoahmpIO)
  end subroutine NoahmpReadLandHeader_fi

  subroutine NoahmpReadLandMain_fi(NoahmpIO_cptr) bind(C, name="NoahmpReadLandMain_fi")
    use, intrinsic :: iso_c_binding 
    implicit none 
    type(NoahmpIO_type_fi), intent(inout) :: NoahmpIO_cptr
    call NoahmpReadLandMain(NoahmpIO)
  end subroutine NoahmpReadLandMain_fi

  subroutine NoahmpDriverMain_fi(NoahmpIO_cptr) bind(C, name="NoahmpDriverMain_fi")
    use, intrinsic :: iso_c_binding 
    implicit none 
    type(NoahmpIO_type_fi), intent(inout) :: NoahmpIO_cptr
    call NoahmpDriverMain(NoahmpIO)
  end subroutine NoahmpDriverMain_fi

end module NoahmpIO_fi
