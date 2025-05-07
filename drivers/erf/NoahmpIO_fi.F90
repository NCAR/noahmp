module NoahmpIO_fi

  use NoahmpIOVarType, ONLY: NoahmpIO_type
  use NoahmpIOVarInitMod, ONLY: NoahmpIOVarInitDefault
  use NoahmpInitMainMod, ONLY: NoahmpInitMain
  use NoahmpReadNamelistMod, ONLY: NoahmpReadNamelist
  use NoahmpReadTableMod, ONLY: NoahmpReadTable
  use NoahmpReadLandMod, ONLY: NoahmpReadLandHeader, NoahmpReadLandMain
  use NoahmpDriverMainMod, ONLY: NoahmpDriverMain

  use  iso_c_binding

  implicit none

  ! ---------------------------------------------------------------------------
  ! Public variable of NoahmpIO_type
  ! --------------------------------------------------------------------------- 
  type(NoahmpIO_type), save, target, public, allocatable, dimension(:) :: NoahmpIO

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
    type(C_PTR) :: llanduse, rank, blkid
    type(C_PTR) :: XLAT, WSLAKEXY
    type(C_PTR) :: U_PHY, T_PHY, V_PHY, QV_CURR
    type(C_PTR) :: SHBXY, EVBXY
  end type NoahmpIO_type_fi

contains

  subroutine NoahmpIOScalarInitDefault_fi(NoahmpIO_cptr) bind(C, name="NoahmpIOScalarInitDefault_fi")
    use  iso_c_binding, only : C_INT
    implicit none
    type(NoahmpIO_type_fi), intent(inout) :: NoahmpIO_cptr
    integer(C_INT), pointer :: bid
    
    call C_F_POINTER(NoahmpIO_cptr%BLKID, bid)

    call C_F_POINTER(NoahmpIO_cptr%XSTART, NoahmpIO(bid)%XSTART)
    call C_F_POINTER(NoahmpIO_cptr%XEND, NoahmpIO(bid)%XEND)
    call C_F_POINTER(NoahmpIO_cptr%YSTART, NoahmpIO(bid)%YSTART)
    call C_F_POINTER(NoahmpIO_cptr%YEND, NoahmpIO(bid)%YEND)

    call C_F_POINTER(NoahmpIO_cptr%NSOIL, NoahmpIO(bid)%NSOIL)
    call C_F_POINTER(NoahmpIO_cptr%NSNOW, NoahmpIO(bid)%NSNOW)

    call C_F_POINTER(NoahmpIO_cptr%IDS, NoahmpIO(bid)%IDS)
    call C_F_POINTER(NoahmpIO_cptr%IDE, NoahmpIO(bid)%IDE)
    call C_F_POINTER(NoahmpIO_cptr%JDS, NoahmpIO(bid)%JDS)
    call C_F_POINTER(NoahmpIO_cptr%JDE, NoahmpIO(bid)%JDE)
    call C_F_POINTER(NoahmpIO_cptr%KDS, NoahmpIO(bid)%KDS)
    call C_F_POINTER(NoahmpIO_cptr%KDE, NoahmpIO(bid)%KDE)

    call C_F_POINTER(NoahmpIO_cptr%IMS, NoahmpIO(bid)%IMS)
    call C_F_POINTER(NoahmpIO_cptr%IME, NoahmpIO(bid)%IME)
    call C_F_POINTER(NoahmpIO_cptr%JMS, NoahmpIO(bid)%JMS)
    call C_F_POINTER(NoahmpIO_cptr%JME, NoahmpIO(bid)%JME)
    call C_F_POINTER(NoahmpIO_cptr%KMS, NoahmpIO(bid)%KMS)
    call C_F_POINTER(NoahmpIO_cptr%KME, NoahmpIO(bid)%KME)

    call C_F_POINTER(NoahmpIO_cptr%ITS, NoahmpIO(bid)%ITS)
    call C_F_POINTER(NoahmpIO_cptr%ITE, NoahmpIO(bid)%ITE)
    call C_F_POINTER(NoahmpIO_cptr%JTS, NoahmpIO(bid)%JTS)
    call C_F_POINTER(NoahmpIO_cptr%JTE, NoahmpIO(bid)%JTE)
    call C_F_POINTER(NoahmpIO_cptr%KTS, NoahmpIO(bid)%KTS)
    call C_F_POINTER(NoahmpIO_cptr%KTE, NoahmpIO(bid)%KTE)

    call C_F_POINTER(NoahmpIO_cptr%ITIMESTEP, NoahmpIO(bid)%ITIMESTEP)
    call C_F_POINTER(NoahmpIO_cptr%NTIME, NoahmpIO(bid)%NTIME)

    call C_F_POINTER(NoahmpIO_cptr%LLANDUSE, NoahmpIO(bid)%LLANDUSE)
    NoahmpIO(bid)%LLANDUSE => NoahmpIO(bid)%LLANDUSE(1:256)

    call C_F_POINTER(NoahmpIO_cptr%RANK, NoahmpIO(bid)%RANK)
  end subroutine NoahmpIOScalarInitDefault_fi

  subroutine NoahmpIOVarInitDefault_fi(NoahmpIO_cptr) bind(C, name="NoahmpIOVarInitDefault_fi")
    use  iso_c_binding, only : C_INT
    implicit none
    type(NoahmpIO_type_fi), intent(inout) :: NoahmpIO_cptr
    integer(C_INT), pointer :: bid
    
    call C_F_POINTER(NoahmpIO_cptr%BLKID, bid)

    call NoahmpIOVarInitDefault(NoahmpIO(bid))

    NoahmpIO_cptr%XLAT = C_LOC(NoahmpIO(bid)%XLAT)
    NoahmpIO_cptr%WSLAKEXY = C_LOC(NoahmpIO(bid)%WSLAKEXY)
    NoahmpIO_cptr%T_PHY = C_LOC(NoahmpIO(bid)%T_PHY)
    NoahmpIO_cptr%U_PHY = C_LOC(NoahmpIO(bid)%U_PHY)
    NoahmpIO_cptr%V_PHY = C_LOC(NoahmpIO(bid)%V_PHY)
    NoahmpIO_cptr%QV_CURR = C_LOC(NoahmpIO(bid)%QV_CURR)
    NoahmpIO_cptr%SHBXY = C_LOC(NoahmpIO(bid)%SHBXY)
    NoahmpIO_cptr%EVBXY = C_LOC(NoahmpIO(bid)%EVBXY)
  end subroutine NoahmpIOVarInitDefault_fi

  subroutine NoahmpInitMain_fi(NoahmpIO_cptr) bind(C, name="NoahmpInitMain_fi")
    use  iso_c_binding, only : C_INT
    implicit none 
    type(NoahmpIO_type_fi), intent(inout) :: NoahmpIO_cptr
    integer(C_INT), pointer :: bid    
    call C_F_POINTER(NoahmpIO_cptr%BLKID, bid)
    call NoahmpInitMain(NoahmpIO(bid))
  end subroutine NoahmpInitMain_fi

  subroutine NoahmpReadTable_fi(NoahmpIO_cptr) bind(C, name="NoahmpReadTable_fi")
    use  iso_c_binding, only : C_INT
    implicit none 
    type(NoahmpIO_type_fi), intent(inout) :: NoahmpIO_cptr
    integer(C_INT), pointer :: bid    
    call C_F_POINTER(NoahmpIO_cptr%BLKID, bid)
    call NoahmpReadTable(NoahmpIO(bid))
  end subroutine NoahmpReadTable_fi

  subroutine NoahmpReadNamelist_fi(NoahmpIO_cptr) bind(C, name="NoahmpReadNamelist_fi")
    use  iso_c_binding, only : C_INT 
    implicit none 
    type(NoahmpIO_type_fi), intent(inout) :: NoahmpIO_cptr
    integer(C_INT), pointer :: bid
    call C_F_POINTER(NoahmpIO_cptr%BLKID, bid)
    call NoahmpReadNamelist(NoahmpIO(bid))
  end subroutine NoahmpReadNamelist_fi

  subroutine NoahmpReadLandHeader_fi(NoahmpIO_cptr) bind(C, name="NoahmpReadLandHeader_fi")
    use  iso_c_binding, only : C_INT 
    implicit none 
    type(NoahmpIO_type_fi), intent(inout) :: NoahmpIO_cptr
    integer(C_INT), pointer :: bid
    call C_F_POINTER(NoahmpIO_cptr%BLKID, bid)
    call NoahmpReadLandHeader(NoahmpIO(bid))
  end subroutine NoahmpReadLandHeader_fi

  subroutine NoahmpReadLandMain_fi(NoahmpIO_cptr) bind(C, name="NoahmpReadLandMain_fi")
    use  iso_c_binding, only : C_INT 
    implicit none 
    type(NoahmpIO_type_fi), intent(inout) :: NoahmpIO_cptr
    integer(C_INT), pointer :: bid
    call C_F_POINTER(NoahmpIO_cptr%BLKID, bid)
    call NoahmpReadLandMain(NoahmpIO(bid))
  end subroutine NoahmpReadLandMain_fi

  subroutine NoahmpDriverMain_fi(NoahmpIO_cptr) bind(C, name="NoahmpDriverMain_fi")
    use iso_c_binding, only : C_INT 
    implicit none 
    type(NoahmpIO_type_fi), intent(inout) :: NoahmpIO_cptr
    integer(C_INT), pointer :: bid
    call C_F_POINTER(NoahmpIO_cptr%BLKID, bid)
    call NoahmpDriverMain(NoahmpIO(bid))
  end subroutine NoahmpDriverMain_fi

  subroutine NoahmpIOTypeVectInit_fi(NBlocks) bind(C, name="NoahmpIOTypeVectInit_fi")
    use iso_c_binding, only : C_INT
    implicit none
    integer(C_INT), intent(in) :: NBlocks
    allocate(NoahmpIO(0:NBlocks-1))
  end subroutine NoahmpIOTypeVectInit_fi

end module NoahmpIO_fi
