module NoahmpIO_fi

  use NoahmpIOVarType, ONLY: NoahmpIO_type
  use NoahmpIOVarInitMod, ONLY: NoahmpIOVarInitDefault
  use NoahmpInitMainMod, ONLY: NoahmpInitMain
  use NoahmpReadNamelistMod, ONLY: NoahmpReadNamelist
  use NoahmpReadTableMod, ONLY: NoahmpReadTable
  use NoahmpReadLandMod, ONLY: NoahmpReadLandHeader, NoahmpReadLandMain
  use NoahmpDriverMainMod, ONLY: NoahmpDriverMain
  use NoahmpWriteLandMod, ONLY: NoahmpWriteLand
  use NoahmpWriteRestartMod, ONLY: NoahmpWriteRestart
  use NoahmpReadRestartMod, ONLY: NoahmpReadRestart

  use  iso_c_binding

  implicit none

  ! ---------------------------------------------------------------------------
  ! Public type for NoahmpIO_level
  ! --------------------------------------------------------------------------- 
  type, public :: NoahmpIO_level
    type(NoahmpIO_type), allocatable :: NoahmpIO(:)
  end type NoahmpIO_level

  ! ---------------------------------------------------------------------------
  ! Valid AMR level range for NoahmpIO_vect. Change both bounds here (and the
  ! matching C++ caller) if the number of supported levels changes.
  ! ---------------------------------------------------------------------------
  integer, parameter, public :: NLEVEL_MIN = 0
  integer, parameter, public :: NLEVEL_MAX = 2

  ! ---------------------------------------------------------------------------
  ! Public declaration for NoahmpIO array per level
  ! ---------------------------------------------------------------------------
  type(NoahmpIO_level), save, target, public :: NoahmpIO_vect(NLEVEL_MIN:NLEVEL_MAX)

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
    type(C_PTR) :: rank, blkid, level
    type(C_PTR) :: comm
    type(C_PTR) :: DTBL
    type(C_PTR) :: ZLVL
    type(C_PTR) :: XLAT, WSLAKEXY
    type(C_PTR) :: U_PHY, T_PHY, V_PHY, QV_CURR
    type(C_PTR) :: HFX, LH
    type(C_PTR) :: SWDOWN, GLW, TSK, EMISS
    type(C_PTR) :: ALBSFCDIRXY, ALBSFCDIFXY
    type(C_PTR) :: COSZEN, P8W
    type(C_PTR) :: TAU_EW, TAU_NS
  end type NoahmpIO_type_fi

contains

  ! ---------------------------------------------------------------------------
  ! Resolve (level, bid) from the incoming C struct and return a pointer to the
  ! corresponding block. Validates the C pointers and both indices, aborting
  ! with a clear message instead of indexing module-global state out of bounds
  ! (which is silent memory corruption in a release build). Every bind(C) entry
  ! point routes through this helper.
  ! ---------------------------------------------------------------------------
  subroutine resolve_block(NoahmpIO_cptr, blk, level, bid, require_init)
    use iso_c_binding, only : C_INT, c_associated
    implicit none
    type(NoahmpIO_type_fi), intent(in)           :: NoahmpIO_cptr
    type(NoahmpIO_type),    pointer, intent(out)  :: blk
    integer(C_INT),         intent(out)           :: level, bid
    logical,                intent(in),  optional :: require_init
    integer(C_INT), pointer :: plevel, pbid
    logical :: check_init

    ! By default require that NoahmpIOScalarInitDefault_fi has already wired the
    ! block's scalar pointers to C++ memory. ScalarInitDefault itself is the call
    ! that performs that wiring, so it passes require_init=.false.
    check_init = .true.
    if (present(require_init)) check_init = require_init

    if (.not. c_associated(NoahmpIO_cptr%LEVEL) .or. &
        .not. c_associated(NoahmpIO_cptr%BLKID)) then
       write(0,*) "NoahmpIO_fi: LEVEL or BLKID C pointer is not associated"
       error stop
    end if

    call C_F_POINTER(NoahmpIO_cptr%LEVEL, plevel)
    call C_F_POINTER(NoahmpIO_cptr%BLKID, pbid)
    level = plevel
    bid   = pbid

    if (level < NLEVEL_MIN .or. level > NLEVEL_MAX) then
       write(0,*) "NoahmpIO_fi: level out of range [", NLEVEL_MIN, ",", &
                  NLEVEL_MAX, "]:", level
       error stop
    end if
    if (.not. allocated(NoahmpIO_vect(level)%NoahmpIO)) then
       write(0,*) "NoahmpIO_fi: level not initialized (call NoahmpIOTypeVectInit "// &
                  "first), level:", level
       error stop
    end if
    if (bid < lbound(NoahmpIO_vect(level)%NoahmpIO, 1) .or. &
        bid > ubound(NoahmpIO_vect(level)%NoahmpIO, 1)) then
       write(0,*) "NoahmpIO_fi: blkid out of range [", &
                  lbound(NoahmpIO_vect(level)%NoahmpIO, 1), ",", &
                  ubound(NoahmpIO_vect(level)%NoahmpIO, 1), "] for level", &
                  level, ", blkid:", bid
       error stop
    end if

    blk => NoahmpIO_vect(level)%NoahmpIO(bid)

    ! DTBL is associated together with every other scalar pointer at the end of
    ! NoahmpIOScalarInitDefault_fi, so it is a reliable proxy for "scalars wired".
    ! This is safe because the pointer components default to => null() in the type
    ! definition, giving them well-defined association status before that call.
    if (check_init .and. .not. associated(blk%DTBL)) then
       write(0,*) "NoahmpIO_fi: block scalars not initialized (call "// &
                  "NoahmpIOScalarInitDefault first) for level", level, ", blkid", bid
       error stop
    end if
  end subroutine resolve_block

  subroutine NoahmpIOScalarInitDefault_fi(NoahmpIO_cptr) bind(C, name="NoahmpIOScalarInitDefault_fi")
    use  iso_c_binding, only : C_INT
    implicit none
    type(NoahmpIO_type_fi), intent(inout) :: NoahmpIO_cptr
    type(NoahmpIO_type), pointer :: blk
    integer(C_INT) :: level, bid

    ! This call performs the scalar wiring, so do not require it to be done yet.
    call resolve_block(NoahmpIO_cptr, blk, level, bid, require_init=.false.)

    call C_F_POINTER(NoahmpIO_cptr%BLKID, blk%blkid)
    call C_F_POINTER(NoahmpIO_cptr%LEVEL, blk%level)

    call C_F_POINTER(NoahmpIO_cptr%XSTART, blk%XSTART)
    call C_F_POINTER(NoahmpIO_cptr%XEND, blk%XEND)
    call C_F_POINTER(NoahmpIO_cptr%YSTART, blk%YSTART)
    call C_F_POINTER(NoahmpIO_cptr%YEND, blk%YEND)

    call C_F_POINTER(NoahmpIO_cptr%NSOIL, blk%NSOIL)
    call C_F_POINTER(NoahmpIO_cptr%NSNOW, blk%NSNOW)

    call C_F_POINTER(NoahmpIO_cptr%IDS, blk%IDS)
    call C_F_POINTER(NoahmpIO_cptr%IDE, blk%IDE)
    call C_F_POINTER(NoahmpIO_cptr%JDS, blk%JDS)
    call C_F_POINTER(NoahmpIO_cptr%JDE, blk%JDE)
    call C_F_POINTER(NoahmpIO_cptr%KDS, blk%KDS)
    call C_F_POINTER(NoahmpIO_cptr%KDE, blk%KDE)

    call C_F_POINTER(NoahmpIO_cptr%IMS, blk%IMS)
    call C_F_POINTER(NoahmpIO_cptr%IME, blk%IME)
    call C_F_POINTER(NoahmpIO_cptr%JMS, blk%JMS)
    call C_F_POINTER(NoahmpIO_cptr%JME, blk%JME)
    call C_F_POINTER(NoahmpIO_cptr%KMS, blk%KMS)
    call C_F_POINTER(NoahmpIO_cptr%KME, blk%KME)

    call C_F_POINTER(NoahmpIO_cptr%ITS, blk%ITS)
    call C_F_POINTER(NoahmpIO_cptr%ITE, blk%ITE)
    call C_F_POINTER(NoahmpIO_cptr%JTS, blk%JTS)
    call C_F_POINTER(NoahmpIO_cptr%JTE, blk%JTE)
    call C_F_POINTER(NoahmpIO_cptr%KTS, blk%KTS)
    call C_F_POINTER(NoahmpIO_cptr%KTE, blk%KTE)

    call C_F_POINTER(NoahmpIO_cptr%ITIMESTEP, blk%ITIMESTEP)
    call C_F_POINTER(NoahmpIO_cptr%NTIME, blk%NTIME)

    call C_F_POINTER(NoahmpIO_cptr%RANK, blk%RANK)
    call C_F_POINTER(NoahmpIO_cptr%COMM, blk%COMM)

    call C_F_POINTER(NoahmpIO_cptr%DTBL, blk%DTBL)
    call C_F_POINTER(NoahmpIO_cptr%ZLVL, blk%zlvl)

  end subroutine NoahmpIOScalarInitDefault_fi

  subroutine NoahmpIOVarInitDefault_fi(NoahmpIO_cptr) bind(C, name="NoahmpIOVarInitDefault_fi")
    use  iso_c_binding, only : C_INT
    implicit none
    type(NoahmpIO_type_fi), intent(inout) :: NoahmpIO_cptr
    type(NoahmpIO_type), pointer :: blk
    integer(C_INT) :: level, bid

    call resolve_block(NoahmpIO_cptr, blk, level, bid)

    call NoahmpIOVarInitDefault(blk)

    NoahmpIO_cptr%XLAT = C_LOC(blk%XLAT)
    NoahmpIO_cptr%WSLAKEXY = C_LOC(blk%WSLAKEXY)
    NoahmpIO_cptr%T_PHY = C_LOC(blk%T_PHY)
    NoahmpIO_cptr%U_PHY = C_LOC(blk%U_PHY)
    NoahmpIO_cptr%V_PHY = C_LOC(blk%V_PHY)
    NoahmpIO_cptr%QV_CURR = C_LOC(blk%QV_CURR)
    NoahmpIO_cptr%HFX = C_LOC(blk%HFX)
    NoahmpIO_cptr%LH = C_LOC(blk%LH)
    NoahmpIO_cptr%SWDOWN = C_LOC(blk%SWDOWN)
    NoahmpIO_cptr%GLW = C_LOC(blk%GLW)
    NoahmpIO_cptr%TSK = C_LOC(blk%TSK)
    NoahmpIO_cptr%EMISS = C_LOC(blk%EMISS)
    NoahmpIO_cptr%ALBSFCDIRXY = C_LOC(blk%ALBSFCDIRXY)
    NoahmpIO_cptr%ALBSFCDIFXY = C_LOC(blk%ALBSFCDIFXY)
    NoahmpIO_cptr%COSZEN = C_LOC(blk%COSZEN)
    NoahmpIO_cptr%P8W = C_LOC(blk%P8W)
    NoahmpIO_cptr%TAU_EW = C_LOC(blk%TAU_EW)
    NoahmpIO_cptr%TAU_NS = C_LOC(blk%TAU_NS)
  end subroutine NoahmpIOVarInitDefault_fi

  subroutine NoahmpInitMain_fi(NoahmpIO_cptr) bind(C, name="NoahmpInitMain_fi")
    use  iso_c_binding, only : C_INT
    implicit none 
    type(NoahmpIO_type_fi), intent(inout) :: NoahmpIO_cptr
    type(NoahmpIO_type), pointer :: blk
    integer(C_INT) :: level, bid
    call resolve_block(NoahmpIO_cptr, blk, level, bid)
    call NoahmpInitMain(blk)
  end subroutine NoahmpInitMain_fi

  subroutine NoahmpReadTable_fi(NoahmpIO_cptr) bind(C, name="NoahmpReadTable_fi")
    use  iso_c_binding, only : C_INT
    implicit none 
    type(NoahmpIO_type_fi), intent(inout) :: NoahmpIO_cptr
    type(NoahmpIO_type), pointer :: blk
    integer(C_INT) :: level, bid
    call resolve_block(NoahmpIO_cptr, blk, level, bid)
    call NoahmpReadTable(blk)
  end subroutine NoahmpReadTable_fi

  subroutine NoahmpReadNamelist_fi(NoahmpIO_cptr) bind(C, name="NoahmpReadNamelist_fi")
    use  iso_c_binding, only : C_INT 
    implicit none 
    type(NoahmpIO_type_fi), intent(inout) :: NoahmpIO_cptr
    type(NoahmpIO_type), pointer :: blk
    integer(C_INT) :: level, bid
    call resolve_block(NoahmpIO_cptr, blk, level, bid)
    call NoahmpReadNamelist(blk)
  end subroutine NoahmpReadNamelist_fi

  subroutine NoahmpReadLandHeader_fi(NoahmpIO_cptr) bind(C, name="NoahmpReadLandHeader_fi")
    use  iso_c_binding, only : C_INT 
    implicit none 
    type(NoahmpIO_type_fi), intent(inout) :: NoahmpIO_cptr
    type(NoahmpIO_type), pointer :: blk
    integer(C_INT) :: level, bid
    call resolve_block(NoahmpIO_cptr, blk, level, bid)
    call NoahmpReadLandHeader(blk)
  end subroutine NoahmpReadLandHeader_fi

  subroutine NoahmpReadLandMain_fi(NoahmpIO_cptr) bind(C, name="NoahmpReadLandMain_fi")
    use  iso_c_binding, only : C_INT 
    implicit none 
    type(NoahmpIO_type_fi), intent(inout) :: NoahmpIO_cptr
    type(NoahmpIO_type), pointer :: blk
    integer(C_INT) :: level, bid
    call resolve_block(NoahmpIO_cptr, blk, level, bid)
    call NoahmpReadLandMain(blk)
  end subroutine NoahmpReadLandMain_fi

  subroutine NoahmpDriverMain_fi(NoahmpIO_cptr) bind(C, name="NoahmpDriverMain_fi")
    use iso_c_binding, only : C_INT 
    implicit none 
    type(NoahmpIO_type_fi), intent(inout) :: NoahmpIO_cptr
    type(NoahmpIO_type), pointer :: blk
    integer(C_INT) :: level, bid
    call resolve_block(NoahmpIO_cptr, blk, level, bid)
    call NoahmpDriverMain(blk)
  end subroutine NoahmpDriverMain_fi

  subroutine NoahmpWriteLand_fi(NoahmpIO_cptr, filenum) bind(C, name="NoahmpWriteLand_fi")
    use iso_c_binding, only : C_INT 
    implicit none 
    type(NoahmpIO_type_fi), intent(inout) :: NoahmpIO_cptr
    integer(C_INT), intent(in) :: filenum
    type(NoahmpIO_type), pointer :: blk
    integer(C_INT) :: level, bid
    call resolve_block(NoahmpIO_cptr, blk, level, bid)
    call NoahmpWriteLand(blk, filenum, SIZE(NoahmpIO_vect(level)%NoahmpIO))
  end subroutine NoahmpWriteLand_fi

  subroutine NoahmpWriteRestart_fi(NoahmpIO_cptr, dir_cptr, dir_len) bind(C, name="NoahmpWriteRestart_fi")
    use iso_c_binding, only : C_INT, C_CHAR
    implicit none
    type(NoahmpIO_type_fi), intent(inout) :: NoahmpIO_cptr
    character(kind=C_CHAR), intent(in) :: dir_cptr(*)
    integer(C_INT), intent(in) :: dir_len
    type(NoahmpIO_type), pointer :: blk
    integer(C_INT) :: level, bid
    integer :: ic
    if (dir_len < 0) then
       write(0,*) "NoahmpWriteRestart_fi: negative dir_len (likely an int "// &
                  "overflow of the C++ string length):", dir_len
       error stop
    end if
    block
      character(len=dir_len) :: dir
      do ic = 1, dir_len
         dir(ic:ic) = dir_cptr(ic)
      end do
      call resolve_block(NoahmpIO_cptr, blk, level, bid)
      call NoahmpWriteRestart(blk, dir, SIZE(NoahmpIO_vect(level)%NoahmpIO))
    end block
  end subroutine NoahmpWriteRestart_fi

  subroutine NoahmpReadRestart_fi(NoahmpIO_cptr, dir_cptr, dir_len) bind(C, name="NoahmpReadRestart_fi")
    use iso_c_binding, only : C_INT, C_CHAR
    implicit none
    type(NoahmpIO_type_fi), intent(inout) :: NoahmpIO_cptr
    character(kind=C_CHAR), intent(in) :: dir_cptr(*)
    integer(C_INT), intent(in) :: dir_len
    type(NoahmpIO_type), pointer :: blk
    integer(C_INT) :: level, bid
    integer :: ic
    if (dir_len < 0) then
       write(0,*) "NoahmpReadRestart_fi: negative dir_len (likely an int "// &
                  "overflow of the C++ string length):", dir_len
       error stop
    end if
    block
      character(len=dir_len) :: dir
      do ic = 1, dir_len
         dir(ic:ic) = dir_cptr(ic)
      end do
      call resolve_block(NoahmpIO_cptr, blk, level, bid)
      call NoahmpReadRestart(blk, dir, SIZE(NoahmpIO_vect(level)%NoahmpIO))
    end block
  end subroutine NoahmpReadRestart_fi

  subroutine NoahmpIOTypeVectInit_fi(level, NBlocks) bind(C, name="NoahmpIOTypeVectInit_fi")
    use iso_c_binding, only : C_INT
    implicit none
    integer(C_INT), intent(in) :: level
    integer(C_INT), intent(in) :: NBlocks
    integer :: ierr
    character(len=256) :: emsg

    if (level < NLEVEL_MIN .or. level > NLEVEL_MAX) then
       write(0,*) "NoahmpIOTypeVectInit_fi: level out of range [", NLEVEL_MIN, &
                  ",", NLEVEL_MAX, "]:", level
       error stop
    end if
    if (NBlocks < 1) then
       write(0,*) "NoahmpIOTypeVectInit_fi: NBlocks must be >= 1, got", NBlocks
       error stop
    end if
    if (allocated(NoahmpIO_vect(level)%NoahmpIO)) then
       ! Same size -> idempotent no-op; different size -> refuse (would discard
       ! existing block state). The C++ caller already forbids a second resize,
       ! so this is a defense-in-depth guard.
       if (size(NoahmpIO_vect(level)%NoahmpIO) == NBlocks) return
       write(0,*) "NoahmpIOTypeVectInit_fi: level", level, "already initialized "// &
                  "with", size(NoahmpIO_vect(level)%NoahmpIO), &
                  "blocks; refusing to re-init with", NBlocks
       error stop
    end if

    allocate(NoahmpIO_vect(level)%NoahmpIO(0:NBlocks-1), stat=ierr, errmsg=emsg)
    if (ierr /= 0) then
       write(0,*) "NoahmpIOTypeVectInit_fi: allocate failed: ", trim(emsg)
       error stop
    end if
  end subroutine NoahmpIOTypeVectInit_fi

  ! ---------------------------------------------------------------------------
  ! ABI guard (run time): report sizeof(NoahmpIO_type_fi) as seen by the Fortran
  ! compiler so the C++ side can confirm the two bind(C) struct definitions agree
  ! byte-for-byte (see NoahmpIO_vector::resize). This catches a CHANGE IN MEMBER
  ! COUNT or unexpected padding. It does NOT catch a pure reordering (the struct
  ! is all C_PTR, so any permutation has the same byte size) -- member ORDER is
  ! checked separately by NoahmpIOCheckMemberOrder_fi below.
  ! ---------------------------------------------------------------------------
  function NoahmpIOTypeFiSize_fi() bind(C, name="NoahmpIOTypeFiSize_fi") result(nbytes)
    use iso_c_binding, only : C_SIZE_T, C_SIZEOF
    implicit none
    integer(C_SIZE_T) :: nbytes
    type(NoahmpIO_type_fi) :: probe
    nbytes = C_SIZEOF(probe)
  end function NoahmpIOTypeFiSize_fi

  ! ---------------------------------------------------------------------------
  ! Precision guard (run time): report sizeof(real(c_kind_noahmp)) as seen by the
  ! Fortran compiler. NoahmpIO_type_fi is all C_PTR, so its byte size is the same
  ! whether the coupling reals are single or double precision; NoahmpIOTypeFiSize_fi
  ! therefore CANNOT detect a float/double mismatch between the two languages
  ! (e.g. DOUBLE_PREC applied to only one compiler). The C++ side compares this
  ! against sizeof(noahmp_real) to catch exactly that case (see resize()).
  ! ---------------------------------------------------------------------------
  function NoahmpRealSize_fi() bind(C, name="NoahmpRealSize_fi") result(nbytes)
    use iso_c_binding, only : C_SIZE_T, C_SIZEOF
    use Machine, only : c_kind_noahmp
    implicit none
    integer(C_SIZE_T) :: nbytes
    real(kind=c_kind_noahmp) :: probe
    nbytes = C_SIZEOF(probe)
  end function NoahmpRealSize_fi

  ! ---------------------------------------------------------------------------
  ! Member-ORDER guard (run time). The size/precision guards above cannot detect
  ! a field reordered on only one side, because every member is a pointer and any
  ! permutation has the same byte size -- yet a reorder silently crosses the whole
  ! C<->Fortran pointer map. This handshake catches it:
  !
  !   The C++ caller fills slot i (in C++ declaration order, i = 0..N-1) with a
  !   pointer to an int holding the value i, then passes the struct here. For each
  !   named member, in the agreed canonical order, we dereference it and confirm
  !   it sees its own canonical index. If the C++ and Fortran field orders agree,
  !   member M (canonical index k) sits at the same offset on both sides and reads
  !   k. If either side is reordered, M reads a different slot's value -> mismatch.
  !
  ! Returns the canonical index of the FIRST member that fails, or -1 if every
  ! member is in agreement. The canonical order below MUST match the C++ struct
  ! declaration order in NoahmpIO.H (and this Fortran type's order); all three are
  ! updated together when a field is added.
  ! ---------------------------------------------------------------------------
  function NoahmpIOCheckMemberOrder_fi(probe) bind(C, name="NoahmpIOCheckMemberOrder_fi") result(bad)
    use iso_c_binding, only : C_INT, C_PTR, C_F_POINTER, c_associated
    implicit none
    type(NoahmpIO_type_fi), intent(in) :: probe
    integer(C_INT) :: bad

    bad = -1
    call chk(probe%ids,  0); call chk(probe%ide,  1); call chk(probe%jds,  2)
    call chk(probe%jde,  3); call chk(probe%kds,  4); call chk(probe%kde,  5)
    call chk(probe%ims,  6); call chk(probe%ime,  7); call chk(probe%jms,  8)
    call chk(probe%jme,  9); call chk(probe%kms, 10); call chk(probe%kme, 11)
    call chk(probe%its, 12); call chk(probe%ite, 13); call chk(probe%jts, 14)
    call chk(probe%jte, 15); call chk(probe%kts, 16); call chk(probe%kte, 17)
    call chk(probe%xstart, 18); call chk(probe%xend, 19)
    call chk(probe%ystart, 20); call chk(probe%yend, 21)
    call chk(probe%nsoil, 22); call chk(probe%nsnow, 23)
    call chk(probe%itimestep, 24); call chk(probe%ntime, 25)
    call chk(probe%rank, 26); call chk(probe%blkid, 27); call chk(probe%level, 28)
    call chk(probe%comm, 29)
    call chk(probe%DTBL, 30); call chk(probe%ZLVL, 31)
    call chk(probe%XLAT, 32); call chk(probe%WSLAKEXY, 33)
    call chk(probe%U_PHY, 34); call chk(probe%T_PHY, 35)
    call chk(probe%V_PHY, 36); call chk(probe%QV_CURR, 37)
    call chk(probe%HFX, 38); call chk(probe%LH, 39)
    call chk(probe%SWDOWN, 40); call chk(probe%GLW, 41)
    call chk(probe%TSK, 42); call chk(probe%EMISS, 43)
    call chk(probe%ALBSFCDIRXY, 44); call chk(probe%ALBSFCDIFXY, 45)
    call chk(probe%COSZEN, 46); call chk(probe%P8W, 47)
    call chk(probe%TAU_EW, 48); call chk(probe%TAU_NS, 49)

  contains

    subroutine chk(p, idx)
      type(C_PTR),    intent(in) :: p
      integer(C_INT), intent(in) :: idx
      integer(C_INT), pointer    :: ip
      if (bad >= 0) return                 ! already found a mismatch; stop at the first
      if (.not. c_associated(p)) then
         bad = idx
         return
      end if
      call C_F_POINTER(p, ip)
      if (ip /= idx) bad = idx
    end subroutine chk

  end function NoahmpIOCheckMemberOrder_fi

end module NoahmpIO_fi
