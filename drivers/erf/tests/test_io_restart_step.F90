! ===========================================================================
! test_io_restart_step -- checkpoint/restart REPRODUCIBILITY across time steps,
! for each major Noah-MP subsystem.
!
! test_io_restart only proves the writer/reader round-trip the 66 serialized
! arrays bit-for-bit. It cannot catch a prognostic variable that the driver
! evolves but the restart set omits: that variable round-trips its (unwritten)
! self trivially, yet a real restarted trajectory silently diverges. This test
! drives NoahmpDriverMain and compares trajectories:
!
!   A (continuous):  cold-init -> step 1..N
!   B (checkpoint):  cold-init -> step 1..N-1 -> WriteRestart
!   C (restart):     cold-init -> ReadRestart(B) -> step N
!
! A and B are identical legs, so C enters step N with B's post-step-(N-1) state
! restored from disk. If the 66-var restart set fully determines the step-N
! evolution (identical forcing + itimestep), then C(N) == A(N) bit-for-bit.
!
! A single scenario would leave most prognostic arrays at their cold-init value,
! making their comparison vacuous (0==0). Each selector below instead drives a
! distinct subsystem to LIFE and then (a) asserts, via an activity guard, that
! the subsystem's signature array actually evolved in leg A, and (b) requires
! the full state to match bit-for-bit across the restart:
!
!   land        warm/dry              -> soil temperature + surface fluxes
!   rain        warm + liquid precip  -> infiltration / runoff / canopy water
!   snow        cold + snowfall       -> snowpack, snow layers, albedo aging
!   carbon      DVEG=2, growing       -> dynamic-vegetation carbon pools
!   glacier     ice point + snowfall  -> glacier energy/mass path
!   groundwater subsurface + recharge -> aquifer / water-table state
!
! All scenarios run in ONE process: the library keeps block state in a saved
! module global (NoahmpIO_vect) with a re-init guard and no teardown API, so
! rather than tearing down between scenarios we allocate one fresh block triplet
! per scenario (blocks 3k,3k+1,3k+2) up front and never re-init. argv(1) selects
! a single scenario for debugging; the default ("all") runs the whole suite.
!
! Generic like test_io_driver: synthesizes its own small wrfinput and owns the
! MPI world. Needs staged namelist.erf + NoahmpTable.TBL (CONFIG_INPUTS).
! Exit code 0 = pass.
! ===========================================================================
program test_io_restart_step

  use iso_c_binding,        only : C_INT
  use netcdf
  use Machine,              only : kind_noahmp, c_kind_noahmp
  use NoahmpIO_fi,          only : NoahmpIO_vect, NoahmpIOTypeVectInit_fi
  use NoahmpIOVarType,      only : NoahmpIO_type
  use NoahmpIOVarInitMod,   only : NoahmpIOVarInitDefault
  use NoahmpReadNamelistMod,only : NoahmpReadNamelist
  use NoahmpReadTableMod,   only : NoahmpReadTable
  use NoahmpReadLandMod,    only : NoahmpReadLandHeader, NoahmpReadLandMain
  use NoahmpInitMainMod,    only : NoahmpInitMain
  use NoahmpDriverMainMod,  only : NoahmpDriverMain
  use NoahmpWriteRestartMod,only : NoahmpWriteRestart
  use NoahmpReadRestartMod, only : NoahmpReadRestart
  use NoahmpTestIOSupport

  implicit none

  ! Per-scenario configuration (forcing + option/surface overrides). A negative
  ! override field means "leave the namelist/cold-init value untouched".
  type :: scenario_t
     character(len=16)        :: name
     integer                  :: nsteps
     real(kind=c_kind_noahmp) :: tair, qv, swdown, glw, coszen
     real(kind=c_kind_noahmp) :: rainbl, sr          ! total precip [mm/step], snow ratio
     integer                  :: ivgtyp              ! surface-type override (<0 none)
     integer                  :: iopt_dveg           ! dynamic-veg option override (<0 none)
     integer                  :: iopt_runsub         ! subsurface-runoff override (<0 none)
     character(len=8)         :: signature           ! array that MUST evolve in leg A
  end type scenario_t

  integer, parameter :: NX = 4, NY = 4, NSOIL = 4
  character(len=16), parameter :: ALL_SCENARIOS(6) = &
       [ character(len=16) :: "land", "rain", "snow", "carbon", "glacier", "groundwater" ]

  integer           :: comm, k, nscn
  character(len=32) :: which
  type(scenario_t)  :: scn

  call get_command_argument(1, which)
  if (len_trim(which) == 0) which = "all"

  call tio_reset()
  call tio_mpi_init(comm)

  if (trim(which) == "all") then
     ! One fresh block triplet per scenario (blocks 3k,3k+1,3k+2); the library
     ! has no teardown, so we never re-init -- we index into distinct blocks.
     nscn = size(ALL_SCENARIOS)
     call NoahmpIOTypeVectInit_fi(0_C_INT, int(3*nscn, C_INT))
     do k = 1, nscn
        scn = configure(trim(ALL_SCENARIOS(k)))
        call run_scenario(scn, 3*(k-1), comm)
     end do
  else
     scn = configure(trim(which))
     call NoahmpIOTypeVectInit_fi(0_C_INT, 3_C_INT)
     call run_scenario(scn, 0, comm)
  end if

  call tio_mpi_finalize()
  call tio_finish("test_io_restart_step")

contains

  ! Run one scenario on the block triplet starting at base index k0 (A/B/C =
  ! k0/k0+1/k0+2): continuous 1..N, checkpoint at N-1, restart, step N, then the
  ! activity guard + full bit-exact comparison.
  subroutine run_scenario(s, k0, bcomm)
    type(scenario_t), intent(in) :: s
    integer,          intent(in) :: k0, bcomm

    type(NoahmpIO_type), pointer :: A, B, C
    integer :: it
    character(len=64) :: path, dir
    ! cold-init snapshots of leg A's signature candidates (activity guard).
    real(kind=kind_noahmp), allocatable :: tslb0(:,:,:)
    real(kind=kind_noahmp), allocatable :: smois0(:,:), snow0(:,:), lf0(:,:), zwt0(:,:)

    write(*,'(/,A,A)') "==== restart-step scenario: ", trim(s%name)//" ===="

    path = "wrfin_" // trim(s%name) // ".nc"
    call make_wrfinput(trim(path), NX, NY, NSOIL)

    A => NoahmpIO_vect(0)%NoahmpIO(k0)
    B => NoahmpIO_vect(0)%NoahmpIO(k0+1)
    C => NoahmpIO_vect(0)%NoahmpIO(k0+2)

    call init_block(A, trim(path), bcomm, s)
    call init_block(B, trim(path), bcomm, s)
    call init_block(C, trim(path), bcomm, s)

    call tio_expect(size(A%TSK,1) == NX .and. size(A%TSK,2) == NY, &
                    trim(s%name)//": state arrays span the full domain (non-empty)")

    ! Snapshot leg A cold-init state for the activity guard.
    tslb0  = A%TSLB
    smois0 = A%SMOIS(:,:,1)
    snow0  = A%SNOW
    lf0    = A%LFMASSXY
    zwt0   = A%ZWTXY

    ! leg A: continuous steps 1..N
    do it = 1, s%nsteps
       call step(A, it, s)
    end do

    ! leg B: steps 1..N-1 then checkpoint
    do it = 1, s%nsteps - 1
       call step(B, it, s)
    end do
    dir = "rst_" // trim(s%name)
    call NoahmpWriteRestart(B, trim(dir), 1)

    ! leg C: restart the checkpoint into a fresh block, then step N
    call NoahmpReadRestart(C, trim(dir), 1)
    call step(C, s%nsteps, s)

    ! activity guard: the subsystem must have actually fired in leg A, else a
    ! bit-exact "pass" below would be vacuous.
    call require_active(s, A, tslb0, smois0, snow0, lf0, zwt0)

    ! A(N) vs C(N): bit-exact or the restart set is incomplete.
    write(*,'(A)') "Comparing continuous vs restarted state at final step:"
    ! surface fluxes / skin -- the most sensitive detectors of any input drift
    call cmp2d(s, "HFX",     A%HFX,     C%HFX)
    call cmp2d(s, "LH",      A%LH,      C%LH)
    call cmp2d(s, "QFX",     A%QFX,     C%QFX)
    call cmp2d(s, "GRDFLX",  A%GRDFLX,  C%GRDFLX)
    call cmp2d(s, "TSK",     A%TSK,     C%TSK)
    call cmp2d(s, "EMISS",   A%EMISS,   C%EMISS)
    call cmp2d(s, "ALBEDO",  A%ALBEDO,  C%ALBEDO)
    ! soil (3D)
    call cmp3d(s, "TSLB",    A%TSLB,    C%TSLB)
    call cmp3d(s, "SMOIS",   A%SMOIS,   C%SMOIS)
    call cmp3d(s, "SH2O",    A%SH2O,    C%SH2O)
    ! snow layers (3D) + active-layer count
    call cmp3d(s, "TSNOXY",  A%TSNOXY,  C%TSNOXY)
    call cmp3d(s, "SNICEXY", A%SNICEXY, C%SNICEXY)
    call cmp3d(s, "SNLIQXY", A%SNLIQXY, C%SNLIQXY)
    call cmp3d(s, "ZSNSOXY", A%ZSNSOXY, C%ZSNSOXY)
    call cmp2di(s,"ISNOWXY", A%ISNOWXY, C%ISNOWXY)
    ! snowpack scalars
    call cmp2d(s, "SNOW",    A%SNOW,    C%SNOW)
    call cmp2d(s, "SNOWH",   A%SNOWH,   C%SNOWH)
    call cmp2d(s, "SNOWC",   A%SNOWC,   C%SNOWC)
    call cmp2d(s, "QSNOWXY", A%QSNOWXY, C%QSNOWXY)
    call cmp2d(s, "QRAINXY", A%QRAINXY, C%QRAINXY)
    ! canopy / vegetation
    call cmp2d(s, "TVXY",    A%TVXY,    C%TVXY)
    call cmp2d(s, "TGXY",    A%TGXY,    C%TGXY)
    call cmp2d(s, "CANWAT",  A%CANWAT,  C%CANWAT)
    call cmp2d(s, "CANLIQXY",A%CANLIQXY,C%CANLIQXY)
    call cmp2d(s, "CANICEXY",A%CANICEXY,C%CANICEXY)
    call cmp2d(s, "EAHXY",   A%EAHXY,   C%EAHXY)
    call cmp2d(s, "TAHXY",   A%TAHXY,   C%TAHXY)
    call cmp2d(s, "FWETXY",  A%FWETXY,  C%FWETXY)
    call cmp2d(s, "QSFC",    A%QSFC,    C%QSFC)
    call cmp2d(s, "LAI",     A%LAI,     C%LAI)
    call cmp2d(s, "XSAIXY",  A%XSAIXY,  C%XSAIXY)
    ! exchange coefficients (carried between steps)
    call cmp2d(s, "CMXY",    A%CMXY,    C%CMXY)
    call cmp2d(s, "CHXY",    A%CHXY,    C%CHXY)
    ! albedo history
    call cmp2d(s, "SNEQVOXY",A%SNEQVOXY,C%SNEQVOXY)
    call cmp2d(s, "ALBOLDXY",A%ALBOLDXY,C%ALBOLDXY)
    call cmp2d(s, "TAUSSXY", A%TAUSSXY, C%TAUSSXY)
    ! aquifer / groundwater
    call cmp2d(s, "ZWTXY",     A%ZWTXY,     C%ZWTXY)
    call cmp2d(s, "WAXY",      A%WAXY,      C%WAXY)
    call cmp2d(s, "WTXY",      A%WTXY,      C%WTXY)
    call cmp2d(s, "SMCWTDXY",  A%SMCWTDXY,  C%SMCWTDXY)
    call cmp2d(s, "DEEPRECHXY",A%DEEPRECHXY,C%DEEPRECHXY)
    call cmp2d(s, "RECHXY",    A%RECHXY,    C%RECHXY)
    ! accumulators / integrated state
    call cmp2d(s, "SFCRUNOFF", A%SFCRUNOFF, C%SFCRUNOFF)
    call cmp2d(s, "UDRUNOFF",  A%UDRUNOFF,  C%UDRUNOFF)
    call cmp2d(s, "SMSTAV",    A%SMSTAV,    C%SMSTAV)
    call cmp2d(s, "SMSTOT",    A%SMSTOT,    C%SMSTOT)
    call cmp2d(s, "ACSNOM",    A%ACSNOM,    C%ACSNOM)
    call cmp2d(s, "ACSNOW",    A%ACSNOW,    C%ACSNOW)
    ! dynamic-vegetation carbon pools (always allocated; evolve only under DVEG 2/5/6)
    call cmp2d(s, "LFMASSXY", A%LFMASSXY, C%LFMASSXY)
    call cmp2d(s, "RTMASSXY", A%RTMASSXY, C%RTMASSXY)
    call cmp2d(s, "STMASSXY", A%STMASSXY, C%STMASSXY)
    call cmp2d(s, "WOODXY",   A%WOODXY,   C%WOODXY)
    call cmp2d(s, "GDDXY",    A%GDDXY,    C%GDDXY)
  end subroutine run_scenario

  ! ---- scenario table ---------------------------------------------------------
  function configure(sel) result(s)
    character(*), intent(in) :: sel
    type(scenario_t) :: s
    ! defaults (the "land" baseline); scenarios override selectively.
    s%name = sel
    s%nsteps = 3
    s%tair = 290.0_c_kind_noahmp; s%qv = 0.006_c_kind_noahmp
    s%swdown = 400.0_c_kind_noahmp; s%glw = 340.0_c_kind_noahmp
    s%coszen = 0.6_c_kind_noahmp
    s%rainbl = 0.0_c_kind_noahmp; s%sr = 0.0_c_kind_noahmp
    s%ivgtyp = -1; s%iopt_dveg = -1; s%iopt_runsub = -1
    s%signature = "TSLB"

    select case (sel)
    case ("land")
       ! defaults
    case ("rain")
       s%nsteps = 4
       s%tair = 291.0_c_kind_noahmp; s%qv = 0.010_c_kind_noahmp
       s%swdown = 300.0_c_kind_noahmp
       s%rainbl = 8.0_c_kind_noahmp; s%sr = 0.0_c_kind_noahmp
       s%signature = "SMOIS"
    case ("snow")
       s%nsteps = 8
       s%tair = 263.0_c_kind_noahmp; s%qv = 0.002_c_kind_noahmp
       s%swdown = 150.0_c_kind_noahmp; s%glw = 250.0_c_kind_noahmp
       s%coszen = 0.3_c_kind_noahmp
       s%rainbl = 6.0_c_kind_noahmp; s%sr = 1.0_c_kind_noahmp
       s%signature = "SNOW"
    case ("carbon")
       s%nsteps = 6
       s%tair = 298.0_c_kind_noahmp; s%qv = 0.012_c_kind_noahmp
       s%swdown = 600.0_c_kind_noahmp; s%glw = 360.0_c_kind_noahmp
       s%coszen = 0.8_c_kind_noahmp
       s%iopt_dveg = 2                                 ! dynamic veg + carbon on
       s%signature = "LFMASS"
    case ("glacier")
       s%nsteps = 8
       s%tair = 263.0_c_kind_noahmp; s%qv = 0.002_c_kind_noahmp
       s%swdown = 150.0_c_kind_noahmp; s%glw = 250.0_c_kind_noahmp
       s%coszen = 0.3_c_kind_noahmp
       s%rainbl = 6.0_c_kind_noahmp; s%sr = 1.0_c_kind_noahmp
       s%ivgtyp = 15                                   ! ISICE (glacier) per fixture ISICE att
       s%signature = "SNOW"
    case ("groundwater")
       s%nsteps = 6
       s%tair = 291.0_c_kind_noahmp; s%qv = 0.010_c_kind_noahmp
       s%swdown = 300.0_c_kind_noahmp
       s%rainbl = 8.0_c_kind_noahmp; s%sr = 0.0_c_kind_noahmp
       s%iopt_runsub = 1                               ! TOPMODEL + groundwater (SIMGM)
       s%signature = "ZWT"
    case default
       write(0,'(A)') "test_io_restart_step: unknown scenario '"//trim(sel)//"'"
       call exit(2)
    end select
  end function configure

  ! ---- the subsystem must have actually evolved in leg A ----------------------
  subroutine require_active(s, blk, tslb_i, smois_i, snow_i, lf_i, zwt_i)
    type(scenario_t),       intent(in) :: s
    type(NoahmpIO_type),    intent(in) :: blk
    real(kind=kind_noahmp), intent(in) :: tslb_i(:,:,:), smois_i(:,:)
    real(kind=kind_noahmp), intent(in) :: snow_i(:,:), lf_i(:,:), zwt_i(:,:)
    logical :: moved
    select case (trim(s%signature))
    case ("TSLB");   moved = any(blk%TSLB       /= tslb_i)
    case ("SMOIS");  moved = any(blk%SMOIS(:,:,1) /= smois_i)
    case ("SNOW");   moved = any(blk%SNOW       /= snow_i)
    case ("LFMASS"); moved = any(blk%LFMASSXY   /= lf_i)
    case ("ZWT");    moved = any(blk%ZWTXY      /= zwt_i)
    case default;    moved = .false.
    end select
    call tio_expect(moved, "scenario exercised subsystem (signature '"// &
                    trim(s%signature)//"' evolved in leg A)")
  end subroutine require_active

  ! allocate + set an integer(C_INT) coupled pointer scalar
  subroutine seti(p, v)
    integer(C_INT), pointer, intent(out) :: p
    integer,                 intent(in)  :: v
    allocate(p); p = v
  end subroutine seti

  ! allocate + set a real coupled pointer scalar
  subroutine setr(p, v)
    real(kind=c_kind_noahmp), pointer, intent(out) :: p
    real(kind=c_kind_noahmp),          intent(in)  :: v
    allocate(p); p = v
  end subroutine setr

  ! Full ERF per-block cold-init chain (mirrors test_io_driver / ERF_NOAHMP_Init),
  ! plus the scenario's option / surface-type overrides at the correct injection
  ! points (options before InitMain; IVGTYP after the land read, before InitMain).
  subroutine init_block(blk, fpath, bcomm, s)
    type(NoahmpIO_type), intent(inout) :: blk
    character(*),        intent(in)    :: fpath
    integer,             intent(in)    :: bcomm
    type(scenario_t),    intent(in)    :: s

    call seti(blk%XSTART, 0);      call seti(blk%XEND, -1)
    call seti(blk%YSTART, 0);      call seti(blk%YEND, -1)
    call seti(blk%IDS, 0);         call seti(blk%IDE, NX-1)
    call seti(blk%JDS, 0);         call seti(blk%JDE, NY-1)
    call seti(blk%KDS, 1);         call seti(blk%KDE, 2)
    call seti(blk%ITS, 0);         call seti(blk%ITE, NX-1)
    call seti(blk%JTS, 0);         call seti(blk%JTE, NY-1)
    call seti(blk%KTS, 1);         call seti(blk%KTE, 2)
    call seti(blk%IMS, 0);         call seti(blk%IME, NX-1)
    call seti(blk%JMS, 0);         call seti(blk%JME, NY-1)
    call seti(blk%KMS, 1);         call seti(blk%KME, 2)
    call seti(blk%NSOIL, NSOIL);   call seti(blk%NSNOW, 4)
    call seti(blk%NUMRAD, 2);      call seti(blk%RANK, 0)
    call seti(blk%BLKID, 0);       call seti(blk%LEVEL, 0)
    call seti(blk%COMM, bcomm)
    allocate(blk%ITIMESTEP)
    call seti(blk%NTIME, 0)
    call setr(blk%DTBL, 3600.0_c_kind_noahmp)
    call setr(blk%ZLVL, 10.0_c_kind_noahmp)

    blk%YR     = 2023
    blk%JULIAN = 229.0_kind_noahmp

    call NoahmpReadNamelist(blk)

    ! Physics-option overrides (safe post-namelist; carbon/aquifer arrays are
    ! always allocated so no re-allocation is needed).
    if (s%iopt_dveg   >= 0) blk%IOPT_DVEG   = s%iopt_dveg
    if (s%iopt_runsub >= 0) blk%IOPT_RUNSUB = s%iopt_runsub

    blk%erf_setup_file_lev = trim(fpath)
    blk%erf_setup_file_01  = trim(fpath)
    call NoahmpReadLandHeader(blk)

    ! Real tile/domain/memory bounds span the whole block (set AFTER the header
    ! read, exactly as ERF does; the 0/-1 placeholders would size arrays to zero).
    blk%XSTART = 0; blk%XEND = NX-1; blk%YSTART = 0; blk%YEND = NY-1
    blk%IDS = 0; blk%IDE = NX-1; blk%JDS = 0; blk%JDE = NY-1
    blk%ITS = 0; blk%ITE = NX-1; blk%JTS = 0; blk%JTE = NY-1
    blk%IMS = 0; blk%IME = NX-1; blk%JMS = 0; blk%JME = NY-1

    call NoahmpIOVarInitDefault(blk)
    call NoahmpReadTable(blk)
    call NoahmpReadLandMain(blk)
    ! Surface-type override (e.g. glacier) after the land read, before init uses it.
    if (s%ivgtyp >= 0) blk%IVGTYP = s%ivgtyp
    call NoahmpInitMain(blk)
  end subroutine init_block

  ! Re-apply the scenario's forcing, then advance one step. Re-applying every step
  ! keeps the forcing seen at step N identical across legs even if the driver
  ! mutates the forcing arrays in place.
  subroutine step(blk, it_in, s)
    type(NoahmpIO_type), intent(inout) :: blk
    integer,             intent(in)    :: it_in
    type(scenario_t),    intent(in)    :: s
    blk%T_PHY(:,1,:)   = s%tair
    blk%QV_CURR(:,1,:) = s%qv
    blk%U_PHY(:,1,:)   = 3.0_c_kind_noahmp
    blk%V_PHY(:,1,:)   = 1.0_c_kind_noahmp
    blk%P8W(:,1,:)     = 1.0e5_c_kind_noahmp
    blk%SWDOWN         = s%swdown
    blk%GLW            = s%glw
    blk%COSZEN         = s%coszen
    ! precip: RAINBL is the total; SR routes it to snow (1) or rain (0).
    blk%RAINBL         = s%rainbl
    blk%SR             = s%sr
    blk%MP_RAINNC      = 0.0_c_kind_noahmp
    blk%MP_SNOW        = 0.0_c_kind_noahmp
    blk%MP_GRAUP       = 0.0_c_kind_noahmp
    blk%MP_HAIL        = 0.0_c_kind_noahmp
    blk%ITIMESTEP      = it_in
    call NoahmpDriverMain(blk)
  end subroutine step

  ! Assumed-shape so negative-indexed snow arrays re-index 1..size identically.
  ! The scenario name qualifies every check so failures are unambiguous in the
  ! merged run.
  subroutine cmp2d(s, name, a2, b2)
    type(scenario_t),       intent(in) :: s
    character(*),           intent(in) :: name
    real(kind=kind_noahmp), intent(in) :: a2(:,:), b2(:,:)
    real(kind=kind_noahmp) :: d
    d = maxval(abs(a2 - b2))
    write(*,'(A,A10,A,ES13.6)') "  cmp ", name, " max|diff| = ", d
    call tio_expect(d == 0.0_kind_noahmp, &
                    trim(s%name)//":"//trim(name)//" bit-exact across restart")
  end subroutine cmp2d

  subroutine cmp3d(s, name, a3, b3)
    type(scenario_t),       intent(in) :: s
    character(*),           intent(in) :: name
    real(kind=kind_noahmp), intent(in) :: a3(:,:,:), b3(:,:,:)
    real(kind=kind_noahmp) :: d
    d = maxval(abs(a3 - b3))
    write(*,'(A,A10,A,ES13.6)') "  cmp ", name, " max|diff| = ", d
    call tio_expect(d == 0.0_kind_noahmp, &
                    trim(s%name)//":"//trim(name)//" bit-exact across restart")
  end subroutine cmp3d

  subroutine cmp2di(s, name, a2, b2)
    type(scenario_t), intent(in) :: s
    character(*),     intent(in) :: name
    integer,          intent(in) :: a2(:,:), b2(:,:)
    write(*,'(A,A10,A,I0)') "  cmp ", name, " max|diff| = ", maxval(abs(a2 - b2))
    call tio_expect(all(a2 == b2), &
                    trim(s%name)//":"//trim(name)//" bit-exact across restart")
  end subroutine cmp2di

end program test_io_restart_step
