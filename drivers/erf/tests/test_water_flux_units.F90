! ===========================================================================
! test_water_flux_units -- water output-transfer units/accumulation contract.
!
! Guards the units contract questioned in the "exposed variables" bug report.
! By the time WaterVarOutTransfer runs, SoilWaterMain has already multiplied
! the runoff/tile-drain rates by SoilTimeStep (src/SoilWaterMainMod.F90:262-264),
! so the flux members it reads are DEPTHS [mm per soil timestep], not rates.
! This test drives the real WaterVarOutTransfer and pins that behavior so a
! well-meaning "fix" (e.g. multiplying QTDRAIN by DTBL, or dividing the runoff
! back out) is caught:
!
!   * land soil point: QTDRAIN / SFCRUNOFF / UDRUNOFF ACCUMULATE the depth as-is
!     (two calls of depth d give 2*d, NOT d*dt and not d), while RUNSFXY /
!     RUNSBXY are per-call SNAPSHOTS (stay d, do not accumulate).
!   * glacier point: the sentinel zeroes TileDrain (so QTDRAIN never grows) and
!     scales the runoff rate by MainTimeStep before it is written, matching the
!     "(glacier points: per main timestep)" interface note.
!
! Pure Fortran, no NetCDF, no domain fixture -- fits the Tier-1 test style.
! Exit code 0 = pass.
! ===========================================================================
program test_water_flux_units

  use iso_c_binding,      only : C_INT
  use Machine,            only : c_kind_noahmp
  use NoahmpIO_fi,        only : NoahmpIO_vect, NoahmpIOTypeVectInit_fi
  use NoahmpIOVarType,    only : NoahmpIO_type
  use NoahmpIOVarInitMod, only : NoahmpIOVarInitDefault
  use NoahmpVarType,      only : noahmp_type
  use WaterVarInitMod,    only : WaterVarInitDefault
  use WaterVarOutTransferMod, only : WaterVarOutTransfer

  implicit none

  integer(C_INT) :: level, nblocks
  type(NoahmpIO_type), pointer :: blk
  type(noahmp_type)            :: noahmp
  integer :: nfail
  integer, parameter :: XS=1, XE=4, YS=1, YE=3, KMS=1, KME=3, NSOIL=4, NSNOW=3, NUMRAD=2

  ! Known flux values already reduced to depth [mm per soil timestep] on the land
  ! path, or still a rate [mm/s] on the glacier path (scaled by MainTimeStep).
  real(kind=c_kind_noahmp), parameter :: TD  = 2.0_c_kind_noahmp   ! tile-drain depth
  real(kind=c_kind_noahmp), parameter :: RS  = 3.0_c_kind_noahmp   ! surface runoff depth
  real(kind=c_kind_noahmp), parameter :: SS  = 5.0_c_kind_noahmp   ! subsurface runoff depth
  real(kind=c_kind_noahmp), parameter :: DT  = 10.0_c_kind_noahmp  ! MainTimeStep / DTBL [s]

  nfail   = 0
  level   = 0
  nblocks = 1

  ! ---- size + allocate the NoahmpIO block (2-D/3-D output storage) ----
  call NoahmpIOTypeVectInit_fi(level, nblocks)
  blk => NoahmpIO_vect(level)%NoahmpIO(0)

  ! Coupled dimension scalars are C++-owned pointer components; in this pure
  ! Fortran path give each its own target before use (see test_fortran_alloc).
  allocate(blk%XSTART); blk%XSTART = XS
  allocate(blk%XEND);   blk%XEND   = XE
  allocate(blk%YSTART); blk%YSTART = YS
  allocate(blk%YEND);   blk%YEND   = YE
  allocate(blk%KMS);    blk%KMS    = KMS
  allocate(blk%KME);    blk%KME    = KME
  allocate(blk%NSOIL);  blk%NSOIL  = NSOIL
  allocate(blk%NSNOW);  blk%NSNOW  = NSNOW
  allocate(blk%NUMRAD); blk%NUMRAD = NUMRAD
  allocate(blk%ITIMESTEP)               ! initialized by NoahmpIOVarInitDefault
  allocate(blk%DTBL);   blk%DTBL   = DT ! used by ACSNOM / IRELOSS accumulations
  blk%IOPT_SOIL = 1; blk%IOPT_ALB = 1; blk%IOPT_WETLAND = 0; blk%SF_URBAN_PHYSICS = 0
  call NoahmpIOVarInitDefault(blk)      ! zeroes QTDRAIN/SFCRUNOFF/UDRUNOFF

  ! ---- build a minimal noahmp column and allocate its water arrays ----
  noahmp%config%domain%NumSoilLayer            = NSOIL
  noahmp%config%domain%NumSnowLayerMax         = NSNOW
  noahmp%config%domain%NumDensitySnwAgeSnicar  = 1
  noahmp%config%domain%NumTempGradSnwAgeSnicar = 1
  noahmp%config%domain%NumTempSnwAgeSnicar     = 1
  noahmp%config%domain%MainTimeStep            = DT
  noahmp%config%nmlist%OptSnowAlbedo           = 1   ! /= 3 : skip SNICAR block
  noahmp%config%nmlist%OptWetlandModel         = 0
  call WaterVarInitDefault(noahmp)

  ! Neutralize the fluxes that feed unrelated accumulations we do not assert,
  ! so no undefined_real sentinel propagates into ACSNOW/ACSNOM/IR* products.
  noahmp%water%flux%MeltGroundSnow          = 0.0_c_kind_noahmp
  noahmp%water%flux%IrrigationRateSprinkler = 0.0_c_kind_noahmp
  noahmp%water%flux%IrrigationRateMicro     = 0.0_c_kind_noahmp
  noahmp%water%flux%IrrigationRateFlood     = 0.0_c_kind_noahmp
  noahmp%water%flux%EvapIrriSprinkler       = 0.0_c_kind_noahmp
  noahmp%water%state%FrozenPrecipFrac       = 0.0_c_kind_noahmp

  ! =====================================================================
  ! Phase A -- land soil point (IndicatorIceSfc == 0): depth accumulation
  ! =====================================================================
  noahmp%config%domain%GridIndexI     = 1
  noahmp%config%domain%GridIndexJ     = 1
  noahmp%config%domain%IndicatorIceSfc = 0

  call set_land_fluxes()
  call WaterVarOutTransfer(noahmp, blk)

  ! after 1 call: accumulators hold one depth, snapshots hold the depth
  call expect_close(blk%QTDRAIN (1,1), TD, "land QTDRAIN after 1 call")
  call expect_close(blk%SFCRUNOFF(1,1), RS, "land SFCRUNOFF after 1 call")
  call expect_close(blk%UDRUNOFF (1,1), SS, "land UDRUNOFF after 1 call")
  call expect_close(blk%RUNSFXY  (1,1), RS, "land RUNSFXY after 1 call")
  call expect_close(blk%RUNSBXY  (1,1), SS, "land RUNSBXY after 1 call")

  call set_land_fluxes()
  call WaterVarOutTransfer(noahmp, blk)

  ! after 2 calls: QTDRAIN/SFCRUNOFF/UDRUNOFF doubled (depth accumulation).
  ! Critically this is 2*TD, NOT TD*DTBL -- proves a rate is not summed and the
  ! depth is not multiplied by a timestep on the way in.
  call expect_close(blk%QTDRAIN (1,1), 2.0_c_kind_noahmp*TD, "land QTDRAIN accumulates depth")
  call expect_close(blk%SFCRUNOFF(1,1), 2.0_c_kind_noahmp*RS, "land SFCRUNOFF accumulates")
  call expect_close(blk%UDRUNOFF (1,1), 2.0_c_kind_noahmp*SS, "land UDRUNOFF accumulates")
  ! RUNSFXY/RUNSBXY are per-call snapshots: still one depth, not doubled.
  call expect_close(blk%RUNSFXY  (1,1), RS, "land RUNSFXY is a snapshot")
  call expect_close(blk%RUNSBXY  (1,1), SS, "land RUNSBXY is a snapshot")

  ! =====================================================================
  ! Phase B -- glacier point (IndicatorIceSfc == -1): sentinel + scaling
  ! =====================================================================
  noahmp%config%domain%GridIndexI     = 2
  noahmp%config%domain%GridIndexJ     = 1
  noahmp%config%domain%IndicatorIceSfc = -1

  ! Here the fluxes are still rates [mm/s]; the glacier branch scales runoff by
  ! MainTimeStep and forces tile drainage to zero.
  noahmp%water%flux%TileDrain        = TD
  noahmp%water%flux%RunoffSurface    = RS
  noahmp%water%flux%RunoffSubsurface = SS
  call WaterVarOutTransfer(noahmp, blk)

  ! glacier has no tile drainage: QTDRAIN(2,1) stays 0 (was zeroed at init)
  call expect_close(blk%QTDRAIN(2,1), 0.0_c_kind_noahmp, "glacier QTDRAIN stays zero")
  ! runoff diagnostics scaled to depth per main timestep
  call expect_close(blk%RUNSFXY(2,1), RS*DT, "glacier RUNSFXY scaled by MainTimeStep")
  call expect_close(blk%RUNSBXY(2,1), SS*DT, "glacier RUNSBXY scaled by MainTimeStep")

  if (nfail == 0) then
     write(*,'(A)') "PASS: test_water_flux_units"
     call exit(0)
  else
     write(*,'(A,I0,A)') "FAILED: test_water_flux_units (", nfail, " check(s))"
     call exit(1)
  end if

contains

  ! Land path: SoilWaterMain has already reduced rates to depths, so the flux
  ! members are the depths themselves. WaterVarOutTransfer does not mutate them
  ! on the land path, but re-set before each call to be explicit.
  subroutine set_land_fluxes()
    noahmp%water%flux%TileDrain        = TD
    noahmp%water%flux%RunoffSurface    = RS
    noahmp%water%flux%RunoffSubsurface = SS
  end subroutine set_land_fluxes

  subroutine expect_close(got, want, name)
    real(kind=c_kind_noahmp), intent(in) :: got, want
    character(*), intent(in) :: name
    real(kind=c_kind_noahmp) :: tol
    tol = 1.0e-4_c_kind_noahmp * max(1.0_c_kind_noahmp, abs(want))
    if (abs(got - want) > tol) then
       write(0,'(A,A,A,ES14.6,A,ES14.6)') "  FAIL: ", name, " got=", got, " want=", want
       nfail = nfail + 1
    end if
  end subroutine expect_close

end program test_water_flux_units
