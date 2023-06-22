module SnowAlbedoSnicarGlacierMod

!!! Compute snow albedo based on SNICAR scheme (Flanner et al. (2021) GMD)

  use Machine
  use NoahmpVarType
  use ConstantDefineMod
  use SnowFreshRadiusMod,     only : SnowFreshRadius
  use SnowAgeSnicarMod,       only : SnowAgeSnicar
  use SnowRadiationSnicarMod, only : SnowRadiationSnicar

  implicit none

contains

  subroutine SnowAlbedoSnicarGlacier(noahmp)

! ------------------------ Code history -----------------------------------
! Original code: Cenlin He and CTSM team 
! -------------------------------------------------------------------------

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

! local variable
    integer                          :: flg_slr_in          ! flag: 1 for direct-beam incident flux, 2 for diffuse incident flux
    integer                          :: flg_snw_ice         ! flag: 1 is land case, 2 is seaice case
! --------------------------------------------------------------------
    associate(                                                                 &
              NumSwRadBand        => noahmp%config%domain%NumSwRadBand        ,& ! in,  number of solar radiation wave bands
              AlbedoSnowDir       => noahmp%energy%state%AlbedoSnowDir        ,& ! out, snow albedo for direct(1=vis, 2=nir)
              AlbedoSnowDif       => noahmp%energy%state%AlbedoSnowDif        ,& ! out, snow albedo for diffuse(1=vis, 2=nir)
              FracRadSwAbsSnowDir => noahmp%energy%flux%FracRadSwAbsSnowDir   ,& ! out, direct solar flux factor absorbed by snow [frc] (-NumSnowLayerMax+1:1,NumSwRadBand)
              FracRadSwAbsSnowDif => noahmp%energy%flux%FracRadSwAbsSnowDif    & ! out, diffuse solar flux factor absorbed by snow [frc] (-NumSnowLayerMax+1:1,NumSwRadBand)

             )
! ----------------------------------------------------------------------

    ! initialization
    AlbedoSnowDir(1:NumSwRadBand) = 0.0
    AlbedoSnowDif(1:NumSwRadBand) = 0.0

    FracRadSwAbsSnowDir(:,1:NumSwRadBand) = 0.0
    FracRadSwAbsSnowDif(:,1:NumSwRadBand) = 0.0

    ! snow radius
    call SnowFreshRadius(noahmp)
    call SnowAgeSnicar(noahmp)

    flg_snw_ice = 2 !Seaice

    flg_slr_in = 1 !Direct
    call SnowRadiationSnicar(noahmp,flg_slr_in,flg_snw_ice) 

    flg_slr_in = 2 !Diffuse
    call SnowRadiationSnicar(noahmp,flg_slr_in,flg_snw_ice)

    end associate

  end subroutine SnowAlbedoSnicarGlacier

end module SnowAlbedoSnicarGlacierMod
