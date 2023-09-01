module SnowAlbedoSnicarMod

!!! Compute snow albedo based on SNICAR scheme (Flanner et al. (2021) GMD)

  use Machine
  use NoahmpVarType
  use ConstantDefineMod
  use SnowRadiationSnicarMod, only : SnowRadiationSnicar

  implicit none

contains

  subroutine SnowAlbedoSnicar(noahmp)

! ------------------------ Code history -----------------------------------
! code: T.-S. Lin, C. He, et al. (2023)
! -------------------------------------------------------------------------

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

! local variable
    integer                          :: flg_slr_in          ! flag: 1 for direct-beam incident flux, 2 for diffuse incident flux
! --------------------------------------------------------------------
    associate(                                                                 &
              NumSwRadBand        => noahmp%config%domain%NumSwRadBand        ,& ! in,  number of solar radiation wave bands
              AlbedoSnowDir       => noahmp%energy%state%AlbedoSnowDir        ,& ! out, snow albedo for direct(1=vis, 2=nir)
              AlbedoSnowDif       => noahmp%energy%state%AlbedoSnowDif         & ! out, snow albedo for diffuse(1=vis, 2=nir)

             )
! ----------------------------------------------------------------------

    ! initialization
    AlbedoSnowDir(1:NumSwRadBand) = 0.0
    AlbedoSnowDif(1:NumSwRadBand) = 0.0

    flg_slr_in = 1 !Direct
    call SnowRadiationSnicar(noahmp,flg_slr_in) 

    flg_slr_in = 2 !Diffuse
    call SnowRadiationSnicar(noahmp,flg_slr_in)
    
    end associate

  end subroutine SnowAlbedoSnicar

end module SnowAlbedoSnicarMod
