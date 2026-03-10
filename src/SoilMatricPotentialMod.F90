module SoilMatricPotentialMod

!!! compute soil matric potential that will be used for
!!! soil moisture transpiration factor and dynamic root calculations

  use Machine
  use NoahmpVarType

  implicit none

contains

  subroutine SoilMatricPotential(noahmp)

! ------------------------ Code history -----------------------------------
! Original Noah-MP subroutine: None (embedded in ENERGY subroutine)
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! -------------------------------------------------------------------------

    implicit none

! in & out variables
    type(noahmp_type), intent(inout) :: noahmp

! local variables
    integer                          :: IndSoil       ! loop index

! --------------------------------------------------------------------
    associate(                                                                             &
              SurfaceType               => noahmp%config%domain%SurfaceType               ,& ! in,  surface type 1-soil; 2-lake
              NumSoilLayer              => noahmp%config%domain%NumSoilLayer              ,& ! in,  number of soil layers
              SoilMatPotentialWilt      => noahmp%water%param%SoilMatPotentialWilt        ,& ! in,  soil metric potential for wilting point [m]
              SoilMatPotentialSat       => noahmp%water%param%SoilMatPotentialSat         ,& ! in,  saturated soil matric potential [m]
              SoilMoistureSat           => noahmp%water%param%SoilMoistureSat             ,& ! in,  saturated value of soil moisture [m3/m3]
              SoilExpCoeffB             => noahmp%water%param%SoilExpCoeffB               ,& ! in,  soil B parameter
              SoilLiqWater              => noahmp%water%state%SoilLiqWater                ,& ! in,  soil water content [m3/m3]
              SoilMatPotential          => noahmp%water%state%SoilMatPotential             & ! out, soil matrix potential [m]
             )
! ----------------------------------------------------------------------

    ! only for soil point
    if ( SurfaceType ==1 ) then

       do IndSoil = 1, NumSoilLayer
             SoilMatPotential(IndSoil) = max(SoilMatPotentialWilt, -SoilMatPotentialSat(IndSoil) * &
                                            (max(0.01,SoilLiqWater(IndSoil))/SoilMoistureSat(IndSoil)) ** &
                                            (-SoilExpCoeffB(IndSoil)))
       enddo

    endif

    end associate

  end subroutine SoilMatricPotential

end module SoilMatricPotentialMod

