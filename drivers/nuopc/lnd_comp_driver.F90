module lnd_comp_driver

  ! This file contains the NoahMP land surface model driver

  use ESMF            , only: operator(+), operator(-), operator(/)
  use ESMF            , only: ESMF_LogWrite, ESMF_SUCCESS, ESMF_LOGMSG_INFO
  use ESMF            , only: ESMF_Finalize, ESMF_END_ABORT
  use ESMF            , only: ESMF_VM, ESMF_VMCommWaitAll, ESMF_VMGet
  use ESMF            , only: ESMF_GridComp, ESMF_GridCompGet
  use ESMF            , only: ESMF_Clock, ESMF_ClockGet, ESMF_TimeIsLeapYear
  use ESMF            , only: ESMF_Time, ESMF_TimeGet, ESMF_TimeSet
  use ESMF            , only: ESMF_TimeInterval, ESMF_TimeIntervalGet
  use ESMF            , only: ESMF_ClockGetAlarm, ESMF_AlarmIsCreated
  use ESMF            , only: ESMF_AlarmIsRinging, ESMF_Alarm, ESMF_AlarmRingerOff

  use lnd_comp_kind   , only: r8 => shr_kind_r8
  use lnd_comp_kind   , only: cl => shr_kind_cl
  use lnd_comp_types  , only: noahmp_type
  use lnd_comp_types  , only: fldsToLnd, fldsToLnd_num
  use lnd_comp_shr    , only: chkerr
  use lnd_comp_io     , only: read_static, read_initial, read_restart
  use lnd_comp_io     , only: write_mosaic_output
  use lnd_comp_import_export, only: check_for_connected

  use sfc_diff        , only: stability
  use noahmpdrv       , only: noahmpdrv_init, noahmpdrv_run
  use namelist_soilveg, only: z0_data
  use set_soilveg_mod , only: set_soilveg
  use noahmp_tables   , only: bexp_table, smcmax_table, smcwlt_table
  use noahmp_tables   , only: dwsat_table, dksat_table, psisat_table
  use noahmp_tables   , only: isurban_table
  use funcphys        , only: gpvs
  use physcons        , only: con_hvap, con_cp, con_jcal
  use physcons        , only: con_eps, con_epsm1, con_fvirt
  use physcons        , only: con_rd, con_hfus
  use physcons        , only: rhoh2o => rhowater
  use physcons        , only: p0 => con_p0
  use physcons        , only: cappa => con_rocp
  use physcons        , only: con_g
  use physcons        , only: karman

  implicit none
  private

  public :: drv_init, drv_run

  !--------------------------------------------------------------------------
  ! Private module data
  !--------------------------------------------------------------------------

  integer                :: dbug = 1
  type(ESMF_Time)        :: epoc
  type(ESMF_Time)        :: dummTime, prevTime, nextTime
  character(*),parameter :: modName =  "(lnd_comp_driver)"

  character(len=*) , parameter :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine drv_init(gcomp, noahmp, rc)

    ! input/output variables
    type(ESMF_GridComp), intent(in)    :: gcomp
    type(noahmp_type)  , intent(inout) :: noahmp 
    integer            , intent(out)   :: rc

    ! local variables
    integer            :: me
    integer, parameter :: nlunit = 9999
    integer, parameter :: lsm = 2
    integer, parameter :: lsm_noahmp = 2
    real(r8), pointer  :: pores(:) => null() !< max soil moisture for a given soil type for land surface model
    real(r8), pointer  :: resid(:) => null() !< min soil moisture for a given soil type for land surface model
    type(ESMF_VM)      :: vm
    character(len=*), parameter :: subname=trim(modName)//':(drv_init) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! ----------------------
    ! Query component
    ! ----------------------

    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_VMGet(vm, localPet=me, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------
    ! Set epoc as 1970-01-01 00:00:00
    !----------------------

    noahmp%model%reference_date = "1970-01-01 00:00:00"
    call ESMF_TimeSet(epoc, yy=1970, mm=1, dd=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! ----------------------
    ! Allocate and initialize model data
    ! ----------------------

    call noahmp%Initialize(noahmp%domain%begl, noahmp%domain%endl, &
      noahmp%static%km, noahmp%static%lsnowl)

    ! ---------------------
    ! Read static information
    ! ---------------------

    call read_static(noahmp, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! ---------------------
    ! Read initial condition / restart file
    ! ---------------------

    if (noahmp%nmlist%restart_run) then
       ! read restart file
       call read_restart(noahmp, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       ! read initial conditions
       call read_initial(noahmp, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    !----------------------
    ! initialize variables
    !----------------------

    noahmp%model%dry = .false.
    where(noahmp%domain%mask == 1) noahmp%model%dry = .true.
    noahmp%model%flag_iter = .true.
    noahmp%model%thsfc_loc = .true.
    noahmp%model%emiss(:) = noahmp%nmlist%initial_emiss
    noahmp%model%alb_monthly(:,:) = noahmp%nmlist%initial_albedo

    !----------------------
    ! initialize soil and vegetation
    !----------------------

    call noahmpdrv_init(lsm, lsm_noahmp, me, noahmp%static%isot, noahmp%static%ivegsrc, &
      nlunit, noahmp%model%pores, noahmp%model%resid, noahmp%static%do_mynnsfclay, &
      noahmp%static%do_mynnedmf, noahmp%static%errmsg, noahmp%static%errflg)
    call gpvs()

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine drv_init

  !===============================================================================
  subroutine drv_run(gcomp, noahmp, rc)

    ! input/output variables
    type(ESMF_GridComp), intent(in)    :: gcomp
    type(noahmp_type)  , intent(inout) :: noahmp 
    integer            , intent(out)   :: rc

    ! local variables
    logical, save               :: first_time = .true.
    integer                     :: i, is, step
    integer                     :: year, month, day, hour, minute, second
    real(r8)                    :: now_time
    character(len=cl)           :: filename
    logical                     :: restart_write
    type(ESMF_VM)               :: vm
    type(ESMF_Clock)            :: clock
    type(ESMF_Alarm)            :: alarm
    type(ESMF_Time)             :: startTime
    type(ESMF_Time)             :: currTime
    type(ESMF_TimeInterval)     :: timeStep
    real(r8)                    :: virtfac, tv1, thv1, tvs
    real(r8)                    :: tem1, tem2, z0max, ztmax, czilc
    real(r8)                    :: bexp, ddz, smcmax, smcwlt, dwsat, dksat, psisat
    real(r8), save, allocatable :: zs(:)
    real(r8), parameter         :: tfreeze = 2.7315e+2_r8 
    real(r8), parameter         :: zmin = 1.0e-6_r8
    real(r8), parameter         :: qmin = 1.0e-8_r8
    real(r8), parameter         :: z0lo = 0.1_r8
    real(r8), parameter         :: z0up = 1.0_r8
    real(r8), parameter         :: log01 = log(0.01_r8)
    real(r8), parameter         :: log05 = log(0.05_r8)
    real(r8), parameter         :: log07 = log(0.07_r8)
    character(len=*),parameter  :: subname = trim(modName)//':(drv_run) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    !----------------------
    ! Query clock and set timestep, current time etc. 
    !----------------------

    call ESMF_GridCompGet(gcomp, clock=clock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGet(clock, startTime=startTime, currTime=currTime, timeStep=timeStep, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_TimeIntervalGet(currTime-epoc, s_r8=now_time, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_TimeGet(currTime, yy=year, mm=month, dd=day, &
      h=hour, m=minute, s=second, dayOfYear_r8=noahmp%model%julian, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    noahmp%model%iyrlen = 365
    if (ESMF_TimeIsLeapYear(currTime, rc=rc)) then
       noahmp%model%iyrlen = 366
    end if

    if (first_time) then
       ! use coupling time step and internal time step of model
       call ESMF_TimeIntervalGet(timeStep, s_r8=noahmp%static%delt, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    !----------------------
    ! initialize model variables
    !----------------------

    step = int((currTime-startTime)/timeStep)
    if (.not. noahmp%nmlist%restart_run .and. step == 0) then
       ! transfer common initial conditions for all configurations
       do i = noahmp%domain%begl, noahmp%domain%endl
          if (noahmp%domain%mask(i) == 1) then
             noahmp%model%weasd(i)  = noahmp%init%snow_water_equivalent(i)
             noahmp%model%snwdph(i) = noahmp%init%snow_depth(i)
             noahmp%model%canopy(i) = noahmp%init%canopy_water(i)
             noahmp%model%tskin(i)  = noahmp%init%skin_temperature(i)
             noahmp%model%stc(i,:)  = noahmp%init%soil_temperature(i,:)
             noahmp%model%smc(i,:)  = noahmp%init%soil_moisture(i,:)
             noahmp%model%slc(i,:)  = noahmp%init%soil_liquid(i,:)
          end if
       end do

       ! transfer custom initial conditions based on selected configuration
       if (trim(noahmp%nmlist%ic_type) == 'sfc') then
          where(noahmp%domain%mask(:) == 1)
             noahmp%model%zorl(:) = noahmp%init%surface_roughness(:)
             noahmp%model%ustar1(:) = noahmp%init%friction_velocity(:)
          end where
       else if (trim(noahmp%nmlist%ic_type) == 'custom') then
          ! get initial value of zorl from pre-defined table
          where(noahmp%model%vegtype(:) > 0) noahmp%model%zorl(:) = z0_data(noahmp%model%vegtype(:))*100.0_r8
          ! additional unit conversion for datm configuration
          where(noahmp%domain%mask(:) > 1) noahmp%model%snwdph(:) = noahmp%model%snwdph(:)*1000.0_r8
       end if

       ! initialize model variables
       call noahmp%InitializeStates(noahmp%nmlist, noahmp%static, month)
    end if

    !----------------------
    ! set internal model variables from forcing
    !----------------------

    ! set forcing height
    noahmp%model%zf = noahmp%forc%hgt

    ! set net shortwave radiation
    if (check_for_connected(fldsToLnd, fldsToLnd_num, 'Faxa_swnet')) then
       noahmp%model%snet = noahmp%forc%snet
    end if
    ! overwrite it if it is explicity specified
    ! this also fixed all zero net shortwave provided by CDEPS
    if (noahmp%nmlist%calc_snet) then
       noahmp%model%snet = noahmp%forc%dswsfc*(1.0_r8-noahmp%model%sfalb)
    end if

    ! set surface temperature, it is same with skin temperature over land
    noahmp%model%tsurf(:) = noahmp%model%tskin(:)

    ! set air pressure at surface adjacent layer
    ! cdeps provides Sa_pbot but Sa_prsl is used for coupling with fv3
    if (check_for_connected(fldsToLnd, fldsToLnd_num, 'Sa_pbot') .or. &
        check_for_connected(fldsToLnd, fldsToLnd_num, 'Sa_prsl')) then
       noahmp%model%prsl1 = noahmp%forc%pbot
    else
       ! calculate it from surface pressure, height and temperature
       noahmp%model%prsl1 = noahmp%forc%ps*exp(-1.0_r8*noahmp%model%zf/29.25_r8/noahmp%forc%t1)
    end if

    ! set dimensionless Exner function at surface adjacent layer
    if (check_for_connected(fldsToLnd, fldsToLnd_num, 'Sa_exner')) then
       noahmp%model%prslk1 = noahmp%forc%prslk1
    else
       ! calculate it based on bottom pressure
       noahmp%model%prslk1 = (noahmp%forc%pbot(:)/p0)**cappa
       ! following is used by ufs-land-driver
       !noahmp%model%prslk1 = (exp(noahmp%model%zf/29.25_r8/noahmp%forc%t1))**(2.0_r8/7.0_r8)
    end if

    ! set dimensionless Exner function at the ground surface
    noahmp%model%prsik1 = (noahmp%forc%ps(:)/p0)**cappa

    ! set Exner function ratio between midlayer and interface
    ! following is used by ufs-land-driver
    !noahmp%model%prslki = (exp(noahmp%model%zf/29.25_r8/noahmp%forc%t1))**(2.0_r8/7.0_r8)
    noahmp%model%prslki(:) = noahmp%model%prsik1(:)/noahmp%model%prslk1(:)

    ! set wind forcing to model internal variables
    noahmp%model%u1(:) = noahmp%forc%u1(:)
    noahmp%model%v1(:) = noahmp%forc%v1(:)

    ! NOTE: CCPP has addtional adjustment in wind speed which could lead to minor difference
    ! The modification is done in GFS_surface_generic_pre.F90
    noahmp%forc%wind(:) = sqrt(noahmp%forc%u1(:)**2+noahmp%forc%v1(:)**2)
    noahmp%forc%wind(:) = max(noahmp%forc%wind(:), 1.0_r8)

    ! set snow precipitation
    if (check_for_connected(fldsToLnd, fldsToLnd_num, 'Faxa_snow')) then
       noahmp%model%snow_mp(:) = noahmp%forc%snow(:)
    else
       noahmp%model%snow_mp(:) = 0.0_r8
       if (check_for_connected(fldsToLnd, fldsToLnd_num, 'Faxa_snowc')) then
          noahmp%model%snow_mp(:) = noahmp%forc%snowc(:)
       end if
       if (check_for_connected(fldsToLnd, fldsToLnd_num, 'Faxa_snowl')) then
          noahmp%model%snow_mp(:) = noahmp%model%snow_mp(:)+noahmp%forc%snowl(:)
       end if
    end if

    ! set precipitation, convective and large scale
    if (check_for_connected(fldsToLnd, fldsToLnd_num, 'Faxa_rainc')) then
       noahmp%model%rainc_mp(:) = noahmp%forc%tprcpc(:)
    else
       noahmp%model%rainc_mp(:) = 0.0_r8
    end if
    if (check_for_connected(fldsToLnd, fldsToLnd_num, 'Faxa_rainl')) then
       noahmp%model%rainn_mp(:) = noahmp%forc%tprcpl(:)
    else
       noahmp%model%rainn_mp(:) = 0.0_r8
    end if
    if (check_for_connected(fldsToLnd, fldsToLnd_num, 'Faxa_rain')) then
       noahmp%model%rainn_mp(:) = noahmp%forc%tprcp(:)
    end if

    ! convert mm/s to m
    noahmp%model%rainn_mp(:) = noahmp%model%rainn_mp(:)*noahmp%static%delt/1000.0_r8
    noahmp%model%rainc_mp(:) = noahmp%model%rainc_mp(:)*noahmp%static%delt/1000.0_r8

    ! calculate total precipitation
    noahmp%model%tprcp(:) = noahmp%model%rainn_mp(:)+noahmp%model%rainc_mp(:)

    !----------------------
    ! interpolate monthly data, vegetation fraction and mean sfc diffuse sw albedo (NOT used)
    !----------------------

    if (check_for_connected(fldsToLnd, fldsToLnd_num, 'vfrac')) then
       where(noahmp%forc%vegfrac(:) < 0.01_r8)
          noahmp%model%sigmaf(:) = 0.01_r8
       else where
          noahmp%model%sigmaf(:) = noahmp%forc%vegfrac(:)
       end where
    else
       call interpolate_monthly(currTime, noahmp%static%im, &
          noahmp%model%gvf_monthly, noahmp%model%sigmaf, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if
 
    ! not used by the model but is set internally in the driver to albedo_total
    call interpolate_monthly(currTime, noahmp%static%im, &
       noahmp%model%alb_monthly, noahmp%model%sfalb, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------
    ! calculate solar zenith angle
    !----------------------

    ! TODO: this could be coupling field for active atm
    call calc_cosine_zenith(currTime, noahmp%static%im, &
      noahmp%domain%lats, noahmp%domain%lons, noahmp%model%xcoszin, &
      noahmp%model%julian, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! -----------------------
    ! custom calculations, unit conversion, variable modification etc.
    ! -----------------------

    ! set snow/rain flag for precipitation
    noahmp%model%srflag = 0.0_r8
    where(noahmp%forc%t1 < tfreeze) noahmp%model%srflag = 1.0_r8

    ! unit conversion for precipitation
    ! TODO: do we need for datm coupling?
    ! convert mm/s to m
    !noahmp%forc%tprcp(:) = noahmp%forc%tprcp(:)*noahmp%static%delt/1000.0_r8
    ! datm
    !noahmp%model%rainn_mp = 1000.0_r8*noahmp%forc%tprcp/noahmp%static%delt

    ! they are not defined as coupling fields but CCPP provides those fields
    noahmp%model%graupel_mp = 0.0_r8
    noahmp%model%ice_mp = 0.0_r8

    !----------------------
    ! calculate initial values, ported from GFS_phys_time_vary.fv3 
    !----------------------

    ! set soil layers
    if (.not. allocated(zs)) allocate(zs(noahmp%nmlist%num_soil_levels))
    zs = -1.0_r8*noahmp%nmlist%soil_level_nodes
 
    do i = noahmp%domain%begl, noahmp%domain%endl
       if (noahmp%domain%mask(i) == 1) then
          if (noahmp%model%soiltyp(i) /= 0) then
             bexp   = bexp_table(noahmp%model%soiltyp(i))
             smcmax = smcmax_table(noahmp%model%soiltyp(i))
             smcwlt = smcwlt_table(noahmp%model%soiltyp(i))
             dwsat  = dwsat_table(noahmp%model%soiltyp(i))
             dksat  = dksat_table(noahmp%model%soiltyp(i))
             psisat = -psisat_table(noahmp%model%soiltyp(i))
          endif

          if (noahmp%model%vegtype(i) == isurban_table) then
             smcmax = 0.45_r8
             smcwlt = 0.40_r8
          endif

          if ((bexp > 0.0_r8) .and. (smcmax > 0.0_r8) .and. (-psisat > 0.0_r8)) then
             do is = 1, noahmp%nmlist%num_soil_levels
                if (is == 1)then
                   ddz = -zs(is+1)*0.5_r8
                elseif (is < noahmp%nmlist%num_soil_levels) then
                   ddz = (zs(is-1)-zs(is+1))*0.5_r8
                else
                   ddz = zs(is-1)-zs(is)
                endif
              noahmp%model%smoiseq(i,is) = min(max(find_eq_smc(bexp, dwsat, dksat, ddz, smcmax),1.e-4_r8),smcmax*0.99_r8)
            enddo
          else ! bexp <= 0.0
            noahmp%model%smoiseq(i,:) = smcmax
          endif  
       
          noahmp%model%smcwtdxy(i) = smcmax
       end if
    end do

    !----------------------
    ! call stability
    !----------------------
    do i = noahmp%domain%begl, noahmp%domain%endl
       if (noahmp%domain%mask(i) == 1 .and. noahmp%model%flag_iter(i)) then
          ! set initial value for ztmax
          ztmax = 1.0_r8

          !virtual temperature in middle of lowest layer
          virtfac = 1.0_r8+con_fvirt*max(noahmp%forc%q1(i),qmin)
          tv1 = noahmp%forc%t1(i)*virtfac

          if (noahmp%model%thsfc_loc) then
             ! use local potential temperature
             thv1 = noahmp%forc%t1(i)*noahmp%model%prslki(i)*virtfac
             tvs  = 0.5_r8*(noahmp%model%tsurf(i)+noahmp%model%tskin(i))*virtfac
          else
             ! use potential temperature referenced to 1000 hPa
             thv1 = noahmp%forc%t1(i)/noahmp%model%prslk1(i)*virtfac
             tvs = 0.5_r8*(noahmp%model%tsurf(i)+noahmp%model%tskin(i))/noahmp%model%prsik1(i)*virtfac
          end if

         ! set initial value for zvfun
         noahmp%model%zvfun(i) = 0.0_r8

         ! set initial value for z0max
         z0max = max(zmin, min(0.01_r8 * noahmp%model%zorl(i), noahmp%model%zf(i)))

         tem1 = 1.0_r8-noahmp%model%shdmax(i)
         tem2 = tem1*tem1
         tem1 = 1.0_r8-tem2

         if (noahmp%static%ivegsrc == 1) then
            if (noahmp%model%vegtype(i) == 10) then
               z0max = exp(tem2*log01+tem1*log07)
            elseif (noahmp%model%vegtype(i) == 6) then
               z0max = exp(tem2*log01+tem1*log05)
            elseif (noahmp%model%vegtype(i) == 7) then
               z0max = 0.01_r8
            elseif (noahmp%model%vegtype(i) == 16) then
               z0max = 0.01_r8
            else
               z0max = exp(tem2*log01+tem1*log(z0max))
            endif
         elseif (noahmp%static%ivegsrc == 2) then
            if (noahmp%model%vegtype(i) == 7) then
              z0max = exp(tem2*log01+tem1*log07)
            elseif (noahmp%model%vegtype(i) == 8) then
              z0max = exp(tem2*log01+tem1*log05)
            elseif (noahmp%model%vegtype(i) == 9) then
              z0max = 0.01_r8
            elseif (noahmp%model%vegtype(i) == 11) then
              z0max = 0.01_r8
            else
              z0max = exp(tem2*log01+tem1*log(z0max))
            endif
         endif

         ! TODO: add surface perturbations to z0max over land
         ! z0pert is defined in CCPP GFS_surface_generic_pre.F90
         ! it is all zero for control_p8 configuration

         ! limit z0max
         z0max = max(z0max, zmin)

         czilc = 10.0_r8**(-4.0_r8*z0max) ! Trier et al. (2011,WAF)
         czilc = max(min(czilc, 0.8_r8), 0.08_r8)
         tem1 = 1.0_r8-noahmp%model%sigmaf(i)
         czilc = czilc*tem1*tem1
         ztmax = z0max * exp( - czilc * karman * 258.2_r8 * sqrt(noahmp%model%ustar1(i)*z0max) )

         ! TODO: add surface perturbations to ztmax over land
         ! it is all zero for control_p8 configuration
         ztmax = max(ztmax, zmin)

         ! compute a function of surface roughness & green vegetation fraction (zvfun)
         tem1 = (z0max-z0lo)/(z0up-z0lo)
         tem1 = min(max(tem1, 0.0_r8), 1.0_r8)
         tem2 = max(noahmp%model%sigmaf(i), 0.1_r8)
         noahmp%model%zvfun(i) = sqrt(tem1*tem2)

         ! call stability function
         call stability( &
         !  ---  inputs:
              noahmp%model%zf(i), noahmp%model%zvfun(i) , sqrt(noahmp%domain%garea(i)), &
              tv1               , thv1                  , noahmp%forc%wind(i)         , &
              z0max             , ztmax                 , tvs                         , &
              con_g             , noahmp%model%thsfc_loc,                               &
         !  ---  outputs:
              noahmp%model%rb1(i)  , noahmp%model%fm1(i)      , noahmp%model%fh1(i)   , &
              noahmp%model%fm101(i), noahmp%model%fh21(i)     , noahmp%model%cm(i)    , &
              noahmp%model%ch(i)   , noahmp%model%stress1(i)  , noahmp%model%ustar1(i))
       end if
    end do

    !----------------------
    ! write out initial conditions in case of debugging
    !----------------------

    if (first_time) then
       write(filename, fmt='(a,i4,a1,i2.2,a1,i2.2,a1,i5.5)') &
          trim(noahmp%nmlist%case_name)//'.lnd.ini.', &
          year, '-', month, '-', day, '-', hour*60*60+minute*60+second
       call write_mosaic_output(filename, noahmp, now_time, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       first_time = .false.
    end if

    !----------------------
    ! run model
    !----------------------

    call noahmpdrv_run( &
    !  ---  inputs:
         noahmp%static%im      , noahmp%static%km       , noahmp%static%lsnowl   , &
         noahmp%static%itime   , &
         noahmp%forc%ps        , noahmp%model%u1        , noahmp%model%v1        , &
         noahmp%forc%t1        , noahmp%forc%q1         , noahmp%model%soiltyp   , &
         noahmp%model%vegtype  , noahmp%model%sigmaf    , noahmp%forc%dlwflx     , &
         noahmp%forc%dswsfc    , noahmp%model%snet      , noahmp%static%delt     , &
         noahmp%model%tg3      , noahmp%model%cm        , noahmp%model%ch        , &
         noahmp%model%prsl1    , noahmp%model%prslk1    , noahmp%model%prslki    , &
         noahmp%model%prsik1   , noahmp%model%zf        , noahmp%model%pblh      , &
         noahmp%model%dry      , noahmp%forc%wind       , noahmp%model%slopetyp  , &
         noahmp%model%shdmin   , noahmp%model%shdmax    , noahmp%model%snoalb    , &
         noahmp%model%sfalb    , noahmp%model%flag_iter , con_g                  , &
         noahmp%static%idveg   , noahmp%static%iopt_crs , noahmp%static%iopt_btr , &
         noahmp%static%iopt_run, noahmp%static%iopt_sfc , noahmp%static%iopt_frz , &
         noahmp%static%iopt_inf, noahmp%static%iopt_rad , noahmp%static%iopt_alb , &
         noahmp%static%iopt_snf, noahmp%static%iopt_tbot, noahmp%static%iopt_stc , &
         noahmp%static%iopt_trs, &
         noahmp%model%xlatin   , noahmp%model%xcoszin   , noahmp%model%iyrlen    , &
         noahmp%model%julian   , noahmp%domain%garea    , &
         noahmp%model%rainn_mp , noahmp%model%rainc_mp  , &
         noahmp%model%snow_mp  , noahmp%model%graupel_mp, noahmp%model%ice_mp    , &
         noahmp%model%rhonewsn1, &
         con_hvap              , con_cp                 , con_jcal               , &
         rhoh2o                , con_eps                , con_epsm1              , &
         con_fvirt             , con_rd                 , con_hfus               , &
         noahmp%model%thsfc_loc, &
    !  ---  in/outs:
         noahmp%model%weasd    , noahmp%model%snwdph    , noahmp%model%tskin     , &
         noahmp%model%tprcp    , noahmp%model%srflag    , noahmp%model%smc       , &
         noahmp%model%stc      , noahmp%model%slc       , noahmp%model%canopy    , &
         noahmp%model%trans    , noahmp%model%tsurf     , noahmp%model%zorl      , &
         noahmp%model%rb1      , noahmp%model%fm1       , noahmp%model%fh1       , &
         noahmp%model%ustar1   , noahmp%model%stress1   , noahmp%model%fm101     , &
         noahmp%model%fh21     , &
         noahmp%model%rmol1    , noahmp%model%flhc1     , noahmp%model%flqc1     , &
         noahmp%model%do_mynnsfclay, &
    ! --- Noah MP specific
         noahmp%model%snowxy   , noahmp%model%tvxy      , noahmp%model%tgxy      , &
         noahmp%model%canicexy , noahmp%model%canliqxy  , noahmp%model%eahxy     , &
         noahmp%model%tahxy    , noahmp%model%cmxy      , noahmp%model%chxy      , &
         noahmp%model%fwetxy   , noahmp%model%sneqvoxy  , noahmp%model%alboldxy  , &
         noahmp%model%qsnowxy  , noahmp%model%wslakexy  , noahmp%model%zwtxy     , &
         noahmp%model%waxy     , noahmp%model%wtxy      , noahmp%model%tsnoxy    , &
         noahmp%model%zsnsoxy  , noahmp%model%snicexy   , noahmp%model%snliqxy   , &
         noahmp%model%lfmassxy , noahmp%model%rtmassxy  , noahmp%model%stmassxy  , &
         noahmp%model%woodxy   , noahmp%model%stblcpxy  , noahmp%model%fastcpxy  , &
         noahmp%model%xlaixy   , noahmp%model%xsaixy    , noahmp%model%taussxy   , &
         noahmp%model%smoiseq  , noahmp%model%smcwtdxy  , noahmp%model%deeprechxy, &
         noahmp%model%rechxy   , noahmp%model%albdvis   , noahmp%model%albdnir   , &
         noahmp%model%albivis  , noahmp%model%albinir   , noahmp%model%emiss     , &
    !  ---  outputs:
         noahmp%model%sncovr1  , noahmp%model%qsurf     , noahmp%model%gflux     , &
         noahmp%model%drain    , noahmp%model%evap      , noahmp%model%hflx      , &
         noahmp%model%ep       , noahmp%model%runoff    , noahmp%model%cmm       , &
         noahmp%model%chh      , noahmp%model%evbs      , noahmp%model%evcw      , &
         noahmp%model%sbsno    , noahmp%model%pah       , noahmp%model%ecan      , &
         noahmp%model%etran    , noahmp%model%edir      , &
         noahmp%model%snowc    , noahmp%model%stm       , &
         noahmp%model%snohf    , noahmp%model%smcwlt2   , noahmp%model%smcref2   , &
         noahmp%model%wet1     , noahmp%model%t2mmp     , noahmp%model%q2mp      , &
         noahmp%model%zvfun    , noahmp%model%ztmax     , &
         noahmp%static%errmsg   , noahmp%static%errflg)

    !----------------------
    ! unit conversions
    !----------------------

    noahmp%model%rho = noahmp%model%prsl1/(con_rd*noahmp%forc%t1*(1.0_r8+con_fvirt*noahmp%forc%q1))
    noahmp%model%hflx = noahmp%model%hflx*noahmp%model%rho*con_cp
    noahmp%model%evap = noahmp%model%evap*noahmp%model%rho*con_hvap
    where(noahmp%forc%dswsfc>0.0_r8 .and. noahmp%model%sfalb<0.0_r8) noahmp%forc%dswsfc = 0.0_r8

    !----------------------
    ! write output 
    !---------------------- 

    ! return date to create file name
    call ESMF_TimeGet(currTime+timeStep, yy=year, mm=month, dd=day, &
      h=hour, m=minute, s=second, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! get as second
    call ESMF_TimeIntervalGet(currTime-epoc+timeStep, s_r8=now_time, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! check the output frequency before calling write method
    if (mod(int(now_time), noahmp%nmlist%output_freq) == 0) then
       write(filename, fmt='(a,i4,a1,i2.2,a1,i2.2,a1,i5.5)') &
          trim(noahmp%nmlist%case_name)//'.lnd.out.', &
          year, '-', month, '-', day, '-', hour*60*60+minute*60+second 
       call write_mosaic_output(filename, noahmp, now_time, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    !----------------------
    ! exit if there is an error 
    !----------------------

    if (noahmp%static%errflg /= 0) then
       ! print error message to stdout
       write(*,*) "noahmpdrv_run reporting an error"
       write(*,*) noahmp%static%errmsg

       ! finalize ESMF
       call ESMF_Finalize(endflag=ESMF_END_ABORT)
    end if

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine drv_run

  !===============================================================================
  subroutine interpolate_monthly(currTime, vector_length, monthly_var, interp_var, rc)

    ! input/output variables
    type(ESMF_Time), intent(in) :: currTime
    integer , intent(in)        :: vector_length
    real(r8), intent(in)        :: monthly_var(vector_length,12)
    real(r8), intent(inout)     :: interp_var(vector_length)
    integer, intent(inout)      :: rc

    ! local variables
    integer                     :: iloop, yy, pre_mm, mm, after_mm, dd
    real(r8)                    :: pre_time, now_time, after_time
    real(r8)                    :: pre_wgt, after_wgt
    character(len=*), parameter :: subname=trim(modName)//':(interpolate_monthly) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    !----------------------
    ! query current time to find out previous and next times 
    !---------------------- 

    call ESMF_TimeGet(currTime, yy=yy, mm=mm, dd=dd, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_TimeSet(prevTime, yy=yy, mm=mm, dd=15, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_TimeSet(nextTime, yy=yy, mm=mm, dd=15, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return     

    if (dd < 15) then
       if (mm == 1) then
          call ESMF_TimeSet(prevTime, yy=yy-1, mm=12, dd=15, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       else
          call ESMF_TimeSet(prevTime, yy=yy, mm=mm-1, dd=15, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if
    else
       if (mm == 12) then
          call ESMF_TimeSet(nextTime, yy=yy+1, mm=1, dd=15, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       else
          call ESMF_TimeSet(nextTime, yy=yy, mm=mm+1, dd=15, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if
    end if

    call ESMF_TimeGet(prevTime, mm=pre_mm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_TimeGet(nextTime, mm=after_mm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_TimeIntervalGet(currTime-epoc, s_r8=now_time, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_TimeIntervalGet(prevTime-epoc, s_r8=pre_time, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_TimeIntervalGet(nextTime-epoc, s_r8=after_time, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !---------------------- 
    ! create interpolation weights
    !---------------------- 

    after_wgt = (now_time-pre_time)/(after_time-pre_time)
    pre_wgt = 1.0-after_wgt
    
    !---------------------- 
    ! perform interpolation
    !---------------------- 

    do iloop = 1, vector_length
       interp_var(iloop) = pre_wgt*monthly_var(iloop,pre_mm)+after_wgt*monthly_var(iloop,after_mm)
    end do

  end subroutine interpolate_monthly

  !===============================================================================
  subroutine calc_cosine_zenith(currTime, vector_length, latitude, longitude, cosz, julian, rc)
    implicit none

    ! input/output variables
    type(ESMF_Time), intent(in) :: currTime
    integer , intent(in)        :: vector_length
    real(r8), intent(in)        :: latitude(vector_length)
    real(r8), intent(in)        :: longitude(vector_length)
    real(r8), intent(inout)     :: cosz(vector_length)
    real(r8), intent(inout)     :: julian
    integer, intent(inout)      :: rc

    ! local variables
    real(r8)                    :: sec_since
    real(r8)                    :: obecl, sinob, sxlong, arg, tloctim, hrang, declin
    integer                     :: iloc, iyear, imonth, iday, ihour, iminute, isecond
    real(r8), parameter         :: degrad = 3.14159265/180.0
    real(r8), parameter         :: dpd    = 360.0/365.0
    character(len=*), parameter :: subname=trim(modName)//':(calc_cosine_zenith) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    !----------------------
    ! query current time
    !---------------------- 

    call ESMF_TimeGet(currTime, yy=iyear, mm=imonth, dd=iday, h=ihour, m=iminute, s=isecond, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_TimeSet(dummTime, yy=iyear, mm=1, dd=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_TimeIntervalGet(currTime-dummTime, s_r8=sec_since, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------
    ! obliquity = 23.5 degree.
    !----------------------

    obecl = 23.5*degrad
    sinob = sin(obecl)

    !----------------------
    ! calculate longitude of the sun from vernal equinox
    !----------------------

    julian = sec_since/86400.0
    if (julian >= 80.0) sxlong = dpd*(julian-80.0)*degrad
    if (julian <  80.0) sxlong = dpd*(julian+285.0)*degrad
    arg = sinob*sin(sxlong)
    declin = asin(arg)

    do iloc = 1, vector_length
       tloctim = real(ihour)+real(iminute)/60.0+real(isecond)/3600.0+ longitude(iloc)/15.0 ! local time in hours
       tloctim = mod(tloctim+24.0, 24.0)
       hrang = 15.*(tloctim-12.0)*degrad
       cosz(iloc) = sin(latitude(iloc)*degrad)*sin(declin)+cos(latitude(iloc)*degrad)*cos(declin)*cos(hrang)
    end do

  end subroutine calc_cosine_zenith

  !===============================================================================
  function find_eq_smc(bexp, dwsat, dksat, ddz, smcmax) result(smc)
    !
    ! Use newton-raphson method to find eq soil moisture
    !
    implicit none

    ! input/output variables
    real(r8), intent(in) :: bexp, dwsat, dksat, ddz, smcmax

    ! local variables
    integer  :: iter
    real(r8) :: smc
    real(r8) :: expon, aa, bb, func, dfunc, dx
    character(len=*), parameter :: subname=trim(modName)//':(find_eq_smc) '
    !-------------------------------------------------------------------------------

    expon = bexp + 1.0_r8
    aa    = dwsat / ddz
    bb    = dksat / smcmax ** expon
    smc = 0.5_r8 * smcmax
    
    do iter = 1, 100
       func  = (smc - smcmax) * aa +  bb * smc ** expon
       dfunc = aa + bb * expon * smc ** bexp
       dx    = func / dfunc
       smc   = smc - dx
       if (abs(dx) < 1.e-6_r8) return
    end do

  end function find_eq_smc

end module lnd_comp_driver
