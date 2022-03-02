module lnd_comp_driver

  ! This file contains the NoahMP land surface model driver

  use ESMF            , only: operator(+), operator(-), operator(/)
  use ESMF            , only: ESMF_LogWrite, ESMF_SUCCESS, ESMF_LOGMSG_INFO
  use ESMF            , only: ESMF_Finalize, ESMF_END_ABORT
  use ESMF            , only: ESMF_VM, ESMF_VMCommWaitAll
  use ESMF            , only: ESMF_GridComp, ESMF_GridCompGet
  use ESMF            , only: ESMF_Clock, ESMF_ClockGet, ESMF_TimeIsLeapYear
  use ESMF            , only: ESMF_Time, ESMF_TimeGet, ESMF_TimeSet
  use ESMF            , only: ESMF_TimeInterval, ESMF_TimeIntervalGet
  use ESMF            , only: ESMF_ClockGetAlarm, ESMF_AlarmIsCreated
  use ESMF            , only: ESMF_AlarmIsRinging, ESMF_Alarm, ESMF_AlarmRingerOff

  use lnd_comp_kind   , only: r8 => shr_kind_r8
  use lnd_comp_kind   , only: cl => shr_kind_cl
  use lnd_comp_types  , only: noahmp_type
  use lnd_comp_shr    , only: chkerr
  use lnd_comp_io     , only: read_static, read_initial, read_restart
  use lnd_comp_io     , only: write_mosaic_output

  use noahmpdrv       , only: noahmpdrv_run
  use namelist_soilveg, only: z0_data
  use set_soilveg_mod , only: set_soilveg
  use funcphys        , only: gpvs
  use physcons        , only: con_hvap, con_cp, con_jcal
  use physcons        , only: con_eps, con_epsm1, con_fvirt
  use physcons        , only: con_rd, con_hfus
  use physcons        , only: rhoh2o => rhowater
  use physcons        , only: p0 => con_p0
  use physcons        , only: cappa => con_rocp
  use physcons        , only: con_g

  implicit none
  private

  public :: drv_init, drv_run, drv_finalize

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
    character(len=*), parameter :: subname=trim(modName)//':(drv_init) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

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
    noahmp%model%flag_iter  = .true.
    noahmp%model%thsfc_loc  = .true.

    !----------------------
    ! initialize soil and vegetation
    !----------------------

    call set_soilveg(0, noahmp%static%isot, noahmp%static%ivegsrc, 0)
    call gpvs()

    !----------------------
    ! unit conversion
    !----------------------

    ! at driver level, roughness length in cm
    if (.not. noahmp%nmlist%restart_run) then
       noahmp%model%zorl = 0.0_r8
       where(noahmp%model%vegtype > 0) noahmp%model%zorl = z0_data(noahmp%model%vegtype)*100.0_r8
    end if

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine drv_init

  !===============================================================================
  subroutine drv_run(gcomp, noahmp, rc)

    ! input/output variables
    type(ESMF_GridComp), intent(in)    :: gcomp
    type(noahmp_type)  , intent(inout) :: noahmp 
    integer            , intent(out)   :: rc

    ! local variables
    integer, save               :: first_time = .true.
    integer                     :: i, step
    integer                     :: year, month, day, hour, minute, second
    real(r8)                    :: now_time
    real(r8), save, allocatable :: rho(:)
    character(len=cl)           :: filename
    logical                     :: restart_write
    type(ESMF_VM)               :: vm
    type(ESMF_Clock)            :: clock
    type(ESMF_Alarm)            :: alarm
    type(ESMF_Time)             :: startTime
    type(ESMF_Time)             :: currTime
    type(ESMF_TimeInterval)     :: timeStep
    real(r8), parameter         :: tfreeze = 2.7315e+2_r8 
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
       ! initialize model variables 
       call noahmp%InitializeStates(noahmp%nmlist, noahmp%static, month)

       ! transfer initial conditions
       do i = noahmp%domain%begl, noahmp%domain%endl
          if (noahmp%domain%mask(i) == 1) then
             noahmp%model%weasd(i)  = noahmp%init%snow_water_equivalent(i)
             noahmp%model%snwdph(i) = noahmp%init%snow_depth(i)*1000.0_r8
             noahmp%model%canopy(i) = noahmp%init%canopy_water(i)
             noahmp%model%tskin(i)  = noahmp%init%skin_temperature(i)
             noahmp%model%stc(i,:)  = noahmp%init%soil_temperature(i,:)
             noahmp%model%smc(i,:)  = noahmp%init%soil_moisture(i,:)
             noahmp%model%slc(i,:)  = noahmp%init%soil_liquid(i,:)
          end if
       end do
    end if

    ! TODO: CDEPS data atmosphere Sa_z (noahmp%forc%hgt) is 30 meters but UFS land driver uses 10 meters?
    ! There could be option in nems.configure to overwrite Sa_z
    noahmp%model%zf = 10.0_r8

    !----------------------
    ! allocate required temporary variable
    !----------------------

    if (.not. allocated(rho)) allocate(rho(noahmp%domain%begl:noahmp%domain%endl))

    !----------------------
    ! interpolate monthly data, vegetation fraction and mean sfc diffuse sw albedo (NOT used)
    !----------------------

    call interpolate_monthly(currTime, noahmp%static%im, &
       noahmp%model%gvf_monthly, noahmp%model%sigmaf, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
 
    ! not used by the model but is set internally in the driver to albedo_total
    call interpolate_monthly(currTime, noahmp%static%im, &
       noahmp%model%alb_monthly, noahmp%model%sfalb, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------
    ! calculate solar zenith angle
    !----------------------

    call calc_cosine_zenith(currTime, noahmp%static%im, &
      noahmp%domain%lats, noahmp%domain%lons, noahmp%model%xcoszin, &
      noahmp%model%julian, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------
    ! set variables from forcing  
    !----------------------

    ! since wind direction is not important for NoahMP, just provide speed
    noahmp%model%u1 = noahmp%forc%wind
    noahmp%model%v1 = 0.0_r8

    noahmp%model%snet = noahmp%forc%dswsfc*(1.0_r8-noahmp%model%sfalb)
    noahmp%model%srflag = 0.0_r8
    where(noahmp%forc%t1 < tfreeze) noahmp%model%srflag = 1.0_r8
    noahmp%model%prsl1 = noahmp%forc%ps*exp(-1.0_r8*noahmp%model%zf/29.25_r8/noahmp%forc%t1) 
    noahmp%model%prslki = (exp(noahmp%model%zf/29.25_r8/noahmp%forc%t1))**(2.0_r8/7.0_r8)    
    noahmp%model%prslk1 = (exp(noahmp%model%zf/29.25_r8/noahmp%forc%t1))**(2.0_r8/7.0_r8)
    noahmp%model%prsik1 = (exp(noahmp%model%zf/29.25_r8/noahmp%forc%t1))**(2.0_r8/7.0_r8)
    noahmp%model%rainn_mp = 1000.0_r8*noahmp%forc%tprcp/noahmp%static%delt
    noahmp%model%rainc_mp = 0.0_r8
    noahmp%model%snow_mp = 0.0_r8
    noahmp%model%graupel_mp = 0.0_r8
    noahmp%model%ice_mp = 0.0_r8

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
         noahmp%static%im      , noahmp%static%km       , noahmp%static%lsnowl   , &
         noahmp%static%itime   , &
         noahmp%forc%ps        , noahmp%model%u1        , noahmp%model%v1        , &
         noahmp%forc%t1        , noahmp%forc%q1         , noahmp%model%soiltyp   , &
         noahmp%model%vegtype  , noahmp%model%sigmaf    , noahmp%forc%dlwflx     , &
         noahmp%forc%dswsfc    , noahmp%model%snet      , noahmp%static%delt     , &
         noahmp%model%tg3      , noahmp%model%cm        , noahmp%model%ch        , &
         noahmp%model%prsl1    , noahmp%model%prslk1    , noahmp%model%prslki    , &
         noahmp%model%prsik1   , noahmp%model%zf        , &
         noahmp%model%dry      , noahmp%forc%wind       , noahmp%model%slopetyp  , &
         noahmp%model%shdmin   , noahmp%model%shdmax    , noahmp%model%snoalb    , &
         noahmp%model%sfalb    , noahmp%model%flag_iter , con_g                  , &
         noahmp%static%idveg   , noahmp%static%iopt_crs , noahmp%static%iopt_btr , &
         noahmp%static%iopt_run, noahmp%static%iopt_sfc , noahmp%static%iopt_frz , &
         noahmp%static%iopt_inf, noahmp%static%iopt_rad , noahmp%static%iopt_alb , &
         noahmp%static%iopt_snf, noahmp%static%iopt_tbot, noahmp%static%iopt_stc , &
         noahmp%model%xlatin   , noahmp%model%xcoszin   , noahmp%model%iyrlen    , &
         noahmp%model%julian   , noahmp%domain%garea    , & 
         noahmp%model%rainn_mp , noahmp%model%rainc_mp  , &
         noahmp%model%snow_mp  , noahmp%model%graupel_mp, noahmp%model%ice_mp    , &
         con_hvap              , con_cp                 , con_jcal               , &
         rhoh2o                , con_eps                , con_epsm1              , &
         con_fvirt             , con_rd                 , con_hfus               , &
         noahmp%model%thsfc_loc, & 
         noahmp%model%weasd    , noahmp%model%snwdph    , noahmp%model%tskin     , &
         noahmp%forc%tprcp     , noahmp%model%srflag    , noahmp%model%smc       , &
         noahmp%model%stc      , noahmp%model%slc       , noahmp%model%canopy    , &
         noahmp%model%trans    , noahmp%model%tsurf     , noahmp%model%zorl      , &
         noahmp%model%rb1      , noahmp%model%fm1       , noahmp%model%fh1       , &
         noahmp%model%ustar1   , noahmp%model%stress1   , noahmp%model%fm101     , &
         noahmp%model%fh21     , &
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
         noahmp%model%sncovr1  , noahmp%model%qsurf     , noahmp%model%gflux     , &
         noahmp%model%drain    , noahmp%model%evap      , noahmp%model%hflx      , &
         noahmp%model%ep       , noahmp%model%runoff    , noahmp%model%cmm       , &
         noahmp%model%chh      , noahmp%model%evbs      , noahmp%model%evcw      , &
         noahmp%model%sbsno    , noahmp%model%pah       , noahmp%model%ecan      , &
         noahmp%model%etran    , noahmp%model%edir      , &
         noahmp%model%snowc    , noahmp%model%stm       , &
         noahmp%model%snohf    , noahmp%model%smcwlt2   , noahmp%model%smcref2   , &
         noahmp%model%wet1     , noahmp%model%t2mmp     , noahmp%model%q2mp      , &
         noahmp%model%zvfun    , noahmp%static%errmsg   , noahmp%static%errflg)

    !----------------------
    ! unit conversions
    !----------------------

    rho = noahmp%model%prsl1/(con_rd*noahmp%forc%t1*(1.0_r8+con_fvirt*noahmp%forc%q1)) 
    noahmp%model%hflx = noahmp%model%hflx*rho*con_cp
    noahmp%model%evap = noahmp%model%evap*rho*con_hvap
    where(noahmp%forc%dswsfc>0.0_r8 .and. noahmp%model%sfalb<0.0_r8) noahmp%forc%dswsfc = 0.0_r8

    !----------------------
    ! write output 
    !---------------------- 

    ! get as second
    call ESMF_TimeIntervalGet(currTime-epoc, s_r8=now_time, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! check the output frequency before calling write method
    if (mod(hour, noahmp%nmlist%output_freq) == 0) then
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
  subroutine drv_finalize()
  end subroutine drv_finalize

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

end module lnd_comp_driver
