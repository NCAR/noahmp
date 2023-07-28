module lnd_comp_nuopc

  !----------------------------------------------------------------------------
  ! This is the NUOPC cap for NoahMP
  !----------------------------------------------------------------------------

  use ESMF             , only : operator(+), operator(-), operator(==)
  use ESMF             , only : ESMF_GridCompSetEntryPoint, ESMF_GridComp
  use ESMF             , only : ESMF_GridCompGet, ESMF_VMGet, ESMF_Mesh, ESMF_MeshWriteVTK
  use ESMF             , only : ESMF_MethodRemove, ESMF_LogWrite, ESMF_LOGMSG_INFO
  use ESMF             , only : ESMF_METHOD_INITIALIZE, ESMF_SUCCESS, ESMF_FAILURE
  use ESMF             , only : ESMF_State, ESMF_Clock, ESMF_Time, ESMF_VM
  use ESMF             , only : ESMF_Array, ESMF_ArrayRead, ESMF_ArrayGet, ESMF_ArrayDestroy
  use ESMF             , only : ESMF_TimeInterval, ESMF_Alarm, ESMF_ClockGet
  use ESMF             , only : ESMF_ClockGetAlarmList, ESMF_Clock, ESMF_Time
  use ESMF             , only : ESMF_ClockSet, ESMF_TimeIntervalGet, ESMF_ALARMLIST_ALL
  use ESMF             , only : ESMF_AlarmSet, ESMF_ClockAdvance
  use ESMF             , only : ESMF_TimeGet, ESMF_TimeInterval
  use ESMF             , only : ESMF_GEOMTYPE_GRID, ESMF_GEOMTYPE_MESH
  use ESMF             , only : ESMF_MeshCreate, ESMF_Grid, ESMF_GeomType_Flag
  use ESMF             , only : ESMF_VMGet, ESMF_VMGetCurrent
  use NUOPC            , only : NUOPC_CompDerive, NUOPC_CompAttributeGet
  use NUOPC            , only : NUOPC_CompFilterPhaseMap, NUOPC_CompSetEntryPoint 
  use NUOPC            , only : NUOPC_CompSpecialize
  use NUOPC_Model      , only : NUOPC_ModelGet
  use NUOPC_Model      , only : model_routine_SS => SetServices
  use NUOPC_Model      , only : SetVM
  use NUOPC_Model      , only : model_label_Advance => label_Advance
  use NUOPC_Model      , only : model_label_SetRunClock => label_SetRunClock
  use NUOPC_Model      , only : model_label_Finalize => label_Finalize

  use lnd_comp_types   , only : noahmp_type
  use lnd_comp_kind    , only : cl => shr_kind_cl
  use lnd_comp_kind    , only : r8 => shr_kind_r8
  use lnd_comp_shr     , only : chkerr, alarm_init
  use lnd_comp_shr     , only : shr_string_listGetName, read_namelist
  use lnd_comp_domain  , only : lnd_set_decomp_and_domain_from_mosaic
  use lnd_comp_import_export, only : advertise_fields, realize_fields
  use lnd_comp_import_export, only : import_fields, export_fields, state_diagnose
  use lnd_comp_driver  , only : drv_init, drv_run, drv_finalize

  implicit none
  private ! except

  ! Module public routines
  public  :: SetServices, SetVM

  ! Module private routines
  private :: InitializeP0        ! Phase zero of initialization
  private :: InitializeAdvertise ! Advertise the fields that can be passed
  private :: InitializeRealize   ! Realize the list of fields that will be exchanged
  private :: ModelAdvance        ! Advance the model
  private :: ModelSetRunClock    ! Set the run clock

  !--------------------------------------------------------------------------
  ! Private module data
  !--------------------------------------------------------------------------

  integer           :: dbug = 2
  type(noahmp_type) :: noahmp

  character(*),parameter :: modName =  "(lnd_comp_nuopc)"

  character(len=*) , parameter :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine SetServices(gcomp, rc)
    ! Setup the pointers to the function calls for the different models phases (initialize, run, finalize)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    character(len=*),parameter  :: subname=trim(modName)//':(SetServices) '

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! the NUOPC gcomp component will register the generic methods
    call NUOPC_CompDerive(gcomp, model_routine_SS, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! switching to IPD versions
    call ESMF_GridCompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         userRoutine=InitializeP0, phase=0, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! set entry point for methods that require specific implementation
    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         phaseLabelList=(/"IPDv01p1"/), userRoutine=InitializeAdvertise, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         phaseLabelList=(/"IPDv01p3"/), userRoutine=InitializeRealize, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! attach specializing method(s)
    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Advance, &
         specRoutine=ModelAdvance, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_MethodRemove(gcomp, label=model_label_SetRunClock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_SetRunClock, &
         specRoutine=ModelSetRunClock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Finalize, &
         specRoutine=ModelFinalize, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine SetServices

  !===============================================================================
  subroutine InitializeP0(gcomp, importState, exportState, clock, rc)

    ! Phase zero initialization
    ! input/output variables
    type(ESMF_GridComp)   :: gcomp
    type(ESMF_State)      :: importState, exportState
    type(ESMF_Clock)      :: clock
    integer, intent(out)  :: rc
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! Switch to IPDv01 by filtering all other phaseMap entries
    call NUOPC_CompFilterPhaseMap(gcomp, ESMF_METHOD_INITIALIZE, acceptStringList=(/"IPDv01p"/), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine InitializeP0

  !===============================================================================
  subroutine InitializeAdvertise(gcomp, importState, exportState, clock, rc)

    ! Advertise the fields that can be exchanged
    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_VM)               :: vm
    integer                     :: lnd_mpi_comm
    character(len=CL)           :: cvalue
    character(len=*), parameter :: subname=trim(modName)//':(InitializeAdvertise) '
    character(len=*), parameter :: format = "('("//trim(subname)//") :',A)"
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! ---------------------
    ! Advertise fields
    ! ---------------------

    call advertise_fields(gcomp, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine InitializeAdvertise

  !===============================================================================
  subroutine InitializeRealize(gcomp, importState, exportState, clock, rc)

    ! Realize the list of fields that will be exchanged
    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState
    type(ESMF_State)     :: exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! local variables
    integer                    :: n, suffix_sec
    real(r8)                   :: dt
    integer, pointer           :: ptr(:)
    character(len=CL)          :: cvalue, cname, msg
    character(len=CL)          :: meshfile_mask
    character(len=CL)          :: model_meshfile
    character(len=CL)          :: model_mosaicfile
    character(len=CL)          :: input_dir
    logical                    :: isPresent, isSet
    integer                    :: year, month, day
    integer                    :: hour, minute, second
    type(ESMF_Time)            :: currTime
    type(ESMF_TimeInterval)    :: timeStep
    character(len=*),parameter :: subname=trim(modName)//':(InitializeRealize) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! ----------------------
    ! Get model physics related namelist options
    ! ----------------------

    call read_namelist(gcomp, noahmp, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! ---------------------
    ! Query input directory and name of mosaic file
    ! ---------------------

    call NUOPC_CompAttributeGet(gcomp, name='mosaic_file', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (isPresent .and. isSet) then
      noahmp%nmlist%mosaic_file = trim(cvalue)
      call ESMF_LogWrite(trim(subname)//': mosaic file = '//trim(noahmp%nmlist%mosaic_file), ESMF_LOGMSG_INFO)
    else
      call ESMF_LogWrite(trim(subname)//': mosaic_file is required! Exiting ....', ESMF_LOGMSG_INFO)
      rc = ESMF_FAILURE
      return
    end if

    call NUOPC_CompAttributeGet(gcomp, name='input_dir', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (isPresent .and. isSet) then
       noahmp%nmlist%input_dir = trim(cvalue)
    else
       noahmp%nmlist%input_dir = "INPUT/"
    end if
    call ESMF_LogWrite(trim(subname)//': input_dir = '//trim(noahmp%nmlist%input_dir), ESMF_LOGMSG_INFO)

    call NUOPC_CompAttributeGet(gcomp, name='restart_dir', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (isPresent .and. isSet) then
       noahmp%nmlist%restart_dir = trim(cvalue)
    else
       noahmp%nmlist%restart_dir = "RESTART/"
    end if
    call ESMF_LogWrite(trim(subname)//': restart_dir = '//trim(noahmp%nmlist%restart_dir), ESMF_LOGMSG_INFO)

    ! ---------------------
    ! Query ESMF attribute, layout 
    ! ---------------------

    call NUOPC_CompAttributeGet(gcomp, name='layout', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       do n = 1, 2
          call shr_string_listGetName(cvalue, n, cname, rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          read(cname,*) noahmp%domain%layout(n)
          write(msg, fmt='(A,I1,A)') trim(subname)//': layout(',n,') = '//trim(cname)
          call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO)
       end do
    else
       noahmp%domain%layout(:) = -1
    end if

    ! ---------------------
    ! Check run is restart or not?
    ! ---------------------

    call NUOPC_CompAttributeGet(gcomp, name='start_type', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (isPresent .and. isSet) then
       if (trim(cvalue) == 'startup') then
          noahmp%nmlist%restart_run = .false.
       elseif (trim(cvalue) == 'continue') then
          noahmp%nmlist%restart_run = .true.
       end if
    else
       noahmp%nmlist%restart_run = .false. 
    end if
    write(msg, fmt='(A,L)') trim(subname)//': restart_run = ', noahmp%nmlist%restart_run
    call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO)

    if (noahmp%nmlist%restart_run) then
       call NUOPC_CompAttributeGet(gcomp, name='restart_file', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       if (isPresent .and. isSet) then
          noahmp%nmlist%restart_file = trim(cvalue)//'.tile*.nc'
       else
          call ESMF_ClockGet(clock, currTime=currTime, timeStep=timeStep, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          call ESMF_TimeGet(currTime, yy=year, mm=month, dd=day, h=hour, m=minute, s=second, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          call ESMF_TimeIntervalGet(timeStep, s_r8=dt, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! update day if it is required
          ! coupling time (dt) needs to be same in case of restart run
          suffix_sec = int(hour*60*60+minute*60+second)
          if (suffix_sec < 0) then
             day = day-1
             suffix_sec = 86400-abs(suffix_sec)
          end if

          write(noahmp%nmlist%restart_file, fmt='(a,i4,a1,i2.2,a1,i2.2,a1,i5.5,a)') &
             trim(noahmp%nmlist%case_name)//'.lnd.out.', year, '-', month, '-', day, '-', suffix_sec, '.tile*.nc'
       end if

       call ESMF_LogWrite(trim(subname)//': restart_file = '//trim(noahmp%nmlist%restart_file), ESMF_LOGMSG_INFO)
    end if

    ! ---------------------
    ! Option to handle initial conditions, *.initial.* or sfc_data.*
    ! ---------------------

    call NUOPC_CompAttributeGet(gcomp, name='ic_type', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (isPresent .and. isSet) then
       noahmp%nmlist%ic_type = trim(cvalue)
       if (trim(noahmp%nmlist%ic_type) /= 'sfc' .and. trim(noahmp%nmlist%ic_type) /= 'custom') then
          call ESMF_LogWrite(trim(subname)//': '//trim(noahmp%nmlist%ic_type)//' is not a valid! It must be [sfc|custom].', ESMF_LOGMSG_INFO)
          rc = ESMF_FAILURE
          return
       end if
    else
       noahmp%nmlist%ic_type = 'sfc'
    end if

    call ESMF_LogWrite(trim(subname)//': ic_type = '//trim(noahmp%nmlist%ic_type), ESMF_LOGMSG_INFO)

    ! ---------------------
    ! Option to disable export fields
    ! ---------------------

    call NUOPC_CompAttributeGet(gcomp, name='has_export', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    noahmp%nmlist%has_export = .true.
    if (isPresent .and. isSet) then
       if (trim(cvalue) .eq. '.false.' .or. trim(cvalue) .eq. 'false') noahmp%nmlist%has_export = .false.
    end if

    write(msg, fmt='(A,L)') trim(subname)//': has_export = ', noahmp%nmlist%has_export
    call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO)

    ! ---------------------
    ! Option to calculate net shortwave internally using downwelling component and albedo
    ! ---------------------

    call NUOPC_CompAttributeGet(gcomp, name='calc_snet', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    noahmp%nmlist%calc_snet = .false.
    if (isPresent .and. isSet) then
       if (trim(cvalue) .eq. '.true.' .or. trim(cvalue) .eq. 'true') noahmp%nmlist%calc_snet = .true.
    end if

    write(msg, fmt='(A,L)') trim(subname)//': calc_snet = ', noahmp%nmlist%calc_snet
    call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO)

    ! ---------------------
    ! Create mosaic grid and convert it to mesh 
    ! ---------------------

    call lnd_set_decomp_and_domain_from_mosaic(gcomp, noahmp, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! check mesh for debugging purposes
    if (dbug > 2) then
       call ESMF_MeshWriteVTK(noahmp%domain%mesh, "lnd_mesh", rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    !----------------------
    ! Initialize NoahMP
    !----------------------

    call drv_init(gcomp, noahmp, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! ---------------------
    ! Realize the actively coupled fields
    ! ---------------------

    call realize_fields(importState, exportState, noahmp%domain%mesh, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! ---------------------
    ! Create land export state
    ! ---------------------

    if (noahmp%nmlist%has_export) then
       call export_fields(gcomp, noahmp, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    ! ---------------------
    ! Diagnostics
    ! ---------------------

    if (dbug > 0) then
       call state_diagnose(exportState, subname//': ExportState ',rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    endif

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine InitializeRealize

  !===============================================================================
  subroutine ModelAdvance(gcomp, rc)

    !------------------------
    ! Run NoahMP
    !------------------------

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Clock)           :: clock
    type(ESMF_Time)            :: currTime
    type(ESMF_Time)            :: nextTime
    type(ESMF_State)           :: importState, exportState
    character(len=*),parameter :: subname=trim(modName)//':(ModelAdvance) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    !-----------------------
    ! Query the Component for its clock, importState and exportState
    !-----------------------

    call NUOPC_ModelGet(gcomp, modelClock=clock, importState=importState, exportState=exportState, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !-----------------------
    ! import state
    !-----------------------

    call import_fields(gcomp, noahmp, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------
    ! diagnostics
    !----------------------

    if (dbug > 1) then
       call state_diagnose(importState, subname//': ImportState ',rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    endif

    !----------------------
    ! run NoahMP
    !----------------------
    call drv_run(gcomp, noahmp, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------
    ! export state
    !----------------------

    if (noahmp%nmlist%has_export) then
       call export_fields(gcomp, noahmp, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    !----------------------
    ! diagnostics
    !----------------------

    if (dbug > 1) then
       call state_diagnose(exportState, subname//': ExportState ',rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    endif

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine ModelAdvance

  !===============================================================================
  subroutine ModelSetRunClock(gcomp, rc)

    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Clock)         :: mclock, dclock
    type(ESMF_Time)          :: mcurrtime, dcurrtime
    type(ESMF_Time)          :: mstoptime
    type(ESMF_TimeInterval)  :: mtimestep, dtimestep
    character(len=256)       :: cvalue
    character(len=256)       :: restart_option ! Restart option units
    integer                  :: restart_n      ! Number until restart interval
    integer                  :: restart_ymd    ! Restart date (YYYYMMDD)
    type(ESMF_Alarm)         :: restart_alarm
    character(len=256)       :: stop_option    ! Stop option units
    integer                  :: stop_n         ! Number until stop interval
    integer                  :: stop_ymd       ! Stop date (YYYYMMDD)
    type(ESMF_Alarm)         :: stop_alarm
    character(len=128)       :: name
    integer                  :: alarmcount
    character(len=*),parameter :: subname=trim(modName)//':(ModelSetRunClock) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! query the Component for its clocks
    call NUOPC_ModelGet(gcomp, driverClock=dclock, modelClock=mclock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGet(dclock, currTime=dcurrtime, timeStep=dtimestep, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGet(mclock, currTime=mcurrtime, timeStep=mtimestep, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! force model clock currtime and timestep to match driver and set stoptime
    !--------------------------------

    mstoptime = mcurrtime + dtimestep
    call ESMF_ClockSet(mclock, currTime=dcurrtime, timeStep=dtimestep, stopTime=mstoptime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! set restart and stop alarms
    !--------------------------------

    call ESMF_ClockGetAlarmList(mclock, alarmlistflag=ESMF_ALARMLIST_ALL, alarmCount=alarmCount, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (alarmCount == 0) then
       call ESMF_GridCompGet(gcomp, name=name, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_LogWrite(subname//'setting alarms for ' // trim(name), ESMF_LOGMSG_INFO)

       !----------------
       ! Restart alarm
       !----------------
       call NUOPC_CompAttributeGet(gcomp, name="restart_option", value=restart_option, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call NUOPC_CompAttributeGet(gcomp, name="restart_n", value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) restart_n

       call NUOPC_CompAttributeGet(gcomp, name="restart_ymd", value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) restart_ymd

       call alarm_init(mclock, restart_alarm, restart_option, &
            opt_n   = restart_n,           &
            opt_ymd = restart_ymd,         &
            RefTime = mcurrTime,           &
            alarmname = 'alarm_restart', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_AlarmSet(restart_alarm, clock=mclock, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       !----------------
       ! Stop alarm
       !----------------
       call NUOPC_CompAttributeGet(gcomp, name="stop_option", value=stop_option, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call NUOPC_CompAttributeGet(gcomp, name="stop_n", value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) stop_n

       call NUOPC_CompAttributeGet(gcomp, name="stop_ymd", value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) stop_ymd

       call alarm_init(mclock, stop_alarm, stop_option, &
            opt_n   = stop_n,           &
            opt_ymd = stop_ymd,         &
            RefTime = mcurrTime,           &
            alarmname = 'alarm_stop', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_AlarmSet(stop_alarm, clock=mclock, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    !--------------------------------
    ! Advance model clock to trigger alarms then reset model clock back to currtime
    !--------------------------------

    call ESMF_ClockAdvance(mclock,rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockSet(mclock, currTime=dcurrtime, timeStep=dtimestep, stopTime=mstoptime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine ModelSetRunClock

  !===============================================================================
  subroutine ModelFinalize(gcomp, rc)

    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    character(len=*),parameter :: subname=trim(modName)//':(ModelFinalize) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! call model finalize routine
    call drv_finalize(gcomp, noahmp, rc)

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine ModelFinalize

end module lnd_comp_nuopc
