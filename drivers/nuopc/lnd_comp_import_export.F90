module lnd_comp_import_export

  ! NoahMP import and export fields exchanged with the coupler
  use ESMF          , only : ESMF_GridComp, ESMF_State, ESMF_LOGMSG_INFO
  use ESMF          , only : ESMF_LogWrite, ESMF_Finalize, ESMF_Mesh
  use ESMF          , only : ESMF_SUCCESS, ESMF_END_ABORT, ESMF_LOGMSG_ERROR
  use ESMF          , only : ESMF_Field, ESMF_GeomType_Flag, ESMF_FieldStatus_Flag
  use ESMF          , only : ESMF_StateGet, ESMF_FAILURE, ESMF_MAXSTR, ESMF_FieldGet
  use ESMF          , only : ESMF_FIELDSTATUS_COMPLETE, ESMF_GEOMTYPE_GRID
  use ESMF          , only : ESMF_FieldWriteVTK, ESMF_MeshGet, ESMF_GEOMTYPE_MESH
  use ESMF          , only : operator(==), operator(/=)
  use ESMF          , only : ESMF_StateItem_Flag, ESMF_STATEITEM_FIELD
  use NUOPC         , only : NUOPC_Advertise, NUOPC_IsConnected
  use NUOPC_Model   , only : NUOPC_ModelGet
  use lnd_comp_shr  , only : ChkErr
  use lnd_comp_kind , only : r8 => shr_kind_r8
  use lnd_comp_types, only : noahmp_type
  use lnd_comp_types, only : fld_list_type, fldsMax
  use lnd_comp_types, only : fldsToLnd, fldsToLnd_num
  use lnd_comp_types, only : fldsFrLnd, fldsFrLnd_num

  implicit none
  private ! except

  public :: advertise_fields
  public :: realize_fields
  public :: export_fields
  public :: import_fields
  public :: state_diagnose
  public :: check_for_connected

  !--------------------------------------------------------------------------
  ! Private module data
  !--------------------------------------------------------------------------

  integer                :: dbug = 2
  character(*),parameter :: modName =  "(lnd_comp_import_export)"
  character(*),parameter :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine advertise_fields(gcomp, rc)

    ! input/output variables
    type(ESMF_GridComp), intent(in)  :: gcomp
    integer,             intent(out) :: rc

    ! local variables
    type(ESMF_State)  :: importState
    type(ESMF_State)  :: exportState
    integer           :: n
    character(len=*), parameter :: subname=trim(modName)//':(advertise_fields)'
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    call NUOPC_ModelGet(gcomp, importState=importState, exportState=exportState, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! Advertise export fields
    !--------------------------------

    ! export to med
    call fldlist_add(fldsFrLnd_num, fldsFrlnd, 'Sl_lfrin')

    ! export to atm
    call fldlist_add(fldsFrLnd_num, fldsFrlnd, 'Sl_sfrac')
    call fldlist_add(fldsFrLnd_num, fldsFrlnd, 'Fall_lat')
    call fldlist_add(fldsFrLnd_num, fldsFrlnd, 'Fall_sen')
    call fldlist_add(fldsFrLnd_num, fldsFrlnd, 'Fall_evap')
    call fldlist_add(fldsFrLnd_num, fldsFrlnd, 'Sl_tref')
    call fldlist_add(fldsFrLnd_num, fldsFrlnd, 'Sl_qref')
    call fldlist_add(fldsFrLnd_num, fldsFrlnd, 'Sl_q')
    call fldlist_add(fldsFrLnd_num, fldsFrlnd, 'Fall_gflx')
    call fldlist_add(fldsFrLnd_num, fldsFrlnd, 'Fall_roff')
    call fldlist_add(fldsFrLnd_num, fldsFrlnd, 'Fall_soff')
    call fldlist_add(fldsFrLnd_num, fldsFrlnd, 'Sl_cmm')
    call fldlist_add(fldsFrLnd_num, fldsFrlnd, 'Sl_chh')
    call fldlist_add(fldsFrLnd_num, fldsFrlnd, 'Sl_zvfun')

    ! Now advertise above export fields
    do n = 1,fldsFrLnd_num
       call NUOPC_Advertise(exportState, standardName=fldsFrLnd(n)%stdname, &
            TransferOfferGeomObject='will provide', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    enddo

    !--------------------------------
    ! Advertise import fields
    !--------------------------------

    ! import from atm
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Sa_z')
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Sa_tbot')
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Sa_ta')
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Sa_tskn')
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Sa_pslv')
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Sa_prsl')
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Sa_pbot')
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Sa_shum')
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Sa_qa')
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Sa_u')
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Sa_v')
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Sa_ua')
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Sa_va')
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Sa_exner')
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Sa_ustar')
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Faxa_swdn')
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Faxa_lwdn')
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Faxa_swnet')
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Faxa_rainc')
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Faxa_rainl')
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Faxa_rain')
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Faxa_snow')
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Faxa_snowc')
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Faxa_snowl')
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'vfrac')
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'zorl')

    ! Now advertise import fields
    do n = 1,fldsToLnd_num
       call NUOPC_Advertise(importState, standardName=fldsToLnd(n)%stdname, &
            TransferOfferGeomObject='will provide', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    enddo

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine advertise_fields

  !===============================================================================
  subroutine fldlist_add(num, fldlist, stdname, ungridded_lbound, ungridded_ubound)

    ! input/output variables
    integer,                    intent(inout) :: num
    type(fld_list_type),        intent(inout) :: fldlist(:)
    character(len=*),           intent(in)    :: stdname
    integer,          optional, intent(in)    :: ungridded_lbound
    integer,          optional, intent(in)    :: ungridded_ubound

    ! local variables
    character(len=*), parameter :: subname=trim(modName)//':(fldlist_add)'
    !-------------------------------------------------------------------------------

    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! Set up a list of field information
    num = num + 1
    if (num > fldsMax) then
       call ESMF_LogWrite(trim(subname)//": ERROR num > fldsMax "//trim(stdname), &
            ESMF_LOGMSG_ERROR, line=__LINE__, file=__FILE__)
       call ESMF_Finalize(endflag=ESMF_END_ABORT)
    endif
    fldlist(num)%stdname = trim(stdname)

    if (present(ungridded_lbound) .and. present(ungridded_ubound)) then
       fldlist(num)%ungridded_lbound = ungridded_lbound
       fldlist(num)%ungridded_ubound = ungridded_ubound
    end if

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine fldlist_add

  !===============================================================================
  subroutine realize_fields(importState, exportState, Emesh, rc)

    ! input/output variables
    type(ESMF_State) , intent(inout) :: importState
    type(ESMF_State) , intent(inout) :: exportState
    type(ESMF_Mesh)  , intent(in)    :: Emesh
    integer          , intent(out)   :: rc

    ! local variables
    character(len=*), parameter :: subname=trim(modName)//':(realize_fields)'
    !---------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    call fldlist_realize( &
         state=ExportState, &
         fldList=fldsFrLnd, &
         numflds=fldsFrLnd_num, &
         tag=subname//':NoahMP_Export',&
         mesh=Emesh, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call fldlist_realize( &
         state=importState, &
         fldList=fldsToLnd, &
         numflds=fldsToLnd_num, &
         tag=subname//':NoahMP_Import',&
         mesh=Emesh, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine realize_fields

  !===============================================================================
  subroutine fldlist_realize(state, fldList, numflds, mesh, tag, rc)

    use NUOPC , only : NUOPC_IsConnected, NUOPC_Realize
    use ESMF  , only : ESMF_MeshLoc_Element, ESMF_FieldCreate, ESMF_TYPEKIND_R8
    use ESMF  , only : ESMF_MAXSTR, ESMF_Field, ESMF_State, ESMF_Mesh, ESMF_StateRemove
    use ESMF  , only : ESMF_LogFoundError, ESMF_LOGMSG_INFO, ESMF_SUCCESS
    use ESMF  , only : ESMF_LogWrite, ESMF_LOGMSG_ERROR, ESMF_LOGERR_PASSTHRU

    ! input/output variables
    type(ESMF_State)    , intent(inout) :: state
    type(fld_list_type) , intent(inout) :: fldList(:)
    integer             , intent(in)    :: numflds
    character(len=*)    , intent(in)    :: tag
    type(ESMF_Mesh)     , intent(in)    :: mesh
    integer             , intent(inout) :: rc

    ! local variables
    integer                :: n
    type(ESMF_Field)       :: field
    character(len=80)      :: stdname
    character(len=*),parameter  :: subname=trim(modName)//':fldlist_realize)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    do n = 1, numflds
       stdname = trim(fldList(n)%stdname)
       if (NUOPC_IsConnected(state, fieldName=stdname)) then
          ! Create the field
          if (fldlist(n)%ungridded_lbound > 0 .and. fldlist(n)%ungridded_ubound > 0) then
             field = ESMF_FieldCreate(mesh, ESMF_TYPEKIND_R8, name=stdname, meshloc=ESMF_MESHLOC_ELEMENT, &
                  ungriddedLbound=(/fldlist(n)%ungridded_lbound/), &
                  ungriddedUbound=(/fldlist(n)%ungridded_ubound/), &
                  gridToFieldMap=(/2/), rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          else
             field = ESMF_FieldCreate(mesh, ESMF_TYPEKIND_R8, name=stdname, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          end if
          call ESMF_LogWrite(trim(subname)//trim(tag)//" Field = "//trim(stdname)//" is connected using mesh", &
             ESMF_LOGMSG_INFO)

          ! NOW call NUOPC_Realize
          call NUOPC_Realize(state, field=field, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! Set flag for connected fields
          fldList(n)%connected = .true.
       else
          call ESMF_LogWrite(subname // trim(tag) // " Field = "// trim(stdname) // " is not connected.", &
             ESMF_LOGMSG_INFO)
          call ESMF_StateRemove(state, (/stdname/), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if
    end do

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine fldlist_realize

  !===============================================================================
  subroutine import_fields(gcomp, noahmp, rc)

    ! input/output variables
    type(ESMF_GridComp), intent(in)    :: gcomp
    type(noahmp_type),   intent(inout) :: noahmp 
    integer,             intent(out)   :: rc

    ! local variables
    type(ESMF_State)            :: importState
    character(len=*), parameter :: subname=trim(modName)//':(import_fields)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! Get export state
    call NUOPC_ModelGet(gcomp, importState=importState, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! -----------------------
    ! atm input fields
    ! -----------------------

    call state_getimport_1d(importState, 'Sa_z'      , noahmp%forc%hgt, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport_1d(importState, 'Sa_ta'     , noahmp%forc%t1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport_1d(importState, 'Sa_tbot'   , noahmp%forc%t1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport_1d(importState, 'Sa_tskn'   , noahmp%forc%tskin, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport_1d(importState, 'Sa_prsl'   , noahmp%forc%pbot, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport_1d(importState, 'Sa_pbot'   , noahmp%forc%pbot, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport_1d(importState, 'Sa_pslv'   , noahmp%forc%ps, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport_1d(importState, 'Sa_shum'   , noahmp%forc%q1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport_1d(importState, 'Sa_qa'     , noahmp%forc%q1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport_1d(importState, 'Faxa_swdn' , noahmp%forc%dswsfc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport_1d(importState, 'Faxa_lwdn' , noahmp%forc%dlwflx, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport_1d(importState, 'Faxa_swnet', noahmp%forc%snet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport_1d(importState, 'Sa_u'      , noahmp%forc%u1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport_1d(importState, 'Sa_ua'     , noahmp%forc%u1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport_1d(importState, 'Sa_v'      , noahmp%forc%v1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport_1d(importState, 'Sa_va'     , noahmp%forc%v1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport_1d(importState, 'Sa_exner'  , noahmp%forc%prslk1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport_1d(importState, 'Sa_ustar'  , noahmp%forc%ustar1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport_1d(importState, 'Faxa_rain' , noahmp%forc%tprcp, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport_1d(importState, 'Faxa_rainc', noahmp%forc%tprcpc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport_1d(importState, 'Faxa_rainl', noahmp%forc%tprcpl, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport_1d(importState, 'Faxa_snow' , noahmp%forc%snow, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport_1d(importState, 'Faxa_snowc', noahmp%forc%snowc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport_1d(importState, 'Faxa_snowl', noahmp%forc%snowl, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport_1d(importState, 'vfrac'     , noahmp%forc%vegfrac, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport_1d(importState, 'zorl'      , noahmp%forc%zorl, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine import_fields

  !===============================================================================
  subroutine export_fields(gcomp, noahmp, rc)

    ! input/output variables
    type(ESMF_GridComp), intent(in) :: gcomp
    type(noahmp_type), intent(in)   :: noahmp
    integer, intent(out)            :: rc

    ! local variables
    type(ESMF_State) :: exportState
    character(len=*), parameter :: subname=trim(modName)//':(export_fields)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! Get export state
    call NUOPC_ModelGet(gcomp, exportState=exportState, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! -----------------------
    ! output to mediator
    ! -----------------------

    call state_setexport_1d(exportState, 'Sl_lfrin', noahmp%domain%frac, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! -----------------------
    ! output to atm 
    ! -----------------------

    call state_setexport_1d(exportState, 'Sl_sfrac', noahmp%model%sncovr1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(exportState, 'Fall_lat', noahmp%model%evap, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(exportState, 'Fall_sen', noahmp%model%hflx, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(exportState, 'Fall_evap', noahmp%model%ep, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(exportState, 'Sl_tref', noahmp%model%t2mmp, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(exportState, 'Sl_qref', noahmp%model%q2mp, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(exportState, 'Sl_q', noahmp%model%qsurf, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(exportState, 'Fall_gflx', noahmp%model%gflux, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(exportState, 'Fall_roff', noahmp%model%runoff, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(exportState, 'Fall_soff', noahmp%model%drain, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(exportState, 'Sl_cmm', noahmp%model%cmm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(exportState, 'Sl_chh', noahmp%model%chh, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(exportState, 'Sl_zvfun', noahmp%model%zvfun, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine export_fields

  !===============================================================================
  subroutine state_getimport_1d(state, fldname, arr1d, rc)

    ! fill in noahmp import data for 1d field

    use ESMF, only : ESMF_LOGERR_PASSTHRU, ESMF_END_ABORT, ESMF_LogFoundError
    use ESMF, only : ESMF_Finalize

    ! input/output variabes
    type(ESMF_State) , intent(in)    :: state
    character(len=*) , intent(in)    :: fldname
    real(r8)         , intent(inout) :: arr1d(:)
    integer          , intent(out)   :: rc

    ! local variables
    real(r8), pointer           :: fldPtr1d(:)
    type(ESMF_StateItem_Flag)   :: itemType
    character(len=*), parameter :: subname=trim(modName)//':(state_getimport_1d)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    call ESMF_StateGet(state, itemName=trim(fldname), itemType=itemType, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (itemType == ESMF_STATEITEM_FIELD) then
       call state_getfldptr(State, trim(fldname), fldptr1d=fldptr1d, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       arr1d(:) = fldptr1d(:)
       call check_for_nans(arr1d, trim(fldname), 1, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       call ESMF_LogWrite(subname//' '//trim(fldname)//' is not in the state!', ESMF_LOGMSG_INFO)
    end if

  end subroutine state_getimport_1d

  !===============================================================================
  subroutine state_setexport_1d(state, fldname, arr1d, minus, rc)

    ! fill in noahmp export data for 1d field

    use ESMF, only : ESMF_LOGERR_PASSTHRU, ESMF_END_ABORT, ESMF_LogFoundError
    use ESMF, only : ESMF_Finalize

    ! input/output variabes
    type(ESMF_State) , intent(in) :: state
    character(len=*) , intent(in) :: fldname
    real(r8)         , intent(in) :: arr1d(:)
    logical, optional, intent(in) :: minus
    integer          , intent(out):: rc

    ! local variables
    logical :: l_minus ! local version of minus
    real(r8), pointer :: fldPtr1d(:)
    integer           :: g
    character(len=*), parameter :: subname='(lnd_export_export:state_setexport_1d)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    if (present(minus)) then
       l_minus = minus
    else
       l_minus = .false.
    end if

    if (.not. NUOPC_IsConnected(state, fieldName=trim(fldname))) return

    call state_getfldptr(state, trim(fldname), fldptr1d=fldptr1d, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    fldptr1d(:) = 0._r8
    if (l_minus) then
       fldptr1d(:) = -arr1d(:)
    else
       fldptr1d(:) = arr1d(:)
    end if
    call check_for_nans(arr1d, trim(fldname), 1, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine state_setexport_1d

  !===============================================================================
  subroutine state_getfldptr(state, fldname, fldptr1d, fldptr2d, rc)

    ! ----------------------------------------------
    ! Get pointer to a state field
    ! ----------------------------------------------

    use ESMF, only : ESMF_State, ESMF_Field, ESMF_Mesh, ESMF_FieldStatus_Flag
    use ESMF, only : ESMF_StateGet, ESMF_FieldGet, ESMF_MeshGet
    use ESMF, only : ESMF_FIELDSTATUS_COMPLETE, ESMF_FAILURE

    ! input/output variables
    type(ESMF_State),            intent(in)    :: State
    character(len=*),            intent(in)    :: fldname
    real(R8), pointer, optional, intent(out)   :: fldptr1d(:)
    real(R8), pointer, optional, intent(out)   :: fldptr2d(:,:)
    integer,                     intent(out)   :: rc

    ! local variables
    type(ESMF_FieldStatus_Flag) :: status
    type(ESMF_Field)            :: lfield
    character(len=*), parameter :: subname=trim(modName)//':(state_getfldptr)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    call ESMF_StateGet(State, itemName=trim(fldname), field=lfield, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (present(fldptr1d)) then
       call ESMF_FieldGet(lfield, farrayPtr=fldptr1d, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else if (present(fldptr2d)) then
       call ESMF_FieldGet(lfield, farrayPtr=fldptr2d, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       call ESMF_LogWrite(trim(subname)//": either fldptr1d or fldptr2d must be an input argument", ESMF_LOGMSG_ERROR)
       rc = ESMF_FAILURE
       return
    end if

  end subroutine state_getfldptr

  !===============================================================================
  subroutine state_diagnose(state, string, rc)

    ! ----------------------------------------------
    ! Diagnose status of State
    ! ----------------------------------------------

    type(ESMF_State), intent(in)  :: state
    character(len=*), intent(in)  :: string
    integer         , intent(out) :: rc

    ! local variables
    integer                         :: i,j,n
    type(ESMF_Field)                :: lfield
    integer                         :: fieldCount, lrank
    character(ESMF_MAXSTR) ,pointer :: lfieldnamelist(:)
    real(r8), pointer               :: dataPtr1d(:)
    real(r8), pointer               :: dataPtr2d(:,:)
    character(len=1024)             :: msgString
    character(len=*),parameter      :: subname='(state_diagnose)'
    ! ----------------------------------------------

    call ESMF_StateGet(state, itemCount=fieldCount, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    allocate(lfieldnamelist(fieldCount))

    call ESMF_StateGet(state, itemNameList=lfieldnamelist, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    do n = 1, fieldCount

       call ESMF_StateGet(state, itemName=lfieldnamelist(n), field=lfield, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       call field_getfldptr(lfield, rc=rc, fldptr1=dataPtr1d, fldptr2=dataPtr2d, rank=lrank)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       if (lrank == 0) then
          ! no local data
       elseif (lrank == 1) then
          if (size(dataPtr1d) > 0) then
             write(msgString,'(A,3g14.7,i8)') trim(string)//': '//trim(lfieldnamelist(n)), &
                  minval(dataPtr1d), maxval(dataPtr1d), sum(dataPtr1d), size(dataPtr1d)
          else
             write(msgString,'(A,a)') trim(string)//': '//trim(lfieldnamelist(n))," no data"
          endif
       elseif (lrank == 2) then
          if (size(dataPtr2d) > 0) then
             write(msgString,'(A,3g14.7,i8)') trim(string)//': '//trim(lfieldnamelist(n)), &
                  minval(dataPtr2d), maxval(dataPtr2d), sum(dataPtr2d), size(dataPtr2d)
          else
             write(msgString,'(A,a)') trim(string)//': '//trim(lfieldnamelist(n))," no data"
          endif
       else
          call ESMF_LogWrite(trim(subname)//": ERROR rank not supported ", ESMF_LOGMSG_ERROR)
          rc = ESMF_FAILURE
          return
       endif
       call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO)
    enddo

    deallocate(lfieldnamelist)

  end subroutine state_diagnose

  !===============================================================================
  subroutine field_getfldptr(field, rc, fldptr1, fldptr2, rank, abort)

    ! ----------------------------------------------
    ! for a field, determine rank and return fldptr1 or fldptr2
    ! abort is true by default and will abort if fldptr is not yet allocated in field
    ! rank returns 0, 1, or 2.  0 means fldptr not allocated and abort=false
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_Field)  , intent(in)              :: field
    integer           , intent(out)             :: rc
    real(r8), pointer , intent(inout), optional :: fldptr1(:)
    real(r8), pointer , intent(inout), optional :: fldptr2(:,:)
    integer           , intent(out)  , optional :: rank
    logical           , intent(in)   , optional :: abort

    ! local variables
    type(ESMF_GeomType_Flag)    :: geomtype
    type(ESMF_FieldStatus_Flag) :: status
    type(ESMF_Mesh)             :: lmesh
    integer                     :: lrank, nnodes, nelements
    logical                     :: labort
    character(len=*), parameter :: subname='(field_getfldptr)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    labort = .true.
    if (present(abort)) then
       labort = abort
    endif
    lrank = -99

    call ESMF_FieldGet(field, status=status, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if (status /= ESMF_FIELDSTATUS_COMPLETE) then
       lrank = 0
       if (labort) then
          call ESMF_LogWrite(trim(subname)//": ERROR data not allocated ", ESMF_LOGMSG_INFO, rc=rc)
          rc = ESMF_FAILURE
          return
       else
          call ESMF_LogWrite(trim(subname)//": WARNING data not allocated ", ESMF_LOGMSG_INFO, rc=rc)
       endif
    else

       call ESMF_FieldGet(field, geomtype=geomtype, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       if (geomtype == ESMF_GEOMTYPE_GRID) then
          call ESMF_FieldGet(field, rank=lrank, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
       elseif (geomtype == ESMF_GEOMTYPE_MESH) then
          call ESMF_FieldGet(field, rank=lrank, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call ESMF_FieldGet(field, mesh=lmesh, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call ESMF_MeshGet(lmesh, numOwnedNodes=nnodes, numOwnedElements=nelements, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          if (nnodes == 0 .and. nelements == 0) lrank = 0
       else  
          call ESMF_LogWrite(trim(subname)//": ERROR geomtype not supported ", &
               ESMF_LOGMSG_INFO, rc=rc)
          rc = ESMF_FAILURE
          return
       endif ! geomtype

       if (lrank == 0) then
          call ESMF_LogWrite(trim(subname)//": no local nodes or elements ", &
               ESMF_LOGMSG_INFO)
       elseif (lrank == 1) then
          if (.not.present(fldptr1)) then
             call ESMF_LogWrite(trim(subname)//": ERROR missing rank=1 array ", &
                  ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u)
             rc = ESMF_FAILURE
             return
          endif
          call ESMF_FieldGet(field, farrayPtr=fldptr1, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
       elseif (lrank == 2) then
          if (.not.present(fldptr2)) then
             call ESMF_LogWrite(trim(subname)//": ERROR missing rank=2 array ", &
                  ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u)
             rc = ESMF_FAILURE
             return
          endif
          call ESMF_FieldGet(field, farrayPtr=fldptr2, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
       else
          call ESMF_LogWrite(trim(subname)//": ERROR in rank ", &
               ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u)
          rc = ESMF_FAILURE
          return
       endif

    endif  ! status

    if (present(rank)) then
       rank = lrank
    endif

  end subroutine field_getfldptr

  !=============================================================================
  subroutine check_for_nans(array, fname, begg, rc)

    ! input/output variables
    real(r8)        , intent(in)  :: array(:)
    character(len=*), intent(in)  :: fname
    integer         , intent(in)  :: begg
    integer         , intent(out) :: rc

    ! local variables
    integer :: i
    character(len=*), parameter :: subname='(check_for_nans)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! Check if any input from mediator or output to mediator is NaN
    if (any(isnan(array))) then
       write(*,*) '# of NaNs = ', count(isnan(array))
       write(*,*) 'Which are NaNs = ', isnan(array)
       do i = 1, size(array)
          if (isnan(array(i))) then
             write(*,*) "NaN found in field ", trim(fname), ' at gridcell index ',begg+i-1
          end if
       end do
       call ESMF_LogWrite(trim(subname)//": one or more of the output from NoahMP to the coupler are NaN", ESMF_LOGMSG_ERROR)
       rc = ESMF_FAILURE
       return
    end if

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine check_for_nans

  !=============================================================================

  logical function check_for_connected(fldList, numflds, fname)

    ! input/output variables
    type(fld_list_type) , intent(inout) :: fldList(:)
    integer             , intent(in)    :: numflds
    character(len=*)    , intent(in)    :: fname

    ! local variables
    integer :: n
    character(len=*), parameter :: subname='(check_for_connected)'
    ! ----------------------------------------------

    check_for_connected = .false.
    do n = 1, numflds
       if (trim(fname) == trim(fldList(n)%stdname)) then
          check_for_connected = fldList(n)%connected 
          exit
       end if
    end do

  end function check_for_connected

end module lnd_comp_import_export
