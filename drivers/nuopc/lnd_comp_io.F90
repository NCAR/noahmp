#define NC_ERR_STOP(status) \
    if (status /= nf90_noerr) write(0,*) "line ", __LINE__, trim(nf90_strerror(status)); \
    if (status /= nf90_noerr) call ESMF_Finalize(endflag=ESMF_END_ABORT)

module lnd_comp_io

  ! This file contains I/O routines for the NoahMP land surface model

  use ESMF             , only : ESMF_VM, ESMF_VMGet, ESMF_VMGetCurrent
  use ESMF             , only : ESMF_TYPEKIND_R4, ESMF_KIND_R4, ESMF_INDEX_GLOBAL
  use ESMF             , only : ESMF_Grid, ESMF_Mesh, ESMF_MESHLOC_ELEMENT
  use ESMF             , only : ESMF_Field, ESMF_FieldGet, ESMF_FieldCreate, ESMF_FieldDestroy
  use ESMF             , only : ESMF_RouteHandle, ESMF_FieldRegridStore, ESMF_FieldRedist
  use ESMF             , only : ESMF_FieldWriteVTK, ESMF_STAGGERLOC_CENTER
  use ESMF             , only : ESMF_ArraySpec, ESMF_ArraySpecSet, ESMF_END_ABORT
  use ESMF             , only : ESMF_SUCCESS, ESMF_LOGMSG_INFO, ESMF_LogWrite
  use ESMF             , only : ESMF_TYPEKIND_R8, ESMF_KIND_R8, ESMF_Finalize

  use mpp_mod          , only : mpp_pe, mpp_error, mpp_sync, FATAL
  use mpp_io_mod       , only : mpp_get_info, mpp_get_fields, mpp_get_atts
  use mpp_io_mod       , only : mpp_def_dim, mpp_write, mpp_write_meta, axistype
  use mpp_io_mod       , only : fieldtype, mpp_open, mpp_read, mpp_close
  use mpp_io_mod       , only : MPP_RDONLY, MPP_NETCDF, MPP_SINGLE, MPP_MULTI, MPP_OVERWR
  use mpp_domains_mod  , only : mpp_get_compute_domain, mpp_get_domain_components, domain1D
  use mpp_parameter_mod, only : MPP_FILL_DOUBLE

  use lnd_comp_types   , only : noahmp_type
  use lnd_comp_kind    , only : cl => shr_kind_cl
  use lnd_comp_kind    , only : r8 => shr_kind_r8
  use lnd_comp_shr     , only : chkerr

  use mpi

  implicit none
  private

  type fields
     character(len=cl) :: short_name = ""       ! short name
     character(len=cl) :: units = ""            ! unit
     character(len=cl) :: long_name = ""        ! long name
     character(len=cl) :: zaxis = ""            ! name of z axis
     integer           :: nlev                  ! number of layers for 3d fields
  end type fields

  public :: read_initial
  public :: read_restart
  public :: read_static
  public :: read_tiled_file
  public :: write_mosaic_output

  !--------------------------------------------------------------------------
  ! Private module data
  !--------------------------------------------------------------------------

  integer, parameter                :: dbug = 1
  integer, parameter                :: iswater = 17
  integer, parameter                :: max_num_variables = 200
  integer                           :: max_indx = 0
  type(fields)                      :: flds(max_num_variables)
  character(*), parameter           :: modName = "(lnd_comp_io)"

  character(len=*) , parameter :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine read_initial(noahmp, rc)

    ! input/output variables
    type(noahmp_type), intent(inout) :: noahmp
    integer          , intent(inout) :: rc

    ! local variables
    integer                          :: nt
    character(len=cl)                :: filename
    real(ESMF_KIND_R8), pointer      :: ptr(:,:,:)
    type(ESMF_Field)                 :: field
    character(len=*), parameter      :: subname=trim(modName)//':(read_initial) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    if (noahmp%nmlist%ic_type == 'sfc') then
       !----------------------
       ! Set file name for initial conditions
       !----------------------

       filename = trim(noahmp%nmlist%input_dir)//'sfc_data.tile'
       call ESMF_LogWrite(subname//' called for '//trim(filename), ESMF_LOGMSG_INFO)

       !----------------------
       ! Read snow water equivalent
       !----------------------

       call read_tiled_file(filename, 'sheleg', noahmp, field, numrec=1, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       noahmp%init%snow_water_equivalent(:) = ptr(:,1,1)
       nullify(ptr)
       call ESMF_FieldDestroy(field, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       !----------------------
       ! Read snow depth
       !----------------------

       call read_tiled_file(filename, 'snwdph', noahmp, field, numrec=1, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       noahmp%init%snow_depth(:) = ptr(:,1,1)
       nullify(ptr)
       call ESMF_FieldDestroy(field, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       !----------------------
       ! Read canopy surface water
       !----------------------

       call read_tiled_file(filename, 'canopy', noahmp, field, numrec=1, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       noahmp%init%canopy_water(:) = ptr(:,1,1)
       nullify(ptr)
       call ESMF_FieldDestroy(field, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       !----------------------
       ! Read surface skin temperature
       !----------------------

       call read_tiled_file(filename, 'tsea', noahmp, field, numrec=1, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       noahmp%init%skin_temperature(:) = ptr(:,1,1)
       nullify(ptr)
       call ESMF_FieldDestroy(field, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       !----------------------
       ! Read surface soil temperature
       !----------------------

       call read_tiled_file(filename, 'stc', noahmp, field, numrec=1, numlev=noahmp%nmlist%num_soil_levels, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       noahmp%init%soil_temperature(:,:) = ptr(:,:,1)
       nullify(ptr)
       call ESMF_FieldDestroy(field, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       !----------------------
       ! Read surface soil moisture
       !----------------------

       call read_tiled_file(filename, 'smc', noahmp, field, numrec=1, numlev=noahmp%nmlist%num_soil_levels, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       noahmp%init%soil_moisture(:,:) = ptr(:,:,1)
       nullify(ptr)
       call ESMF_FieldDestroy(field, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       !----------------------
       ! Read surface soil liquid
       !----------------------

       call read_tiled_file(filename, 'slc', noahmp, field, numrec=1, numlev=noahmp%nmlist%num_soil_levels, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       noahmp%init%soil_liquid(:,:) = ptr(:,:,1)
       nullify(ptr)
       call ESMF_FieldDestroy(field, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       !----------------------
       ! Set file name for initial conditions
       !----------------------

       write(filename, fmt="(A,I0,A)") trim(noahmp%nmlist%input_dir)//'C', maxval(noahmp%domain%nit), '.initial.tile'
       call ESMF_LogWrite(subname//' called for '//trim(filename), ESMF_LOGMSG_INFO)

       !----------------------
       ! Read snow water equivalent
       !----------------------

       call read_tiled_file(filename, 'snow_water_equivalent', noahmp, field, numrec=1, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       noahmp%init%snow_water_equivalent(:) = ptr(:,1,1)
       nullify(ptr)
       call ESMF_FieldDestroy(field, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       !----------------------
       ! Read snow depth
       !----------------------

       call read_tiled_file(filename, 'snow_depth', noahmp, field, numrec=1, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       noahmp%init%snow_depth(:) = ptr(:,1,1)
       nullify(ptr)
       call ESMF_FieldDestroy(field, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       !----------------------
       ! Read canopy surface water
       !----------------------

       call read_tiled_file(filename, 'canopy_water', noahmp, field, numrec=1, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       noahmp%init%canopy_water(:) = ptr(:,1,1)
       nullify(ptr)
       call ESMF_FieldDestroy(field, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       !----------------------
       ! Read surface skin temperature
       !----------------------

       call read_tiled_file(filename, 'skin_temperature', noahmp, field, numrec=1, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       noahmp%init%skin_temperature(:) = ptr(:,1,1)
       nullify(ptr)
       call ESMF_FieldDestroy(field, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       !----------------------
       ! Read surface soil temperature
       !----------------------

       call read_tiled_file(filename, 'soil_temperature', noahmp, field, numrec=1, numlev=noahmp%nmlist%num_soil_levels, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       noahmp%init%soil_temperature(:,:) = ptr(:,:,1)
       nullify(ptr)
       call ESMF_FieldDestroy(field, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       !----------------------
       ! Read surface soil moisture
       !----------------------

       call read_tiled_file(filename, 'soil_moisture', noahmp, field, numrec=1, numlev=noahmp%nmlist%num_soil_levels, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       noahmp%init%soil_moisture(:,:) = ptr(:,:,1)
       nullify(ptr)
       call ESMF_FieldDestroy(field, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       !----------------------
       ! Read surface soil liquid
       !----------------------

       call read_tiled_file(filename, 'soil_liquid', noahmp, field, numrec=1, numlev=noahmp%nmlist%num_soil_levels, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       noahmp%init%soil_liquid(:,:) = ptr(:,:,1)
       nullify(ptr)
       call ESMF_FieldDestroy(field, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    call ESMF_LogWrite(subname//' done for '//trim(filename), ESMF_LOGMSG_INFO)

  end subroutine read_initial

  !===============================================================================
  subroutine read_restart(noahmp, rc)

    ! input/output variables
    type(noahmp_type), intent(inout) :: noahmp
    integer          , intent(inout) :: rc

    ! local variables
    integer                          :: nt
    character(len=cl)                :: filename
    real(ESMF_KIND_R8), pointer      :: ptr(:,:,:)
    type(ESMF_Field)                 :: field
    character(len=*), parameter      :: subname=trim(modName)//':(read_restart) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    !----------------------
    ! Set file name for restart
    !----------------------

    filename = trim(noahmp%nmlist%restart_dir)//trim(noahmp%nmlist%restart_file)
    call ESMF_LogWrite(subname//' called for '//trim(filename)//'*', ESMF_LOGMSG_INFO)

    !----------------------
    ! zonal wind at lowest model layer
    !----------------------

    call read_tiled_file(filename, 'u1', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%u1(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! meridional wind at lowest model layer
    !----------------------

    call read_tiled_file(filename, 'v1', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%v1(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! soil type
    !----------------------

    call read_tiled_file(filename, 'soiltyp', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%soiltyp(:) = int(ptr(:,1,1))
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! vegetation type
    !----------------------

    call read_tiled_file(filename, 'vegtype', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%vegtype(:) = int(ptr(:,1,1))
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! green vegetation fraction
    !----------------------

    call read_tiled_file(filename, 'sigmaf', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%sigmaf(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! forcing net shortwave flux
    !----------------------

    call read_tiled_file(filename, 'snet', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%snet(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! deep soil temperature
    !----------------------

    call read_tiled_file(filename, 'tg3', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%tg3(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! surface exchange coeff for momentum
    !----------------------

    call read_tiled_file(filename, 'cm', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%cm(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! surface exchange coeff for heat and moisture
    !----------------------

    call read_tiled_file(filename, 'ch', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%ch(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! mean pressure at lowest model layer
    !----------------------

    call read_tiled_file(filename, 'prsl1', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%prsl1(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! dimensionless Exner function at the lowest model layer
    !----------------------

    call read_tiled_file(filename, 'prslk1', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%prslk1(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! Exner function ratio bt midlayer and interface at 1st layer
    !----------------------

    call read_tiled_file(filename, 'prslki', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%prslki(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! dimensionless Exner function at the ground surface
    !----------------------

    call read_tiled_file(filename, 'prsik1', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%prsik1(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! height of bottom layer
    !----------------------

    call read_tiled_file(filename, 'zf', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%zf(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! flag for a point with any land
    !----------------------

    call read_tiled_file(filename, 'dry', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%dry(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! class of sfc slope
    !----------------------

    call read_tiled_file(filename, 'slopetyp', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%slopetyp(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! min fractional coverage of green veg
    !----------------------

    call read_tiled_file(filename, 'shdmin', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%shdmin(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! max fractional coverage of green veg
    !----------------------

    call read_tiled_file(filename, 'shdmax', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%shdmax(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! upper bound on max albedo over deep snow
    !----------------------

    call read_tiled_file(filename, 'snoalb', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%snoalb(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! mean sfc diffuse sw albedo
    !----------------------

    call read_tiled_file(filename, 'sfalb', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%sfalb(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! latitude
    !----------------------

    call read_tiled_file(filename, 'xlatin', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%xlatin(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! cosine of zenith angle
    !----------------------

    call read_tiled_file(filename, 'xcoszin', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%xcoszin(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! microphysics non-convective precipitation
    !----------------------

    call read_tiled_file(filename, 'rainn_mp', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%rainn_mp(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! microphysics convective precipitation
    !----------------------

    call read_tiled_file(filename, 'rainc_mp', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%rainc_mp(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! microphysics snow
    !----------------------

    call read_tiled_file(filename, 'snow_mp', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%snow_mp(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! 
    !----------------------

    call read_tiled_file(filename, 'graupel_mp', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%graupel_mp(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! microphysics ice/hail
    !----------------------

    call read_tiled_file(filename, 'ice_mp', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%ice_mp(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! water equivalent accumulated snow depth
    !----------------------

    call read_tiled_file(filename, 'weasd', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%weasd(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! snow depth (water equiv) over land
    !----------------------

    call read_tiled_file(filename, 'snwdph', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%snwdph(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! ground surface skin temperature
    !----------------------

    call read_tiled_file(filename, 'tskin', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%tskin(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! snow/rain flag for precipitation
    !----------------------

    call read_tiled_file(filename, 'srflag', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%srflag(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! canopy moisture content
    !----------------------

    call read_tiled_file(filename, 'canopy', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%canopy(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! total plant transpiration
    !----------------------

    call read_tiled_file(filename, 'trans', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%trans(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! surface skin temperature (after iteration)
    !----------------------

    call read_tiled_file(filename, 'tsurf', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%tsurf(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! surface roughness
    !----------------------

    call read_tiled_file(filename, 'zorl', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%zorl(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! bulk Richardson number at the surface over land
    !----------------------

    call read_tiled_file(filename, 'rb1', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%rb1(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! Monin-Obukhov similarity function for momentum over land
    !----------------------

    call read_tiled_file(filename, 'fm1', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%fm1(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! Monin-Obukhov similarity function for heat over land
    !----------------------

    call read_tiled_file(filename, 'fh1', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%fh1(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! surface friction velocity over land
    !----------------------

    call read_tiled_file(filename, 'ustar1', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%ustar1(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! surface wind stress over land
    !----------------------

    call read_tiled_file(filename, 'stress1', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%stress1(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! Monin-Obukhov similarity parameter for momentum at 10m over land
    !----------------------

    call read_tiled_file(filename, 'fm101', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%fm101(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! Monin-Obukhov similarity parameter for heat at 2m over land
    !----------------------

    call read_tiled_file(filename, 'fh21', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%fh21(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! actual no. of snow layers
    !----------------------

    call read_tiled_file(filename, 'snowxy', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%snowxy(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! vegetation leaf temperature
    !----------------------

    call read_tiled_file(filename, 'tvxy', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%tvxy(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! bulk ground surface temperature
    !----------------------

    call read_tiled_file(filename, 'tgxy', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%tgxy(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! canopy-intercepted ice
    !----------------------

    call read_tiled_file(filename, 'canicexy', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%canicexy(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! canopy-intercepted liquid water
    !----------------------

    call read_tiled_file(filename, 'canliqxy', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%canliqxy(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! canopy air vapor pressure
    !----------------------

    call read_tiled_file(filename, 'eahxy', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%eahxy(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! canopy air temperature
    !----------------------

    call read_tiled_file(filename, 'tahxy', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%tahxy(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! bulk momentum drag coefficient
    !----------------------

    call read_tiled_file(filename, 'cmxy', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%cmxy(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! bulk sensible heat exchange coefficient
    !----------------------

    call read_tiled_file(filename, 'chxy', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%chxy(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! wetted or snowed fraction of the canopy
    !----------------------

    call read_tiled_file(filename, 'fwetxy', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%fwetxy(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! snow mass at last time step
    !----------------------

    call read_tiled_file(filename, 'sneqvoxy', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%sneqvoxy(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! snow albedo at last time step
    !----------------------

    call read_tiled_file(filename, 'alboldxy', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%alboldxy(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! snowfall on the ground
    !----------------------

    call read_tiled_file(filename, 'qsnowxy', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%qsnowxy(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! lake water storage
    !----------------------

    call read_tiled_file(filename, 'wslakexy', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%wslakexy(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! water table depth
    !----------------------

    call read_tiled_file(filename, 'zwtxy', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%zwtxy(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! water in the aquifer
    !----------------------

    call read_tiled_file(filename, 'waxy', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%waxy(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! groundwater storage
    !----------------------

    call read_tiled_file(filename, 'wtxy', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%wtxy(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! leaf mass
    !----------------------

    call read_tiled_file(filename, 'lfmassxy', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%lfmassxy(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! mass of fine roots
    !----------------------

    call read_tiled_file(filename, 'rtmassxy', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%rtmassxy(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! stem mas
    !----------------------

    call read_tiled_file(filename, 'stmassxy', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%stmassxy(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! mass of wood incl woody roots
    !----------------------

    call read_tiled_file(filename, 'woodxy', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%woodxy(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! stable carbon in deep soil
    !----------------------

    call read_tiled_file(filename, 'stblcpxy', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%stblcpxy(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! short-lived carbon, shallow soil
    !----------------------

    call read_tiled_file(filename, 'fastcpxy', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%fastcpxy(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! leaf area index
    !----------------------

    call read_tiled_file(filename, 'xlaixy', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%xlaixy(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! stem area index
    !----------------------

    call read_tiled_file(filename, 'xsaixy', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%xsaixy(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! snow age factor
    !----------------------

    call read_tiled_file(filename, 'taussxy', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%taussxy(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! soil moisture content in the layer to the water table when deep
    !----------------------

    call read_tiled_file(filename, 'smcwtdxy', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%smcwtdxy(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! recharge to the water table when deep
    !----------------------

    call read_tiled_file(filename, 'deeprechxy', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%deeprechxy(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! recharge to the water table (diagnostic)
    !----------------------

    call read_tiled_file(filename, 'rechxy', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%rechxy(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! albedo - direct visible
    !----------------------

    call read_tiled_file(filename, 'albdvis', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%albdvis(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! albedo - direct NIR
    !----------------------

    call read_tiled_file(filename, 'albdnir', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%albdnir(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! albedo - diffuse visible
    !----------------------

    call read_tiled_file(filename, 'albivis', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%albivis(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! albedo - diffuse NIR
    !----------------------

    call read_tiled_file(filename, 'albinir', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%albinir(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! surface emissivity
    !----------------------

    call read_tiled_file(filename, 'emiss', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%emiss(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! snow cover over land
    !----------------------

    call read_tiled_file(filename, 'sncovr1', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%sncovr1(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! specific humidity at sfc
    !----------------------

    call read_tiled_file(filename, 'qsurf', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%qsurf(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! soil heat flux
    !----------------------

    call read_tiled_file(filename, 'gflux', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%gflux(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! subsurface runoff
    !----------------------

    call read_tiled_file(filename, 'drain', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%drain(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! evaporation from latent heat flux
    !----------------------

    call read_tiled_file(filename, 'evap', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%evap(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! sensible heat flux
    !----------------------

    call read_tiled_file(filename, 'hflx', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%hflx(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! potential evaporation
    !----------------------

    call read_tiled_file(filename, 'ep', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%ep(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! surface runoff
    !----------------------

    call read_tiled_file(filename, 'runoff', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%runoff(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! cm * rho
    !----------------------

    call read_tiled_file(filename, 'cmm', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%cmm(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! ch * rho
    !----------------------

    call read_tiled_file(filename, 'chh', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%chh(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! direct soil evaporation
    !----------------------

    call read_tiled_file(filename, 'evbs', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%evbs(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! canopy water evaporation
    !----------------------

    call read_tiled_file(filename, 'evcw', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%evcw(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! sublimation/deposit from snopack
    !----------------------

    call read_tiled_file(filename, 'sbsno', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%sbsno(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! precipitation advected heat - total
    !----------------------

    call read_tiled_file(filename, 'pah', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%pah(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! evaporation of intercepted water
    !----------------------

    call read_tiled_file(filename, 'ecan', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%ecan(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! transpiration rate
    !----------------------

    call read_tiled_file(filename, 'etran', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%etran(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! soil surface evaporation rate
    !----------------------

    call read_tiled_file(filename, 'edir', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%edir(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! fractional snow cover
    !----------------------

    call read_tiled_file(filename, 'snowc', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%snowc(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! total soil column moisture content
    !----------------------

    call read_tiled_file(filename, 'stm', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%stm(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! snow/freezing-rain latent heat flux
    !----------------------

    call read_tiled_file(filename, 'snohf', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%snohf(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! dry soil moisture threshold
    !----------------------

    call read_tiled_file(filename, 'smcwlt2', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%smcwlt2(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! soil moisture threshold
    !----------------------

    call read_tiled_file(filename, 'smcref2', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%smcref2(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! normalized soil wetness
    !----------------------

    call read_tiled_file(filename, 'wet1', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%wet1(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! combined T2m from tiles
    !----------------------

    call read_tiled_file(filename, 't2mmp', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%t2mmp(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! combined q2m from tiles
    !----------------------

    call read_tiled_file(filename, 'q2mp', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%q2mp(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! function of surface roughness length and green vegetation fraction
    !----------------------

    call read_tiled_file(filename, 'zvfun', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%zvfun(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)

    !----------------------
    ! total soil moisture content
    !----------------------

    call read_tiled_file(filename, 'smc', noahmp, field, numrec=1, numlev=noahmp%nmlist%num_soil_levels, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%smc(:,:) = ptr(:,:,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------
    ! soil temperature
    !----------------------

    call read_tiled_file(filename, 'stc', noahmp, field, numrec=1, numlev=noahmp%nmlist%num_soil_levels, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%stc(:,:) = ptr(:,:,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------
    ! liquid soil moisture
    !----------------------

    call read_tiled_file(filename, 'slc', noahmp, field, numrec=1, numlev=noahmp%nmlist%num_soil_levels, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%slc(:,:) = ptr(:,:,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------
    ! equilibrium soil water content
    !----------------------

    call read_tiled_file(filename, 'smoiseq', noahmp, field, numrec=1, numlev=noahmp%nmlist%num_soil_levels, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%smoiseq(:,:) = ptr(:,:,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------
    ! temperature in surface snow
    !----------------------

    call read_tiled_file(filename, 'tsnoxy', noahmp, field, numrec=1, numlev=abs(noahmp%static%lsnowl)+1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%tsnoxy(:,:) = ptr(:,:,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------
    ! lwe thickness of ice in surface snow
    !----------------------

    call read_tiled_file(filename, 'snicexy', noahmp, field, numrec=1, numlev=abs(noahmp%static%lsnowl)+1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%snicexy(:,:) = ptr(:,:,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------
    ! snow layer liquid water
    !----------------------

    call read_tiled_file(filename, 'snliqxy', noahmp, field, numrec=1, numlev=abs(noahmp%static%lsnowl)+1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%snliqxy(:,:) = ptr(:,:,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------
    ! depth from the top of the snow surface at the bottom of the layer
    !----------------------

    call read_tiled_file(filename, 'zsnsoxy', noahmp, field, numrec=1, numlev=abs(noahmp%static%lsnowl)+noahmp%nmlist%num_soil_levels+1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%zsnsoxy(:,:) = ptr(:,:,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine read_restart

  !===============================================================================
  subroutine read_static(noahmp, rc)

    ! input/output variables
    type(noahmp_type), intent(inout) :: noahmp
    integer          , intent(inout) :: rc

    ! local variables
    integer                          :: nt
    character(len=CL)                :: filename
    real(ESMF_KIND_R8), pointer      :: ptr(:,:,:)
    type(ESMF_Field)                 :: field 
    character(len=*), parameter      :: subname=trim(modName)//':(read_static) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    !----------------------
    ! Set data sources 
    !----------------------

    noahmp%static%isot = 1
    noahmp%static%ivegsrc = 1

    !----------------------
    ! Read latitude 
    !----------------------

    filename = trim(noahmp%nmlist%input_dir)//'oro_data.tile'
    call read_tiled_file(filename, 'geolat', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%xlatin(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------
    ! Read soil type
    !----------------------

    write(filename, fmt="(A,I0,A)") trim(noahmp%nmlist%input_dir)//'C', maxval(noahmp%domain%nit), '.soil_type.tile'
    call read_tiled_file(filename, 'soil_type', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%soiltyp(:) = int(ptr(:,1,1))
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------
    ! Read vegetation type
    !----------------------

    write(filename, fmt="(A,I0,A)") trim(noahmp%nmlist%input_dir)//'C', maxval(noahmp%domain%nit), '.vegetation_type.tile'
    call read_tiled_file(filename, 'vegetation_type', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%vegtype(:) = int(ptr(:,1,1))
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------
    ! Read slope type
    !----------------------

    write(filename, fmt="(A,I0,A)") trim(noahmp%nmlist%input_dir)//'C', maxval(noahmp%domain%nit), '.slope_type.tile'
    call read_tiled_file(filename, 'slope_type', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%slopetyp(:) = int(ptr(:,1,1))
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------
    ! Read deep soil temperature
    !----------------------

    write(filename, fmt="(A,I0,A)") trim(noahmp%nmlist%input_dir)//'C', maxval(noahmp%domain%nit), '.substrate_temperature.tile'
    call read_tiled_file(filename, 'substrate_temperature', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%tg3(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------
    ! Set emissivity
    !----------------------

    ! TODO: this needs to be option in nems.configure
    noahmp%model%emiss(:) = 0.95

    !----------------------
    ! Set albedo
    !----------------------

    ! TODO: this needs to be option in nems.configure
    noahmp%model%alb_monthly(:,:) = 0.25

    !----------------------
    ! Read maximum snow albedo
    !----------------------

    write(filename, fmt="(A,I0,A)") trim(noahmp%nmlist%input_dir)//'C', maxval(noahmp%domain%nit), '.maximum_snow_albedo.tile'
    call read_tiled_file(filename, 'maximum_snow_albedo', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%snoalb(:) = ptr(:,1,1)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------
    ! Read vegetation greenness, monthly average 
    !----------------------

    write(filename, fmt="(A,I0,A)") trim(noahmp%nmlist%input_dir)//'C', maxval(noahmp%domain%nit), '.vegetation_greenness.tile'
    call read_tiled_file(filename, 'vegetation_greenness', noahmp, field, numrec=12, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%gvf_monthly(:,:) = ptr(:,1,:)
    noahmp%model%shdmin(:) = minval(ptr(:,1,:), dim=2)
    noahmp%model%shdmax(:) = maxval(ptr(:,1,:), dim=2)
    nullify(ptr)
    call ESMF_FieldDestroy(field, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------
    ! Set dry 
    !----------------------

    noahmp%model%dry(:) = .false.
    where(noahmp%model%vegtype(:) /= iswater) noahmp%model%dry(:) = .true. 

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine read_static

  !===============================================================================
  subroutine read_tiled_file(filename, varname, noahmp, field, numrec, numlev, rc)

    ! input/output variables
    character(len=*), intent(in)             :: filename
    character(len=*), intent(in)             :: varname 
    type(noahmp_type), intent(inout)         :: noahmp
    type(ESMF_Field), intent(inout)          :: field
    integer, intent(in), optional            :: numrec
    integer, intent(in), optional            :: numlev
    integer, intent(inout), optional         :: rc

    ! local variables
    integer                         :: funit, my_tile, n, i, j, nt, nl
    integer                         :: ndim, nvar, natt, ntime
    integer                         :: isc, iec, jsc, jec
    logical                         :: not_found, is_root_pe
    integer, allocatable            :: dimsizes(:)
    character(len=CL)               :: cname, fname 
    real(ESMF_KIND_R8), pointer     :: ptr(:), ptr3d(:,:,:), ptr4d(:,:,:,:)
    real(ESMF_KIND_R8), allocatable :: rdata(:,:,:,:)
    type(fieldtype), allocatable    :: vars(:)   
    type(ESMF_Field)                :: field_src, field_tmp
    type(ESMF_ArraySpec)            :: arraySpec
    character(len=*), parameter     :: subname=trim(modName)//': (read_tiled_file) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called for '//trim(varname), ESMF_LOGMSG_INFO)

    !----------------------
    ! Define required variables
    !----------------------

    if (present(numrec)) then
       nt = numrec
    else
       nt = 1
    end if

    if (present(numlev)) then
       nl = numlev
    else
       nl = 1
    end if

    my_tile = int(mpp_pe()/(noahmp%domain%layout(1)*noahmp%domain%layout(2)))+1

    is_root_pe = .false.
    if (mpp_pe() == (my_tile-1)*(noahmp%domain%layout(1)*noahmp%domain%layout(2))) is_root_pe = .true.

    !----------------------
    ! Open file and query file attributes
    !----------------------
   
    write(cname, fmt='(A,I1,A)') trim(filename), my_tile, '.nc' 
    call mpp_open(funit, trim(cname), action=MPP_RDONLY, form=MPP_NETCDF, &
                  threading=MPP_MULTI, fileset=MPP_SINGLE, is_root_pe=is_root_pe)
    call mpp_get_info(funit, ndim, nvar, natt, ntime)
    allocate(vars(nvar))
    call mpp_get_fields(funit, vars(:))

    !----------------------
    ! Read requested variable
    !----------------------

    not_found = .true.
    do n = 1, nvar
       ! get variable name
       call mpp_get_atts(vars(n), name=cname)

       ! check variable name
       if (trim(cname) == trim(varname)) then
          ! get array bounds or domain
          call mpp_get_compute_domain(noahmp%domain%mosaic_domain, isc, iec, jsc, jec)

          ! allocate data array and set initial value
          allocate(rdata(isc:iec,jsc:jec,nl,nt))
          rdata(:,:,:,:) = 0.0_r8

          ! read data
          do i = 1, nt
             call mpp_read(funit, vars(n), noahmp%domain%mosaic_domain, rdata, 1)
          end do

          ! set missing values to zero
          where (rdata == 1.0e20)
             rdata(:,:,:,:) = 0.0_r8
          end where
       end if

       not_found = .false.
    end do

    if (not_found) then 
       call mpp_error(FATAL, 'File being read is not the expected one. '//trim(varname)//' is not found.')
    end if

    !----------------------
    ! Move data from ESMF grid to mesh
    !----------------------

    ! set type and rank for ESMF arrayspec
    call ESMF_ArraySpecSet(arraySpec, typekind=ESMF_TYPEKIND_R8, rank=4, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! create source field
    field_src = ESMF_FieldCreate(noahmp%domain%grid, arraySpec, staggerloc=ESMF_STAGGERLOC_CENTER, &
       indexflag=ESMF_INDEX_GLOBAL, ungriddedLBound=(/1,1/), ungriddedUBound=(/nl,nt/), &
       gridToFieldMap=(/1,2/), name=trim(varname), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! get pointer and fill it
    call ESMF_FieldGet(field_src, localDe=0, farrayPtr=ptr4d, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ptr4d(:,:,:,:) = rdata(:,:,:,:)
    nullify(ptr4d)
    if (allocated(rdata)) deallocate(rdata)

    ! create destination field
    field = ESMF_FieldCreate(noahmp%domain%mesh, ESMF_TYPEKIND_R8, name=trim(varname), &
       meshloc=ESMF_MESHLOC_ELEMENT,  ungriddedLbound=(/1,1/), &
       ungriddedUbound=(/nl,nt/), gridToFieldMap=(/1/), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! create routehandle from grid to mesh
    call ESMF_FieldRegridStore(field_src, field, routehandle=noahmp%domain%rh_grid2mesh, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! redist field from ESMF Grid to Mesh
    call ESMF_FieldRedist(field_src, field, noahmp%domain%rh_grid2mesh, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! clean memory
    call ESMF_FieldDestroy(field_src, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------
    ! Output result field for debugging purpose
    !----------------------

    if (dbug > 1) then
       ! TODO: ESMF_FieldWriteVTK() call does not support ungridded dimension
       ! The workaround is implemented in here but it would be nice to extend
       ! ESMF_FieldWriteVTK() call to handle it.  
       field_tmp = ESMF_FieldCreate(noahmp%domain%mesh, ESMF_TYPEKIND_R8, &
          name=trim(varname), meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_FieldGet(field_tmp, localDe=0, farrayPtr=ptr, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr3d, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! write to different file along ungridded dimension
       do i = 1, nl
          do j = 1, nt
            ptr(:) = ptr3d(:,i,j)
            write(fname, fmt='(A,I2.2,A,I2.2)') trim(varname)//'_lev', i, '_time', j
            call ESMF_FieldWriteVTK(field_tmp, trim(fname), rc=rc)
            if (ChkErr(rc,__LINE__,u_FILE_u)) return
          end do
       end do

       ! clean memory
       nullify(ptr)
       nullify(ptr3d)
       call ESMF_FieldDestroy(field_tmp, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine read_tiled_file

  !===============================================================================
  subroutine write_mosaic_output(filename, noahmp, now_time, rc)

    ! input/output variables
    character(len=*)  , intent(in)    :: filename
    type(noahmp_type) , intent(inout) :: noahmp
    real(ESMF_KIND_R8), intent(in)    :: now_time
    integer           , intent(inout) :: rc

    ! local variables
    integer, save                   :: first_time = .true. 
    integer                         :: i, j, id, fid, my_tile
    integer                         :: nx, ny, nz, max_level
    logical                         :: fopen, is_root_pe
    character(len=cl)               :: output_filename, prevar
    real*8, allocatable, save       :: data3d(:,:,:)
    real*8, allocatable, save       :: data4d(:,:,:,:)
    type(axistype)                  :: x, y, z, z1, z2, t
    type(fieldtype)                 :: f, f1, f2
    type(domain1D)                  :: xdom, ydom
    type(ESMF_ArraySpec)            :: arraySpec
    type(ESMF_Field), save          :: field_grid, field_mesh
    real(ESMF_KIND_R8), pointer     :: ptr2d(:,:)
    real(ESMF_KIND_R8), pointer     :: ptr3d(:,:,:)
    character(len=*), parameter     :: subname=trim(modName)//':(write_mosaic_output) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    !----------------------
    ! Define required variables
    !----------------------

    my_tile = int(mpp_pe()/(noahmp%domain%layout(1)*noahmp%domain%layout(2)))+1

    is_root_pe = .false.
    if (mpp_pe() == (my_tile-1)*(noahmp%domain%layout(1)*noahmp%domain%layout(2))) is_root_pe = .true.

    !----------------------
    ! open file
    !----------------------

    ! check unit
    fid = 7
    do
       inquire(unit=fid, opened=fopen)
       if (.not. fopen) exit
       fid = fid+1
       if (fid .eq. 100) call mpp_error(FATAL, 'Unable to locate unit number.')
    end do

    ! open file
    write(output_filename, fmt='(a,i1)') trim(filename)//'.tile', my_tile
    call mpp_open(fid, output_filename, action=MPP_OVERWR, form=MPP_NETCDF, &
                  threading=MPP_SINGLE, domain=noahmp%domain%mosaic_domain, &
                  is_root_pe=is_root_pe)

    !----------------------
    ! add global attributes to file, provenance information
    !----------------------

    call mpp_write_meta(fid, 'delt'     , rval=real(noahmp%static%delt))
    call mpp_write_meta(fid, 'idveg'    , ival=noahmp%static%idveg)
    call mpp_write_meta(fid, 'iopt_crs' , ival=noahmp%static%iopt_crs)
    call mpp_write_meta(fid, 'iopt_btr' , ival=noahmp%static%iopt_btr)
    call mpp_write_meta(fid, 'iopt_run' , ival=noahmp%static%iopt_run)
    call mpp_write_meta(fid, 'iopt_sfc' , ival=noahmp%static%iopt_sfc)
    call mpp_write_meta(fid, 'iopt_frz' , ival=noahmp%static%iopt_frz)
    call mpp_write_meta(fid, 'iopt_inf' , ival=noahmp%static%iopt_inf)
    call mpp_write_meta(fid, 'iopt_rad' , ival=noahmp%static%iopt_rad)
    call mpp_write_meta(fid, 'iopt_alb' , ival=noahmp%static%iopt_alb)
    call mpp_write_meta(fid, 'iopt_snf' , ival=noahmp%static%iopt_snf)
    call mpp_write_meta(fid, 'iopt_tbot', ival=noahmp%static%iopt_tbot)
    call mpp_write_meta(fid, 'iopt_stc' , ival=noahmp%static%iopt_stc)
    call mpp_write_meta(fid, 'iopt_trs' , ival=noahmp%static%iopt_trs)

    !----------------------
    ! define dimensions
    !----------------------

    ! query domain to get x and y components
    call mpp_get_domain_components(noahmp%domain%mosaic_domain, xdom, ydom)

    ! x-axis
    nx = noahmp%domain%nit(my_tile)
    call mpp_write_meta(fid, x, 'xc', 'unitless', 'x-coordinate', cartesian='X', domain=xdom, data=(/(i*1.0,i=1,nx)/))

    ! y-axis
    ny = noahmp%domain%njt(my_tile)
    call mpp_write_meta(fid, y, 'yc', 'unitless', 'y-coordinate', cartesian='Y', domain=ydom, data=(/(i*1.0,i=1,ny)/))

    ! z-axises, soil and snow layers
    call mpp_write_meta(fid, z, 'soil_levels', 'meters', 'soil levels', data=real(noahmp%nmlist%soil_level_nodes))
    call mpp_write_meta(fid, z1,'snow_levels', 'unitless', 'snow_levels', data=(/(i*1.0,i=noahmp%static%lsnowl,0)/))
    call mpp_write_meta(fid, z2,'snso_levels', 'unitless', 'snso_levels', data=(/(i*1.0,i=noahmp%static%lsnowl,noahmp%nmlist%num_soil_levels)/))

    ! time axis
    call mpp_write_meta(fid, t, 'time', "seconds since "//noahmp%model%reference_date, 'time', cartesian='T')

    ! write them to the file
    call mpp_write(fid, x)
    call mpp_write(fid, y)
    call mpp_write(fid, z)
    call mpp_write(fid, z1)
    call mpp_write(fid, z2)

    !----------------------
    ! define fixed variables
    !----------------------

    call mpp_write_meta(fid, f, (/x,y/), 'grid_xt', 'degrees_E', 'T-cell longitude', pack=1)
    call mpp_write(fid, f, noahmp%domain%mosaic_domain, noahmp%domain%lont)
    call mpp_write_meta(fid, f, (/x,y/), 'grid_yt', 'degrees_N', 'T-cell latitude', pack=1)
    call mpp_write(fid, f, noahmp%domain%mosaic_domain, noahmp%domain%latt)

    !----------------------
    ! create routehandle to redist selected fields from mesh to mosaic grid
    !----------------------

    if (first_time) then
       ! set field type
       call ESMF_ArraySpecSet(arraySpec, typekind=ESMF_TYPEKIND_R8, rank=3, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! create field in mesh side
       field_mesh = ESMF_FieldCreate(noahmp%domain%mesh, ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, &
          ungriddedLbound=(/1/), ungriddedUbound=(/max_num_variables/), gridToFieldMap=(/1/), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! create field in grid side
       field_grid = ESMF_FieldCreate(noahmp%domain%grid, arraySpec, staggerloc=ESMF_STAGGERLOC_CENTER, &
          indexflag=ESMF_INDEX_GLOBAL, ungriddedLBound=(/1/), ungriddedUBound=(/max_num_variables/), gridToFieldMap=(/1,2/), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! create routehandle from mesh to grid
       call ESMF_FieldRegridStore(field_mesh, field_grid, routehandle=noahmp%domain%rh_mesh2grid_r8, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    !----------------------
    ! get pionter from field and construct table for fields
    !----------------------

    call ESMF_FieldGet(field_mesh, localDe=0, farrayPtr=ptr2d, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call fld_add("ps"        ,"surface pressure"                                                  ,"Pa"     ,ptr2d=ptr2d, v1r8=noahmp%forc%ps)
    call fld_add("u1"        ,"u-component of wind"                                               ,"m/s"    ,ptr2d=ptr2d, v1r8=noahmp%model%u1)
    call fld_add("v1"        ,"v-component of wind"                                               ,"m/s"    ,ptr2d=ptr2d, v1r8=noahmp%model%v1)
    call fld_add("t1"        ,"forcing air temperature"                                           ,"K"      ,ptr2d=ptr2d, v1r8=noahmp%forc%t1)
    call fld_add("q1"        ,"forcing specific humidity"                                         ,"kg/kg"  ,ptr2d=ptr2d, v1r8=noahmp%forc%q1)
    call fld_add("soiltyp"   ,"soil type"                                                         ,"1"      ,ptr2d=ptr2d, v1i4=noahmp%model%soiltyp)
    call fld_add("vegtype"   ,"vegetation type"                                                   ,"1"      ,ptr2d=ptr2d, v1i4=noahmp%model%vegtype)
    call fld_add("sigmaf"    ,"green vegetation fraction"                                         ,"1"      ,ptr2d=ptr2d, v1r8=noahmp%model%sigmaf) 
    call fld_add("dlwflx"    ,"forcing longwave downward flux"                                    ,"W/m2"   ,ptr2d=ptr2d, v1r8=noahmp%forc%dlwflx)
    call fld_add("dswsfc"    ,"forcing shortwave downward flux"                                   ,"W/m2"   ,ptr2d=ptr2d, v1r8=noahmp%forc%dswsfc)
    call fld_add("snet"      ,"forcing net shortwave flux"                                        ,"W/m2"   ,ptr2d=ptr2d, v1r8=noahmp%model%snet) 
    call fld_add("tg3"       ,"deep soil temperature"                                             ,"K"      ,ptr2d=ptr2d, v1r8=noahmp%model%tg3)
    call fld_add("cm"        ,"surface exchange coeff for momentum"                               ,"m/s"    ,ptr2d=ptr2d, v1r8=noahmp%model%cm)
    call fld_add("ch"        ,"surface exchange coeff for heat and moisture"                      ,"m/s"    ,ptr2d=ptr2d, v1r8=noahmp%model%ch) 
    call fld_add("prsl1"     ,"mean pressure at lowest model layer"                               ,"Pa"     ,ptr2d=ptr2d, v1r8=noahmp%model%prsl1)
    call fld_add("prslk1"    ,"dimensionless Exner function at the lowest model layer"            ,"1"      ,ptr2d=ptr2d, v1r8=noahmp%model%prslk1)
    call fld_add("prslki"    ,"Exner function ratio bt midlayer and interface at 1st layer"       ,"1"      ,ptr2d=ptr2d, v1r8=noahmp%model%prslki)
    call fld_add("prsik1"    ,"dimensionless Exner function at the ground surface"                ,"1"      ,ptr2d=ptr2d, v1r8=noahmp%model%prsik1)
    call fld_add("zf"        ,"height of bottom layer"                                            ,"m"      ,ptr2d=ptr2d, v1r8=noahmp%model%zf)
    call fld_add("dry"       ,"flag for a point with any land"                                    ,"1"      ,ptr2d=ptr2d, v1l =noahmp%model%dry)  
    call fld_add("wind"      ,"wind speed"                                                        ,"m/s"    ,ptr2d=ptr2d, v1r8=noahmp%forc%wind)
    call fld_add("slopetyp"  ,"class of sfc slope"                                                ,"1"      ,ptr2d=ptr2d, v1i4=noahmp%model%slopetyp)
    call fld_add("shdmin"    ,"min fractional coverage of green veg"                              ,"1"      ,ptr2d=ptr2d, v1r8=noahmp%model%shdmin)
    call fld_add("shdmax"    ,"max fractional coverage of green veg"                              ,"1"      ,ptr2d=ptr2d, v1r8=noahmp%model%shdmax)
    call fld_add("snoalb"    ,"upper bound on max albedo over deep snow"                          ,"1"      ,ptr2d=ptr2d, v1r8=noahmp%model%snoalb)
    call fld_add("sfalb"     ,"mean sfc diffuse sw albedo"                                        ,"1"      ,ptr2d=ptr2d, v1r8=noahmp%model%sfalb)
    call fld_add("xlatin"    ,"latitude"                                                          ,"radian" ,ptr2d=ptr2d, v1r8=noahmp%model%xlatin)
    call fld_add("xcoszin"   ,"cosine of zenith angle"                                            ,"degree" ,ptr2d=ptr2d, v1r8=noahmp%model%xcoszin)
    call fld_add("garea"     ,"area of the grid cell"                                             ,"m2"     ,ptr2d=ptr2d, v1r8=noahmp%domain%garea)
    call fld_add("rainn_mp"  ,"microphysics non-convective precipitation"                         ,"mm"     ,ptr2d=ptr2d, v1r8=noahmp%model%rainn_mp)
    call fld_add("rainc_mp"  ,"microphysics convective precipitation"                             ,"mm"     ,ptr2d=ptr2d, v1r8=noahmp%model%rainc_mp)
    call fld_add("snow_mp"   ,"microphysics snow"                                                 ,"mm"     ,ptr2d=ptr2d, v1r8=noahmp%model%snow_mp)
    call fld_add("graupel_mp","microphysics graupel"                                              ,"mm"     ,ptr2d=ptr2d, v1r8=noahmp%model%graupel_mp)
    call fld_add("ice_mp"    ,"microphysics ice/hail"                                             ,"mm"     ,ptr2d=ptr2d, v1r8=noahmp%model%ice_mp)
    call fld_add("weasd"     ,"water equivalent accumulated snow depth"                           ,"mm"     ,ptr2d=ptr2d, v1r8=noahmp%model%weasd)
    call fld_add("snwdph"    ,"snow depth (water equiv) over land"                                ,"m"      ,ptr2d=ptr2d, v1r8=noahmp%model%snwdph)
    call fld_add("tskin"     ,"ground surface skin temperature"                                   ,"K"      ,ptr2d=ptr2d, v1r8=noahmp%model%tskin)
    call fld_add("tprcp"     ,"total precipitation"                                               ,"mm"     ,ptr2d=ptr2d, v1r8=noahmp%forc%tprcp)
    call fld_add("srflag"    ,"snow/rain flag for precipitation"                                  ,"1"      ,ptr2d=ptr2d, v1r8=noahmp%model%srflag)
    call fld_add("smc"       ,"total soil moisture content"                                       ,"m3/m3"  ,ptr2d=ptr2d, v2r8=noahmp%model%smc, zaxis="z")
    call fld_add("stc"       ,"soil temperature"                                                  ,"K"      ,ptr2d=ptr2d, v2r8=noahmp%model%stc, zaxis="z")
    call fld_add("slc"       ,"liquid soil moisture"                                              ,"m3/m3"  ,ptr2d=ptr2d, v2r8=noahmp%model%slc, zaxis="z")
    call fld_add("canopy"    ,"canopy moisture content"                                           ,"m"      ,ptr2d=ptr2d, v1r8=noahmp%model%canopy)
    call fld_add("trans"     ,"total plant transpiration"                                         ,"m/2"    ,ptr2d=ptr2d, v1r8=noahmp%model%trans)
    call fld_add("tsurf"     ,"surface skin temperature (after iteration)"                        ,"K"      ,ptr2d=ptr2d, v1r8=noahmp%model%tsurf)
    call fld_add("zorl"      ,"surface roughness"                                                 ,"m"      ,ptr2d=ptr2d, v1r8=noahmp%model%zorl)
    call fld_add("rb1"       ,"bulk Richardson number at the surface over land"                   ,"1"      ,ptr2d=ptr2d, v1r8=noahmp%model%rb1)
    call fld_add("fm1"       ,"Monin-Obukhov similarity function for momentum over land"          ,"1"      ,ptr2d=ptr2d, v1r8=noahmp%model%fm1)
    call fld_add("fh1"       ,"Monin-Obukhov similarity function for heat over land"              ,"1"      ,ptr2d=ptr2d, v1r8=noahmp%model%fh1)
    call fld_add("ustar1"    ,"surface friction velocity over land"                               ,"m/s"    ,ptr2d=ptr2d, v1r8=noahmp%model%ustar1)
    call fld_add("stress1"   ,"surface wind stress over land"                                     ,"m2/s2"  ,ptr2d=ptr2d, v1r8=noahmp%model%stress1)
    call fld_add("fm101"     ,"Monin-Obukhov similarity parameter for momentum at 10m over land"  ,"1"      ,ptr2d=ptr2d, v1r8=noahmp%model%fm101)
    call fld_add("fh21"      ,"Monin-Obukhov similarity parameter for heat at 2m over land"       ,"1"      ,ptr2d=ptr2d, v1r8=noahmp%model%fh21)
    call fld_add("snowxy"    ,"actual no. of snow layers"                                         ,"1"      ,ptr2d=ptr2d, v1r8=noahmp%model%snowxy)
    call fld_add("tvxy"      ,"vegetation leaf temperature"                                       ,"K"      ,ptr2d=ptr2d, v1r8=noahmp%model%tvxy)
    call fld_add("tgxy"      ,"bulk ground surface temperature"                                   ,"K"      ,ptr2d=ptr2d, v1r8=noahmp%model%tgxy)
    call fld_add("canicexy"  ,"canopy-intercepted ice"                                            ,"mm"     ,ptr2d=ptr2d, v1r8=noahmp%model%canicexy)
    call fld_add("canliqxy"  ,"canopy-intercepted liquid water"                                   ,"mm"     ,ptr2d=ptr2d, v1r8=noahmp%model%canliqxy)
    call fld_add("eahxy"     ,"canopy air vapor pressure"                                         ,"Pa"     ,ptr2d=ptr2d, v1r8=noahmp%model%eahxy)
    call fld_add("tahxy"     ,"canopy air temperature"                                            ,"K"      ,ptr2d=ptr2d, v1r8=noahmp%model%tahxy)
    call fld_add("cmxy"      ,"bulk momentum drag coefficient"                                    ,"m/s"    ,ptr2d=ptr2d, v1r8=noahmp%model%cmxy)
    call fld_add("chxy"      ,"bulk sensible heat exchange coefficient"                           ,"m/s"    ,ptr2d=ptr2d, v1r8=noahmp%model%chxy)
    call fld_add("fwetxy"    ,"wetted or snowed fraction of the canopy"                           ,"1"      ,ptr2d=ptr2d, v1r8=noahmp%model%fwetxy)
    call fld_add("sneqvoxy"  ,"snow mass at last time step"                                       ,"mm"     ,ptr2d=ptr2d, v1r8=noahmp%model%sneqvoxy)
    call fld_add("alboldxy"  ,"snow albedo at last time step"                                     ,"1"      ,ptr2d=ptr2d, v1r8=noahmp%model%alboldxy)
    call fld_add("qsnowxy"   ,"snowfall on the ground"                                            ,"mm/s"   ,ptr2d=ptr2d, v1r8=noahmp%model%qsnowxy)
    call fld_add("wslakexy"  ,"lake water storage"                                                ,"mm"     ,ptr2d=ptr2d, v1r8=noahmp%model%wslakexy)
    call fld_add("zwtxy"     ,"water table depth"                                                 ,"m"      ,ptr2d=ptr2d, v1r8=noahmp%model%zwtxy)
    call fld_add("waxy"      ,"water in the aquifer"                                              ,"mm"     ,ptr2d=ptr2d, v1r8=noahmp%model%waxy)
    call fld_add("wtxy"      ,"groundwater storage"                                               ,"mm"     ,ptr2d=ptr2d, v1r8=noahmp%model%wtxy)
    call fld_add("tsnoxy"    ,"temperature in surface snow"                                       ,"K"      ,ptr2d=ptr2d, v2r8=noahmp%model%tsnoxy , zaxis="z1")
    call fld_add("zsnsoxy"   ,"depth from the top of the snow surface at the bottom of the layer" ,"m"      ,ptr2d=ptr2d, v2r8=noahmp%model%zsnsoxy, zaxis="z2")
    call fld_add("snicexy"   ,"lwe thickness of ice in surface snow"                              ,"mm"     ,ptr2d=ptr2d, v2r8=noahmp%model%snicexy, zaxis="z1")
    call fld_add("snliqxy"   ,"snow layer liquid water"                                           ,"mm"     ,ptr2d=ptr2d, v2r8=noahmp%model%snliqxy, zaxis="z1")
    call fld_add("lfmassxy"  ,"leaf mass"                                                         ,"g/m2"   ,ptr2d=ptr2d, v1r8=noahmp%model%lfmassxy)
    call fld_add("rtmassxy"  ,"mass of fine roots"                                                ,"g/m2"   ,ptr2d=ptr2d, v1r8=noahmp%model%rtmassxy)
    call fld_add("stmassxy"  ,"stem mas"                                                          ,"g/m2"   ,ptr2d=ptr2d, v1r8=noahmp%model%stmassxy)
    call fld_add("woodxy"    ,"mass of wood incl woody roots"                                     ,"g/m2"   ,ptr2d=ptr2d, v1r8=noahmp%model%woodxy)
    call fld_add("stblcpxy"  ,"stable carbon in deep soil"                                        ,"g/m2"   ,ptr2d=ptr2d, v1r8=noahmp%model%stblcpxy)
    call fld_add("fastcpxy"  ,"short-lived carbon, shallow soil"                                  ,"g/m2"   ,ptr2d=ptr2d, v1r8=noahmp%model%fastcpxy)
    call fld_add("xlaixy"    ,"leaf area index"                                                   ,"1"      ,ptr2d=ptr2d, v1r8=noahmp%model%xlaixy)
    call fld_add("xsaixy"    ,"stem area index"                                                   ,"1"      ,ptr2d=ptr2d, v1r8=noahmp%model%xsaixy)
    call fld_add("taussxy"   ,"snow age factor"                                                   ,"1"      ,ptr2d=ptr2d, v1r8=noahmp%model%taussxy)
    call fld_add("smoiseq"   ,"equilibrium soil water content"                                    ,"m3/m3"  ,ptr2d=ptr2d, v2r8=noahmp%model%smoiseq, zaxis="z")
    call fld_add("smcwtdxy"  ,"soil moisture content in the layer to the water table when deep"   ,"mm"     ,ptr2d=ptr2d, v1r8=noahmp%model%smcwtdxy)
    call fld_add("deeprechxy","recharge to the water table when deep"                             ,"1"      ,ptr2d=ptr2d, v1r8=noahmp%model%deeprechxy)
    call fld_add("rechxy"    ,"recharge to the water table (diagnostic)"                          ,"1"      ,ptr2d=ptr2d, v1r8=noahmp%model%rechxy)
    call fld_add("albdvis"   ,"albedo - direct visible"                                           ,"1"      ,ptr2d=ptr2d, v1r8=noahmp%model%albdvis)
    call fld_add("albdnir"   ,"albedo - direct NIR"                                               ,"1"      ,ptr2d=ptr2d, v1r8=noahmp%model%albdnir)
    call fld_add("albivis"   ,"albedo - diffuse visible"                                          ,"1"      ,ptr2d=ptr2d, v1r8=noahmp%model%albivis)
    call fld_add("albinir"   ,"albedo - diffuse NIR"                                              ,"1"      ,ptr2d=ptr2d, v1r8=noahmp%model%albinir)
    call fld_add("emiss"     ,"surface emissivity"                                                ,"1"      ,ptr2d=ptr2d, v1r8=noahmp%model%emiss)
    call fld_add("sncovr1"   ,"snow cover over land"                                              ,"1"      ,ptr2d=ptr2d, v1r8=noahmp%model%sncovr1)
    call fld_add("qsurf"     ,"specific humidity at sfc"                                          ,"kg/kg"  ,ptr2d=ptr2d, v1r8=noahmp%model%qsurf)
    call fld_add("gflux"     ,"soil heat flux"                                                    ,"W/m2"   ,ptr2d=ptr2d, v1r8=noahmp%model%gflux)
    call fld_add("drain"     ,"subsurface runoff"                                                 ,"mm/s"   ,ptr2d=ptr2d, v1r8=noahmp%model%drain)
    call fld_add("evap"      ,"evaporation from latent heat flux"                                 ,"mm/s"   ,ptr2d=ptr2d, v1r8=noahmp%model%evap)
    call fld_add("hflx"      ,"sensible heat flux"                                                ,"W/m2"   ,ptr2d=ptr2d, v1r8=noahmp%model%hflx)
    call fld_add("ep"        ,"potential evaporation"                                             ,"W/m2"   ,ptr2d=ptr2d, v1r8=noahmp%model%ep)
    call fld_add("runoff"    ,"surface runoff"                                                    ,"m/s"    ,ptr2d=ptr2d, v1r8=noahmp%model%runoff)
    call fld_add("cmm"       ,"cm * rho"                                                          ,"m/s"    ,ptr2d=ptr2d, v1r8=noahmp%model%cmm)
    call fld_add("chh"       ,"ch * rho"                                                          ,"kg/m2/s",ptr2d=ptr2d, v1r8=noahmp%model%chh)
    call fld_add("evbs"      ,"direct soil evaporation"                                           ,"m/s"    ,ptr2d=ptr2d, v1r8=noahmp%model%evbs)
    call fld_add("evcw"      ,"canopy water evaporation"                                          ,"m/s"    ,ptr2d=ptr2d, v1r8=noahmp%model%evcw)
    call fld_add("sbsno"     ,"sublimation/deposit from snopack"                                  ,"m/s"    ,ptr2d=ptr2d, v1r8=noahmp%model%sbsno)
    call fld_add("pah"       ,"precipitation advected heat - total"                               ,"W/m2"   ,ptr2d=ptr2d, v1r8=noahmp%model%pah)
    call fld_add("ecan"      ,"evaporation of intercepted water"                                  ,"kg/m2/s",ptr2d=ptr2d, v1r8=noahmp%model%ecan)
    call fld_add("etran"     ,"transpiration rate"                                                ,"kg/m2/s",ptr2d=ptr2d, v1r8=noahmp%model%etran)
    call fld_add("edir"      ,"soil surface evaporation rate"                                     ,"kg/m2/s",ptr2d=ptr2d, v1r8=noahmp%model%edir)
    call fld_add("snowc"     ,"fractional snow cover"                                             ,"1"      ,ptr2d=ptr2d, v1r8=noahmp%model%snowc)
    call fld_add("stm"       ,"total soil column moisture content"                                ,"m"      ,ptr2d=ptr2d, v1r8=noahmp%model%stm)
    call fld_add("snohf"     ,"snow/freezing-rain latent heat flux"                               ,"W/m2"   ,ptr2d=ptr2d, v1r8=noahmp%model%snohf)
    call fld_add("smcwlt2"   ,"dry soil moisture threshold"                                       ,"m3/m3"  ,ptr2d=ptr2d, v1r8=noahmp%model%smcwlt2)
    call fld_add("smcref2"   ,"soil moisture threshold"                                           ,"m3/m3"  ,ptr2d=ptr2d, v1r8=noahmp%model%smcref2)
    call fld_add("wet1"      ,"normalized soil wetness"                                           ,"1"      ,ptr2d=ptr2d, v1r8=noahmp%model%wet1)
    call fld_add("t2mmp"     ,"combined T2m from tiles"                                           ,"K"      ,ptr2d=ptr2d, v1r8=noahmp%model%t2mmp)
    call fld_add("q2mp"      ,"combined q2m from tiles"                                           ,"kg/kg"  ,ptr2d=ptr2d, v1r8=noahmp%model%q2mp)
    call fld_add("zvfun"     ,"function of surface roughness length and green vegetation fraction","1"      ,ptr2d=ptr2d, v1r8=noahmp%model%zvfun)

    !----------------------
    ! masked out data over ocean/inland water/lake
    !----------------------

    do i = 1, max_indx
       where(noahmp%domain%mask < 1) ptr2d(:,i) = 1.0e20
    end do
    nullify(ptr2d)

    !----------------------
    ! redist from mesh to grid and extract pointer to write out
    !----------------------

    call ESMF_FieldRedist(field_mesh, field_grid, noahmp%domain%rh_mesh2grid_r8, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_FieldGet(field_grid, localDe=0, farrayPtr=ptr3d, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------
    ! define temporal variables that will be used to write data through FMS interface
    !----------------------

    if (first_time) then
       ! this is one time allocate, no need to deallocate it
       allocate(data3d(lbound(ptr3d, dim=1):ubound(ptr3d, dim=1),lbound(ptr3d, dim=2):ubound(ptr3d, dim=2),1))
       max_level = noahmp%nmlist%num_soil_levels-noahmp%static%lsnowl+1
       allocate(data4d(lbound(ptr3d, dim=1):ubound(ptr3d, dim=1),lbound(ptr3d, dim=2):ubound(ptr3d, dim=2),max_level,1))
       first_time = .false.
    end if

    prevar = ""
    do i = 1, max_indx
       ! 2d fields with x, y
       if (flds(i)%nlev == 0) then
          ! define variable
          ! TODO: pack = 1 writes double. without pack is writes float by default. make it configurable.
          call mpp_write_meta(fid, f1, (/y,x,t/), trim(flds(i)%short_name), trim(flds(i)%units), trim(flds(i)%long_name), missing=1.0e20, pack=1)
          ! put data to temporary variable
          data3d(:,:,1) = ptr3d(:,:,i)
          ! write to file
          call mpp_write(fid, f1, noahmp%domain%mosaic_domain, data3d, now_time)
       ! 3d fields with x, y, z/z1/z2
       else
          ! control to skip writing data multiple times
          if (trim(prevar) /= trim(flds(i)%short_name)) then
             ! define variable
             if (trim(flds(i)%zaxis) == "z") then
                call mpp_write_meta(fid, f2, (/y,x,z,t/), trim(flds(i)%short_name), trim(flds(i)%units), trim(flds(i)%long_name), missing=1.0e20, pack=1)
             else if (trim(flds(i)%zaxis) == "z1") then
                call mpp_write_meta(fid, f2, (/y,x,z1,t/), trim(flds(i)%short_name), trim(flds(i)%units), trim(flds(i)%long_name), missing=1.0e20, pack=1)
             else if (trim(flds(i)%zaxis) == "z2") then
                call mpp_write_meta(fid, f2, (/y,x,z2,t/), trim(flds(i)%short_name), trim(flds(i)%units), trim(flds(i)%long_name), missing=1.0e20, pack=1)
             else
                call mpp_error(FATAL, 'zaxis can be z, z1 or z2. '//trim(flds(i)%zaxis)//' not recognized for '//trim(flds(i)%short_name))
             end if
             ! put data to temporary variable
             data4d(:,:,1:flds(i)%nlev,1) = ptr3d(:,:,i:i+flds(i)%nlev-1)
             ! write to file
             call mpp_write(fid, f2, noahmp%domain%mosaic_domain, data4d(:,:,1:flds(i)%nlev,:), now_time)
          end if
       end if
       prevar = trim(flds(i)%short_name)
    end do
    nullify(ptr3d)

    !----------------------
    ! close file
    !----------------------

    call mpp_close(fid)

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine write_mosaic_output

  !===============================================================================
  subroutine fld_add(short_name, variable_unit, long_name, ptr2d, v1r8, v1i4, v1l, v2r8, zaxis)

    ! input/output variables
    character(len=*),  intent(in)              :: short_name
    character(len=*),  intent(in)              :: variable_unit
    character(len=*),  intent(in)              :: long_name
    real(r8), pointer, optional, intent(inout) :: ptr2d(:,:)
    real(r8),          optional, intent(in)    :: v1r8(:)
    integer,           optional, intent(in)    :: v1i4(:)
    logical,           optional, intent(in)    :: v1l(:)
    real(r8),          optional, intent(in)    :: v2r8(:,:)
    character(len=*),  optional, intent(in)    :: zaxis

    ! local variables
    integer                     :: i, indx
    logical                     :: found, restart
    character(len=*), parameter :: subname=trim(modName)//': (fld_add) '
    !-------------------------------------------------------------------------------

    ! find out indices
    indx = 0
    found = .false.
    do i = 1, max_num_variables
       if (trim(flds(i)%short_name) == trim(short_name)) then
          indx = i
          found = .true.
          exit
       end if
    end do

    ! if it is a new entry, increment max_indx
    if (.not. found) then
       indx = max_indx+1
       if (present(v2r8)) then
          ! variables with vertical dimension, add each layer as a seperate variable
          max_indx = max_indx+size(v2r8, dim=2)
       else
          ! variables without vertical dimension
          max_indx = max_indx+1
       end if
       if (max_indx > max_num_variables) then
          print*, "max_indx > max_num_variables could not add more variable! increase max_num_variables ..."
          return
       end if
    end if

    ! add field metadata and fill pointer with data
    ! NOTE: if data is not allocated then all present statements return .false. and does not set flds(indx)
    ! 2d variables
    if (present(v1r8) .or. present(v1i4) .or. present(v1l)) then
       flds(indx)%short_name = trim(short_name)
       flds(indx)%units = trim(variable_unit)
       flds(indx)%long_name = trim(long_name)
       flds(indx)%nlev = 0
       if (present(v1r8)) then
          ptr2d(:,indx) = v1r8(:)    
       else if (present(v1i4)) then
          ptr2d(:,indx) = v1i4(:)
       else if (present(v1l)) then
          ptr2d(:,indx) = v1l(:)
       end if
    ! 3d variables
    else if (present(v2r8)) then
       ! add each level as seperate variable
       do i = 1, size(v2r8, dim=2)
          flds(indx+i-1)%short_name = trim(short_name)
          flds(indx+i-1)%units = trim(variable_unit)
          flds(indx+i-1)%long_name = trim(long_name)
          flds(indx+i-1)%zaxis = trim(zaxis)
          flds(indx+i-1)%nlev = size(v2r8, dim=2)
          ptr2d(:,indx+i-1) = v2r8(:,lbound(v2r8, dim=2)+i-1)
       end do
    end if
 
  end subroutine fld_add  

end module lnd_comp_io
