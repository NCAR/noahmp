module lnd_comp_domain

  ! This file contains land domain related routines 

  use ESMF ,  only : operator(/=)
  use ESMF ,  only : ESMF_GridComp, ESMF_GridCompGet, ESMF_Grid, ESMF_VM
  use ESMF ,  only : ESMF_Mesh, ESMF_MeshCreate, ESMF_MeshGet
  use ESMF ,  only : ESMF_Field, ESMF_Array, ESMF_FILEFORMAT_ESMFMESH
  use ESMF ,  only : ESMF_DistGrid, ESMF_DistGridGet, ESMF_Decomp_Flag
  use ESMF ,  only : ESMF_KIND_R4, ESMF_TYPEKIND_R8, ESMF_SUCCESS
  use ESMF ,  only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_FieldGet
  use ESMF ,  only : ESMF_VMGet, ESMF_VMReduce, ESMF_VMBroadcast
  use ESMF ,  only : ESMF_RouteHandle, ESMF_REDUCE_SUM, ESMF_FieldCreate
  use ESMF ,  only : ESMF_FieldRegridStore, ESMF_REGRIDMETHOD_CONSERVE
  use ESMF ,  only : ESMF_NORMTYPE_DSTAREA, ESMF_UNMAPPEDACTION_IGNORE 
  use ESMF ,  only : ESMF_MESHLOC_ELEMENT, ESMF_ArrayCreate
  use ESMF ,  only : ESMF_INDEX_GLOBAL, ESMF_GRIDITEM_MASK, ESMF_GridAddItem
  use ESMF ,  only : ESMF_FieldRegrid, ESMF_TERMORDER_SRCSEQ, ESMF_FieldDestroy
  use ESMF ,  only : ESMF_REGION_TOTAL, ESMF_VMAllReduce, ESMF_MeshSet
  use ESMF ,  only : ESMF_LogWrite, ESMF_LOGMSG_ERROR, ESMF_FAILURE
  use ESMF ,  only : ESMF_DECOMP_SYMMEDGEMAX, ESMF_GridCreateMosaic
  use ESMF ,  only : ESMF_STAGGERLOC_CENTER, ESMF_STAGGERLOC_CORNER
  use ESMF ,  only : ESMF_RouteHandleDestroy, ESMF_GridGet, ESMF_GridGetCoord
  use ESMF ,  only : ESMF_FieldRegridGetArea, ESMF_CoordSys_Flag
  use ESMF ,  only : ESMF_MeshGetFieldBounds, ESMF_COORDSYS_CART, ESMF_KIND_R8
  use NUOPC,  only : NUOPC_CompAttributeGet

  use lnd_comp_kind  , only : r4 => shr_kind_r4
  use lnd_comp_kind  , only : r8 => shr_kind_r8 
  use lnd_comp_kind  , only : cl => shr_kind_cl
  use lnd_comp_types , only : noahmp_type
  use lnd_comp_types , only : field_type
  use lnd_comp_shr   , only : chkerr
  use lnd_comp_io    , only : read_tiled_file

  implicit none
  private

  public :: lnd_set_decomp_and_domain_from_mosaic

  !--------------------------------------------------------------------------
  ! Private module data
  !--------------------------------------------------------------------------

  integer :: dbug = 1
  character(*), parameter :: modName = "(lnd_comp_domain)"

  character(len=*), parameter :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine lnd_set_decomp_and_domain_from_mosaic(gcomp, noahmp, rc)

    ! input/output variables
    type(ESMF_GridComp), intent(in)  :: gcomp
    type(noahmp_type), intent(inout) :: noahmp
    integer, intent(out)             :: rc

    ! local variables
    real(r4), target, allocatable :: tmpr4(:)
    integer                          :: n
    integer                          :: decomptile(2,6)
    integer                          :: maxIndex(2)
    type(ESMF_Decomp_Flag)           :: decompflagPTile(2,6)

    type(field_type)                 :: flds(1)
    integer                          :: numOwnedElements, spatialDim, rank

    integer                          :: tlb(1), tub(1), tc(1)
    real(r8), allocatable            :: ownedElemCoords(:)
    integer, allocatable             :: vegtype(:)
    real(ESMF_KIND_R8), pointer      :: ptr1d(:)
    real(ESMF_KIND_R8), pointer      :: ptr2d(:,:)
    real(ESMF_KIND_R8), pointer      :: ptr3d(:,:,:)
    character(len=CL)                :: msg, filename
    logical                          :: isPresent, isSet
    type(ESMF_Field)                 :: field, farea
    type(ESMF_CoordSys_Flag)         :: coordSys
    real(r8), parameter              :: con_rerth = 6.3712e+6_r8
    character(len=*), parameter      :: subname = trim(modName)//':(lnd_set_decomp_and_domain_from_mosaic) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! ---------------------
    ! Set decomposition and decide it is regional or global
    ! ---------------------

    noahmp%domain%global = .true.

    ! ---------------------
    ! Create ESMF grid 
    ! ---------------------

    if (noahmp%domain%global) then
       ! set number of tiles
       noahmp%domain%ntiles = 6

       ! set decomposition
       do n = 1, noahmp%domain%ntiles
          decomptile(1,n) = noahmp%domain%layout(1)
          decomptile(2,n) = noahmp%domain%layout(2)
          decompflagPTile(:,n) = (/ ESMF_DECOMP_SYMMEDGEMAX, ESMF_DECOMP_SYMMEDGEMAX /)
       end do

       ! create grid
       noahmp%domain%grid = ESMF_GridCreateMosaic(filename=trim(noahmp%nmlist%mosaic_file), &
              regDecompPTile=decomptile, tileFilePath=trim(noahmp%nmlist%input_dir), &
              decompflagPTile=decompflagPTile, &
              staggerlocList=(/ESMF_STAGGERLOC_CENTER, ESMF_STAGGERLOC_CORNER/), &
              indexflag=ESMF_INDEX_GLOBAL, &
              name='lnd_grid', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! retrieve coordinate information, this will be used for model output
       call ESMF_GridGetCoord(noahmp%domain%grid, coordDim=1, farrayPtr=ptr2d, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       if (.not. allocated(noahmp%domain%lont)) then
          allocate(noahmp%domain%lont(lbound(ptr2d,dim=1):ubound(ptr2d,dim=1), lbound(ptr2d,dim=2):ubound(ptr2d,dim=2)))
       end if
       noahmp%domain%lont(:,:) = ptr2d(:,:)
       nullify(ptr2d)

       call ESMF_GridGetCoord(noahmp%domain%grid, coordDim=2, farrayPtr=ptr2d, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       if (.not. allocated(noahmp%domain%latt)) then
          allocate(noahmp%domain%latt(lbound(ptr2d,dim=1):ubound(ptr2d,dim=1),lbound(ptr2d,dim=2):ubound(ptr2d,dim=2)))
       end if
       noahmp%domain%latt(:,:) = ptr2d(:,:)
       nullify(ptr2d)

       ! query grid resolution (96, 384 etc.)
       call ESMF_GridGet(noahmp%domain%grid, 1, ESMF_STAGGERLOC_CENTER, maxIndex=maxIndex, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       noahmp%domain%ni = maxIndex(1)
       noahmp%domain%nj = maxIndex(2)
    else
       ! TODO: need to define grid for regional application such as HAFS
       call ESMF_LogWrite(trim(subname)//": "//' number of tile is 1, regional application is not supported!', ESMF_LOGMSG_ERROR)
       rc = ESMF_FAILURE
       return
    end if

    ! ---------------------
    ! Add mask to the grid
    ! ---------------------

    call ESMF_GridAddItem(noahmp%domain%grid, itemflag=ESMF_GRIDITEM_MASK, &
         staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! ---------------------
    ! Convert ESMF grid to mesh 
    ! ---------------------

    noahmp%domain%mesh = ESMF_MeshCreate(noahmp%domain%grid, 'lnd_mesh_from_grid', rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! ---------------------
    ! Query sizes from mesh
    ! ---------------------

    call ESMF_MeshGetFieldBounds(noahmp%domain%mesh, meshloc=ESMF_MESHLOC_ELEMENT, &
      totalLBound=tlb, totalUBound=tub, totalCount=tc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    noahmp%domain%begl = tlb(1)
    noahmp%domain%endl = tub(1)
    noahmp%static%im = tc(1)
    write(msg, fmt='(A,3I5)') trim(subname)//' : begl, endl, im = ', noahmp%domain%begl, &
      noahmp%domain%endl, noahmp%static%im
    call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO)

    !----------------------
    ! allocate temporary data structures
    !----------------------

    if (.not. allocated(tmpr4)) then
       allocate(tmpr4(noahmp%domain%begl:noahmp%domain%endl))
       tmpr4(:) = 0.0
    end if

    ! ---------------------
    ! Get fraction from orography file
    ! ---------------------

    ! read field
    filename = trim(noahmp%nmlist%input_dir)//'oro_data.tile*.nc'
    flds(1)%short_name = 'land_frac'
    flds(1)%ptr1r4 => tmpr4
    call read_tiled_file(noahmp, filename, flds, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! allocate data
    if (.not. allocated(noahmp%domain%frac)) then
       allocate(noahmp%domain%frac(noahmp%domain%begl:noahmp%domain%endl))
    end if
    noahmp%domain%frac = dble(tmpr4)

    ! ---------------------
    ! Read one of the static files to get mask information. This will be used to fix
    ! land-sea mask inconsistency among the files and it is documented in
    ! following link: https://github.com/ufs-community/ufs-weather-model/issues/1423 
    ! ---------------------

    ! read field
    write(filename, fmt="(A,I0,A)") trim(noahmp%nmlist%input_dir)//'C',noahmp%domain%ni, '.vegetation_type.tile*.nc'
    flds(1)%short_name = 'vegetation_type'
    flds(1)%ptr1r4 => tmpr4
    call read_tiled_file(noahmp, filename, flds, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! allocate data
    if (.not. allocated(vegtype)) then
       allocate(vegtype(noahmp%domain%begl:noahmp%domain%endl))
    end if
    vegtype(:) = int(tmpr4)

    ! ---------------------
    ! Calculate mask from land-sea fraction
    ! ---------------------

    if (.not. allocated(noahmp%domain%mask)) then
       allocate(noahmp%domain%mask(noahmp%domain%begl:noahmp%domain%endl))
    end if

    where (noahmp%domain%frac(:) > 0.0_r8 .and. vegtype(:) >= 0)
       noahmp%domain%mask(:) = 1
    elsewhere
       noahmp%domain%mask(:) = 0
    end where

    if (allocated(vegtype)) deallocate(vegtype)

    ! ---------------------
    ! Set mask in mesh 
    ! ---------------------

    call ESMF_MeshSet(noahmp%domain%mesh, elementMask=noahmp%domain%mask, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! ---------------------
    ! Get height from orography file
    ! ---------------------

    ! read field
    filename = trim(noahmp%nmlist%input_dir)//'oro_data.tile*.nc'
    flds(1)%short_name = 'orog_raw'
    flds(1)%ptr1r4 => tmpr4
    call read_tiled_file(noahmp, filename, flds, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! allocate data
    if (.not. allocated(noahmp%domain%hgt)) then
       allocate(noahmp%domain%hgt(noahmp%domain%begl:noahmp%domain%endl))
    end if
    noahmp%domain%hgt = dble(tmpr4)

    ! ---------------------
    ! Query cell area 
    ! ---------------------

    ! create field in R8 type
    farea = ESMF_FieldCreate(noahmp%domain%mesh, ESMF_TYPEKIND_R8, &
      meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)

    ! get cell area to the field
    call ESMF_FieldRegridGetArea(farea, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (.not. allocated(noahmp%domain%garea)) then
       allocate(noahmp%domain%garea(noahmp%domain%begl:noahmp%domain%endl))
    end if

    ! retrieve pointer and fill area array
    call ESMF_FieldGet(farea, farrayPtr=ptr1d, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    noahmp%domain%garea(:) = ptr1d(:)

    ! make unit conversion from square radians to square meters if it is required
    call ESMF_MeshGet(noahmp%domain%mesh, coordSys=coordSys, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (coordSys /= ESMF_COORDSYS_CART) then
       noahmp%domain%garea(:) = noahmp%domain%garea(:)*(con_rerth**2)
    end if

    ! ---------------------
    ! Query coordiates from ESMF mesh
    ! ---------------------

    ! determine dimensions in mesh
    call ESMF_MeshGet(noahmp%domain%mesh, spatialDim=spatialDim, &
      numOwnedElements=numOwnedElements, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! obtain mesh longitudes and latitudes
    if (.not. allocated(noahmp%domain%lats)) then
       allocate(noahmp%domain%lats(numOwnedElements))
    end if
    if (.not. allocated(noahmp%domain%lons)) then
       allocate(noahmp%domain%lons(numOwnedElements))
    end if
    if (.not. allocated(ownedElemCoords)) then
       allocate(ownedElemCoords(spatialDim*numOwnedElements))
    end if
    call ESMF_MeshGet(noahmp%domain%mesh, ownedElemCoords=ownedElemCoords)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    do n = 1,numOwnedElements
       noahmp%domain%lons(n) = ownedElemCoords(2*n-1)
       noahmp%domain%lats(n) = ownedElemCoords(2*n)
    end do
    deallocate(ownedElemCoords)

    ! ---------------------
    ! Clean memory
    ! ---------------------

    if (allocated(tmpr4)) deallocate(tmpr4)

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine lnd_set_decomp_and_domain_from_mosaic

end module lnd_comp_domain
