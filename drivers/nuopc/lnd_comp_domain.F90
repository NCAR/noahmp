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
  use ESMF ,  only : ESMF_MeshDestroy, ESMF_DistGridCreate, ESMF_VMAllGatherV
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
    real(r4), target, allocatable    :: tmpr4(:)
    real(r8), target, allocatable    :: tmpr8(:)
    integer                          :: n, petCount, localPet
    integer                          :: decomptile(2,6)
    integer                          :: maxIndex(2)
    type(ESMF_Decomp_Flag)           :: decompflagPTile(2,6)
    type(ESMF_VM)                    :: vm

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
    ! Query VM
    ! ---------------------

    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_VMGet(vm, petCount=petCount, localPet=localPet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

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

       ! check user provided layout
       if (petCount /= noahmp%domain%layout(1)*noahmp%domain%layout(2)*noahmp%domain%ntiles) then
          call ESMF_LogWrite(trim(subname)//": ERROR in layout. layout_x * layout_y * 6 != #PETs", ESMF_LOGMSG_INFO)
          rc = ESMF_FAILURE
          return
       end if

       ! use user provided layout to set decomposition
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
    if (.not. allocated(tmpr8)) then
       allocate(tmpr8(noahmp%domain%begl:noahmp%domain%endl))
       tmpr8(:) = 0.0
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

    ! allocate data
    if (.not. allocated(vegtype)) then
       allocate(vegtype(noahmp%domain%begl:noahmp%domain%endl))
    end if

    ! read field
    if (trim(noahmp%nmlist%ic_type) == 'sfc') then
       filename = trim(noahmp%nmlist%input_dir)//'sfc_data.tile*.nc'
       flds(1)%short_name = 'vtype'
       flds(1)%ptr1r8 => tmpr8
       call read_tiled_file(noahmp, filename, flds, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       vegtype(:) = int(tmpr8)
    else
       write(filename, fmt="(A,I0,A)") trim(noahmp%nmlist%input_dir)//'C',noahmp%domain%ni, '.vegetation_type.tile*.nc'
       flds(1)%short_name = 'vegetation_type'
       flds(1)%ptr1r4 => tmpr4
       call read_tiled_file(noahmp, filename, flds, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       vegtype(:) = int(tmpr4)
    end if

    ! ---------------------
    ! Calculate mask from land-sea fraction
    ! ---------------------

    if (.not. allocated(noahmp%domain%mask)) then
       allocate(noahmp%domain%mask(noahmp%domain%begl:noahmp%domain%endl))
    end if

    where (noahmp%domain%frac(:) > 0.0_r8 .and. vegtype(:) > 0)
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
    ! Modify decomposition to evenly distribute land and ocean points
    ! ---------------------

    if (trim(noahmp%nmlist%decomp_type) == 'custom') then
       ! modify decomposition
       call lnd_modify_decomp(gcomp, noahmp, rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if

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
    if (allocated(tmpr8)) deallocate(tmpr8)

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine lnd_set_decomp_and_domain_from_mosaic

  !===============================================================================
  subroutine lnd_modify_decomp(gcomp, noahmp, rc)

    ! input/output variables
    type(ESMF_GridComp), intent(in)  :: gcomp    
    type(noahmp_type), intent(inout) :: noahmp
    integer, intent(out)             :: rc

    ! local variables
    type(ESMF_VM)                    :: vm
    type(ESMF_Mesh)                  :: mesh
    type(ESMF_DistGrid)              :: distgrid, distgrid_new
    type(field_type)                 :: flds(1)    
    integer                          :: n, m, g
    integer                          :: petCount, localPet
    integer                          :: lsize, gsize(1)
    integer                          :: nlnd, nocn
    integer                          :: begl_l, endl_l, begl_o, endl_o
    integer, allocatable             :: nlnd_loc(:)
    integer, allocatable             :: nocn_loc(:)
    integer, allocatable             :: mask_glb(:)
    integer, allocatable             :: gindex_loc(:)
    integer, allocatable             :: gindex_glb(:)
    integer, allocatable             :: gindex_lnd(:)
    integer, allocatable             :: gindex_ocn(:)
    integer, allocatable             :: gindex_new(:)
    integer, allocatable             :: lsize_arr(:)
    real(r4), target, allocatable    :: tmpr4(:)
    character(len=cl)                :: msg, filename
    character(len=*), parameter      :: subname = trim(modName)//':(lnd_modify_decomp) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! ---------------------
    ! Query VM
    ! ---------------------

    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_VMGet(vm, petCount=petCount, localPet=localPet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! ---------------------
    ! Query existing mesh
    ! ---------------------

    ! retrive default distgrid
    call ESMF_MeshGet(noahmp%domain%mesh, elementdistGrid=distgrid, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! get local number of elements
    call ESMF_DistGridGet(distgrid, localDe=0, elementCount=lsize, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! calculate number of elements globally
    call ESMF_VMAllReduce(vm, (/ lsize /), gsize, 1, ESMF_REDUCE_SUM, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    nlnd = count(noahmp%domain%mask(:) > 0, dim=1)
    nocn = lsize-nlnd
    write(msg, fmt='(A,4I8)') trim(subname)//' : lsize, gsize, nlnd, nocn = ', &
      lsize, gsize(1), nlnd, nocn
    call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO)

    ! get default sequence index
    allocate(gindex_loc(lsize))
    call ESMF_DistGridGet(distgrid, 0, seqIndexList=gindex_loc, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! ---------------------
    ! Create global view of indexes and mask
    ! --------------------

    ! collect local sizes
    allocate(lsize_arr(petCount))
    call ESMF_VMAllGatherV(vm, sendData=(/ lsize /), sendCount=1, &
      recvData=lsize_arr, recvCounts=(/ (1, n = 0, petCount-1) /), &
      recvOffsets=(/ (n, n = 0, petCount-1) /), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! global view of indexes
    allocate(gindex_glb(gsize(1)))
    gindex_glb(:) = 0
    call ESMF_VMAllGatherV(vm, sendData=gindex_loc, sendCount=lsize, &
      recvData=gindex_glb, recvCounts=lsize_arr, &
      recvOffsets=(/ (sum(lsize_arr(1:n))-lsize_arr(1), n = 1, petCount) /), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! global view of land sea mask
    allocate(mask_glb(gsize(1)))
    mask_glb(:) = 0
    call ESMF_VMAllGatherV(vm, sendData=noahmp%domain%mask, sendCount=lsize, &
      recvData=mask_glb, recvCounts=lsize_arr, &
      recvOffsets=(/ (sum(lsize_arr(1:n))-lsize_arr(1), n = 1, petCount) /), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! ---------------------    
    ! split global indexes as land and ocean
    ! ---------------------

    nlnd = count(mask_glb > 0, dim=1)
    nocn = gsize(1)-nlnd
    allocate(gindex_lnd(nlnd))
    gindex_lnd = 0
    allocate(gindex_ocn(nocn))
    gindex_ocn = 0

    n = 0
    m = 0
    do g = 1, gsize(1)
       if (mask_glb(g) > 0) then
          n = n+1
          gindex_lnd(n) = gindex_glb(g)
       else
          m = m+1
          gindex_ocn(m) = gindex_glb(g)
       end if
    end do

    ! ---------------------    
    ! create new local indexes
    ! ---------------------

    allocate(nlnd_loc(0:petCount-1))
    nlnd_loc = 0
    allocate(nocn_loc(0:petCount-1))
    nocn_loc = 0
    do n = 0, petCount-1
       nlnd_loc(n) = nlnd/petCount
       nocn_loc(n) = nocn/petCount
       if (n < mod(nlnd, petCount)) then
          nlnd_loc(n) = nlnd_loc(n)+1
       else
          nocn_loc(n) = nocn_loc(n)+1
       end if
    end do
    if (localPet == 0) then
       begl_l = 1
       begl_o = 1
    else
       begl_l = sum(nlnd_loc(0:localPet-1))+1
       begl_o = sum(nocn_loc(0:localPet-1))+1
    end if
    endl_l = sum(nlnd_loc(0:localPet))
    endl_o = sum(nocn_loc(0:localPet))

    write(msg,'(A,6I8)') trim(subname)//' : nlnd_loc, nocn_loc, begl_l, endl_l, begl_o, endl_o = ', &
      nlnd_loc(localPet), nocn_loc(localPet), begl_l, endl_l, begl_o, endl_o 
    call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO)

    allocate(gindex_new(nlnd_loc(localPet)+nocn_loc(localPet)))
    gindex_new(:nlnd_loc(localPet)) = gindex_lnd(begl_l:begl_l)
    gindex_new(nlnd_loc(localPet)+1:) = gindex_ocn(begl_o:endl_o)

    if (noahmp%nmlist%debug_level > 10) then
       do n = 1, nlnd_loc(localPet)+nocn_loc(localPet)
          write(msg,'(A,2I8)') trim(subname)//' : n, gindex_new = ', n, gindex_new(n)
          call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO)
       end do
    end if

    ! ---------------------
    ! Update mesh with new decomposition
    ! ---------------------

    ! create new distgrid with new index
    distgrid_new = ESMF_DistGridCreate(arbSeqIndexList=gindex_new, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! create new mesh with new distgrid
    mesh = ESMF_MeshCreate(noahmp%domain%mesh, elementDistGrid=distgrid_new, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! destroy old and replace with new one
    call ESMF_MeshDestroy(noahmp%domain%mesh, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    noahmp%domain%mesh = mesh

    ! ---------------------
    ! fix mask and fraction to be consistent with new decomposition
    ! ---------------------

    ! mask
    noahmp%domain%mask(:nlnd_loc(localPet)) = 1
    noahmp%domain%mask(nlnd_loc(localPet)+1:) = 0

    ! fraction, read from file again
    allocate(tmpr4(noahmp%domain%begl:noahmp%domain%endl))
    tmpr4(:) = 0.0
    filename = trim(noahmp%nmlist%input_dir)//'oro_data.tile*.nc'
    flds(1)%short_name = 'land_frac'
    flds(1)%ptr1r4 => tmpr4
    call read_tiled_file(noahmp, filename, flds, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    deallocate(tmpr4)

    ! ---------------------
    ! Clean memory
    ! ---------------------    

    deallocate(nlnd_loc)
    deallocate(nocn_loc)
    deallocate(mask_glb)
    deallocate(gindex_loc)
    deallocate(gindex_glb)
    deallocate(gindex_lnd)
    deallocate(gindex_ocn)
    deallocate(gindex_new)
    deallocate(lsize_arr)

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine lnd_modify_decomp

end module lnd_comp_domain
