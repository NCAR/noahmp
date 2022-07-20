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
  use ESMF ,  only : ESMF_COORDSYS_CART, ESMF_KIND_R8
  use NUOPC,  only : NUOPC_CompAttributeGet

  use lnd_comp_kind  , only : r8 => shr_kind_r8 
  use lnd_comp_kind  , only : cl => shr_kind_cl
  use lnd_comp_types , only : noahmp_type
  use lnd_comp_shr   , only : chkerr
  use lnd_comp_io    , only : read_tiled_file

  use fms2_io_mod    , only : FmsNetcdfFile_t, open_file
  use mosaic2_mod    , only : get_mosaic_ntiles, get_mosaic_grid_sizes
  use mosaic2_mod    , only : get_mosaic_contact, get_mosaic_ncontacts
  use mpp_mod        , only : mpp_root_pe
  use mpp_domains_mod, only : mpp_define_mosaic, mpp_domains_init
  use mpp_domains_mod, only : mpp_define_layout

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
    integer                          :: n, numOwnedElements, spatialDim, rank
    integer                          :: decomptile(2,6)
    real(r8), allocatable            :: ownedElemCoords(:)
    real(ESMF_KIND_R8), pointer      :: ptr1d(:)
    real(ESMF_KIND_R8), pointer      :: ptr2d(:,:)
    real(ESMF_KIND_R8), pointer      :: ptr3d(:,:,:)
    character(len=CL)                :: msg, filename
    logical                          :: isPresent, isSet
    type(FmsNetcdfFile_t)            :: mosaic_fileobj
    type(ESMF_Field)                 :: field, farea
    type(ESMF_Decomp_Flag)           :: decompflagPTile(2,6)
    type(ESMF_CoordSys_Flag)         :: coordSys
    real(r8), parameter              :: con_rerth = 6.3712e+6_r8
    character(len=*), parameter      :: subname = trim(modName)//':(lnd_set_decomp_and_domain_from_mosaic) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! ---------------------
    ! Open mosaic file and query some information 
    ! ---------------------

    if (.not. open_file(mosaic_fileobj, trim(noahmp%nmlist%mosaic_file), 'read')) then
       call ESMF_LogWrite(trim(subname)//'error in opening file '//trim(noahmp%nmlist%mosaic_file), ESMF_LOGMSG_ERROR)
       rc = ESMF_FAILURE
       return
    end if

    ! query number of tiles
    noahmp%domain%ntiles = get_mosaic_ntiles(mosaic_fileobj)
    noahmp%domain%global = (noahmp%domain%ntiles > 1)

    ! query domain sizes for each tile
    if (.not. allocated(noahmp%domain%nit)) allocate(noahmp%domain%nit(noahmp%domain%ntiles))
    if (.not. allocated(noahmp%domain%njt)) allocate(noahmp%domain%njt(noahmp%domain%ntiles))
    call get_mosaic_grid_sizes(mosaic_fileobj, noahmp%domain%nit, noahmp%domain%njt)

    ! query number of contacts
    noahmp%domain%ncontacts = get_mosaic_ncontacts(mosaic_fileobj)

    ! allocate required arrays to create FMS domain from mosaic file
    if (.not. allocated(noahmp%domain%tile1)) allocate(noahmp%domain%tile1(noahmp%domain%ncontacts))
    if (.not. allocated(noahmp%domain%tile2)) allocate(noahmp%domain%tile2(noahmp%domain%ncontacts))
    if (.not. allocated(noahmp%domain%istart1)) allocate(noahmp%domain%istart1(noahmp%domain%ncontacts))
    if (.not. allocated(noahmp%domain%iend1)) allocate(noahmp%domain%iend1(noahmp%domain%ncontacts))
    if (.not. allocated(noahmp%domain%jstart1)) allocate(noahmp%domain%jstart1(noahmp%domain%ncontacts))
    if (.not. allocated(noahmp%domain%jend1)) allocate(noahmp%domain%jend1(noahmp%domain%ncontacts))
    if (.not. allocated(noahmp%domain%istart2)) allocate(noahmp%domain%istart2(noahmp%domain%ncontacts))
    if (.not. allocated(noahmp%domain%iend2)) allocate(noahmp%domain%iend2(noahmp%domain%ncontacts))
    if (.not. allocated(noahmp%domain%jstart2)) allocate(noahmp%domain%jstart2(noahmp%domain%ncontacts))
    if (.not. allocated(noahmp%domain%jend2)) allocate(noahmp%domain%jend2(noahmp%domain%ncontacts))

    ! query domain related information
    call get_mosaic_contact(mosaic_fileobj, noahmp%domain%tile1, noahmp%domain%tile2, &
         noahmp%domain%istart1, noahmp%domain%iend1, noahmp%domain%jstart1, noahmp%domain%jend1, &
         noahmp%domain%istart2, noahmp%domain%iend2, noahmp%domain%jstart2, noahmp%domain%jend2)

    do n = 1, noahmp%domain%ncontacts
       write(msg, fmt='(A,I2,A,2I5)') trim(subname)//' : tile1, tile2 (', n ,') = ', &
         noahmp%domain%tile1(n), noahmp%domain%tile2(n)
       call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO)    
       write(msg, fmt='(A,I2,A,4I5)') trim(subname)//' : istart1, iend1, jstart1, jend1 (', n ,') = ', &
         noahmp%domain%istart1(n), noahmp%domain%iend1(n), noahmp%domain%jstart1(n), noahmp%domain%jend1(n)
       call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO)    
       write(msg, fmt='(A,I2,A,4I5)') trim(subname)//' : istart2, iend2, jstart2, jend2 (', n ,') = ', &
         noahmp%domain%istart2(n), noahmp%domain%iend2(n), noahmp%domain%jstart2(n), noahmp%domain%jend2(n)
       call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO)    
    end do

    ! ---------------------
    ! Create FMS domain 
    ! ---------------------

    call lnd_domain_create(gcomp, noahmp, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! ---------------------
    ! Create ESMF grid 
    ! ---------------------

    if (noahmp%domain%global) then
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
    ! Get fraction from orography file 
    ! ---------------------

    ! input file name, the tile will be added to it based on the active PET
    filename = trim(noahmp%nmlist%input_dir)//'oro_data.tile'

    ! read data to ESMF field
    call read_tiled_file(filename, 'land_frac', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! get pointer
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr3d, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%domain%begl = lbound(ptr3d, dim=1)
    noahmp%domain%endl = ubound(ptr3d, dim=1)
    noahmp%static%im = noahmp%domain%endl-noahmp%domain%begl+1

    ! allocate variable and fill it
    if (.not. allocated(noahmp%domain%frac)) then
       allocate(noahmp%domain%frac(noahmp%domain%begl:noahmp%domain%endl))
    end if
    noahmp%domain%frac(:) = ptr3d(:,1,1)
    nullify(ptr3d)

    ! ---------------------
    ! Calculate mask from land-sea fraction
    ! ---------------------

    if (.not. allocated(noahmp%domain%mask)) then
       allocate(noahmp%domain%mask(noahmp%domain%begl:noahmp%domain%endl))
    end if

    where (noahmp%domain%frac(:) > 0.0_r8)
       noahmp%domain%mask(:)  = 1.0_r8
    elsewhere
       noahmp%domain%mask(:)  = 0.0_r8
    end where

    ! ---------------------
    ! Set mask in mesh 
    ! ---------------------

    call ESMF_MeshSet(noahmp%domain%mesh, elementMask=noahmp%domain%mask, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! ---------------------
    ! Get height from orography file
    ! ---------------------

    ! read data to ESMF field
    call read_tiled_file(filename, 'orog_raw', noahmp, field, numrec=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! get pointer
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr3d, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! allocate variable
    if (.not. allocated(noahmp%domain%hgt)) then
       allocate(noahmp%domain%hgt(noahmp%domain%begl:noahmp%domain%endl))
    end if
    noahmp%domain%hgt(:) = ptr3d(:,1,1)

    ! clean memory
    nullify(ptr3d)

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

    call ESMF_FieldDestroy(field, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldDestroy(farea, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine lnd_set_decomp_and_domain_from_mosaic

  !===============================================================================
  subroutine lnd_domain_create(gcomp, noahmp, rc)

    ! input/output variables
    type(ESMF_GridComp), intent(in)  :: gcomp
    type(noahmp_type), intent(inout) :: noahmp
    integer, intent(inout)           :: rc

    ! local variables
    type(ESMF_VM)               :: vm
    integer                     :: n, npet, npes_per_tile
    integer                     :: halo = 0
    integer                     :: global_indices(4,6)
    integer                     :: layout2d(2,6)
    integer, allocatable        :: pe_start(:), pe_end(:)
    character(len=cl)           :: msg
    character(len=*), parameter :: subname=trim(modName)//':(lnd_domain_create) '
    !-------------------------------------------------------------------------------

    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! ---------------------
    ! Query components 
    ! ---------------------

    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_VMGet(vm=vm, petCount=npet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------
    ! Initialize domain 
    !----------------------

    call mpp_domains_init()

    !----------------------
    ! Create domain 
    !----------------------

    ! setup global indices
    do n = 1, noahmp%domain%ntiles
       global_indices(1,n) = 1
       global_indices(2,n) = noahmp%domain%nit(n)
       global_indices(3,n) = 1
       global_indices(4,n) = noahmp%domain%njt(n)
    enddo

    ! check total number of PETs
    if (mod(npet, noahmp%domain%ntiles) /= 0) then
       write(msg, fmt='(A,I5)') trim(subname)//' : nPet should be multiple of 6 to read initial conditions but it is ', npet 
       call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_ERROR)
       rc = ESMF_FAILURE
       return
    end if

    ! calculate layout if it is not provided as configuration option
    if (noahmp%domain%layout(1) < 0 .and. noahmp%domain%layout(2) < 0) then
       npes_per_tile = npet/noahmp%domain%ntiles
       call mpp_define_layout(global_indices(:,1), npes_per_tile, noahmp%domain%layout)
    end if

    ! set layout and print out debug information
    do n = 1, noahmp%domain%ntiles
       layout2d(:,n) = noahmp%domain%layout(:)
       write(msg, fmt='(A,I2,A,2I5)') trim(subname)//' layout (', n ,') = ', layout2d(1,n), layout2d(2,n)
       call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO)
       write(msg, fmt='(A,I2,A,4I5)') trim(subname)//' global_indices (', n,') = ', &
         global_indices(1,n), global_indices(2,n), global_indices(3,n), global_indices(4,n)
       call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO)
    enddo

    !----------------------
    ! Set pe_start, pe_end 
    !----------------------

    allocate(pe_start(noahmp%domain%ntiles))
    allocate(pe_end(noahmp%domain%ntiles))
    do n = 1, noahmp%domain%ntiles
       pe_start(n) = mpp_root_pe()+(n-1)*noahmp%domain%layout(1)*noahmp%domain%layout(2)
       pe_end(n) = mpp_root_pe()+n*noahmp%domain%layout(1)*noahmp%domain%layout(2)-1
       write(msg, fmt='(A,I2,A,2I5)') trim(subname)//' pe_start, pe_end (', n ,') = ', pe_start(n), pe_end(n)
       call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO)
    enddo

    call mpp_define_mosaic(global_indices, layout2d, noahmp%domain%mosaic_domain, &
         noahmp%domain%ntiles, noahmp%domain%ncontacts, noahmp%domain%tile1, noahmp%domain%tile2, &
         noahmp%domain%istart1, noahmp%domain%iend1, noahmp%domain%jstart1, noahmp%domain%jend1, &
         noahmp%domain%istart2, noahmp%domain%iend2, noahmp%domain%jstart2, noahmp%domain%jend2, &
         pe_start, pe_end, symmetry=.true., whalo=halo, ehalo=halo, shalo=halo, nhalo=halo, &
         name='lnd domain')

    !----------------------
    ! Deallocate temporary arrays
    !----------------------

    deallocate(pe_start)
    deallocate(pe_end)

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine lnd_domain_create

end module lnd_comp_domain
