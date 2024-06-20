!WRF:DRIVER_LAYER:DOMAIN_OBJECT

MODULE module_domain

    ! extract from module_domain_type.F from WRF/frame/
    ! PLACEHOLDER for domain in GROUNDWATER_INIT subroutine which is not used but defined.
       TYPE domain

          INTEGER                                             :: id
          INTEGER                                             :: domdesc
          INTEGER                                             :: communicator
          INTEGER                                             :: iocommunicator
          INTEGER,POINTER                                     :: mapping(:,:)
          INTEGER,POINTER                                     :: i_start(:),i_end(:)
          INTEGER,POINTER                                     :: j_start(:),j_end(:)
          INTEGER                                             :: max_tiles
          INTEGER                                             :: num_tiles        ! taken out of namelist 20000908
          INTEGER                                             :: num_tiles_x      ! taken out of namelist 20000908
          INTEGER                                             :: num_tiles_y      ! taken out of namelist 20000908
          INTEGER                                             :: num_tiles_spec   ! place to store number of tiles computed from
                                                                                  ! externally specified params
       END TYPE domain

       TYPE(domain) , POINTER :: head_grid , new_grid , next_grid , old_grid

    CONTAINS

    END MODULE module_domain