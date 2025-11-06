module NoahmpMosaicSortTileCatMod

!!! This module part of NoahMP Mosaic/Subgrid Tiling Scheme
!!! Purpose: To sort and identify most dominant lulc/soiltype/hydro types in a grid
!!!          identify dominant that contributed to > 90% grid area and 
!!!          scale the dominant ones to 100%

! ------------------------ Code history -----------------------------------
! Original code : Prasanth Valayamkunnath (IISER Thiruvananthapuram)
! Date          : July 10, 2025
! -------------------------------------------------------------------------

  use Machine
  use NoahmpIOVarType

  implicit none

contains

!=== sort landuse/soiltype/hydrotype index based on area fraction and identify most dominant types

  subroutine NoahmpMosaicSortTileCat (NoahmpIO)

    implicit none

    type(NoahmpIO_type), intent(inout)  :: NoahmpIO

!---------------------------------------------------------------------
!  Local Variables
!---------------------------------------------------------------------

    integer                                           :: i, j, k, m 
    integer                                           :: max_index 
    integer                                           :: temp_idx
    integer                                           :: n_dominant
    integer,                allocatable, dimension(:) :: sorted_indices
    real(kind=kind_noahmp)                            :: cumulative_sum
    real(kind=kind_noahmp)                            :: temp_val
    real(kind=kind_noahmp), allocatable, dimension(:) :: frac_vec
    real(kind=kind_noahmp), allocatable, dimension(:) :: sorted_values
    real(kind=kind_noahmp), allocatable, dimension(:) :: rescaled_values
! -------------------------------------------------------------------------
    associate(                                                      &
              XSTART               =>  NoahmpIO%XSTART             ,&
              XEND                 =>  NoahmpIO%XEND               ,&
              YSTART               =>  NoahmpIO%YSTART             ,&
              YEND                 =>  NoahmpIO%YEND               ,&
              NumMosaicCat         =>  noahmpio%NumMosaicCat       ,&
              SubGrdFrac           =>  NoahmpIO%SubGrdFrac         ,&
              SubGrdFracRescaled   =>  NoahmpIO%SubGrdFracRescaled ,& 
              SubGrdIndexSorted    =>  NoahmpIO%SubGrdIndexSorted  ,&
              NumberOfTiles        =>  NoahmpIO%NumberOfTiles      ,&
              NTilesMax            =>  NoahmpIO%NTilesMax          ,&
              NTiles_user          =>  NoahmpIO%IOPT_MOSAIC_NTILES  &
             )
! -------------------------------------------------------------------------

!   Step1: Allocate local vectors
    if ( .not. allocated (sorted_indices) )     allocate ( sorted_indices     (NumMosaicCat) )
    if ( .not. allocated (frac_vec)       )     allocate ( frac_vec           (NumMosaicCat) )
    if ( .not. allocated (sorted_values)  )     allocate ( sorted_values      (NumMosaicCat) )   
    if ( .not. allocated (rescaled_values))     allocate ( rescaled_values    (NumMosaicCat) )


 ! initialize temp variables   
    temp_val       = 0.
    NTilesMax      = 0
    cumulative_sum = 0.

!   Step 2: Loop over each grid cell (i,j)
    do i = XSTART, XEND
       do j = YSTART, YEND

          ! Extract fractions and initialize sorting arrays
          do k = 1, NumMosaicCat
             if (k .eq. NoahmpIO%ISWATER) then     ! to skip water fractions in the grid
               SubGrdFrac(i,k,j) = 0.
             endif
             frac_vec(k)         = SubGrdFrac(i,k,j)
             sorted_indices(k)   = k
             sorted_values(k)    = frac_vec(k)
          end do

          ! Sort descending by fraction value (Selection Sort)
          do k = 1, NumMosaicCat - 1
            max_index = k
            do m = k + 1, NumMosaicCat
              if (sorted_values(m) > sorted_values(max_index)) then
                 max_index = m
              end if
            end do

            ! Swap values
            temp_val = sorted_values(k)
            sorted_values(k) = sorted_values(max_index)
            sorted_values(max_index) = temp_val

            temp_idx = sorted_indices(k)
            sorted_indices(k) = sorted_indices(max_index)
            sorted_indices(max_index) = temp_idx
          end do

          ! Step 3: Identify dominant types (sum > 0.9)
          cumulative_sum = 0.0
          n_dominant = 0
          do k = 1, NumMosaicCat
             if (sorted_values(k) == 0) exit
             ! Now limit the number of tiles to user defined number from the namelist
             if (k <= NTiles_user) then 
                n_dominant = k
                cumulative_sum = cumulative_sum + sorted_values(k)  !
                if (cumulative_sum > 0.9) exit
             else 
                exit
             endif
          end do

          ! Step 4: Rescale selected values to sum to 1
          do k = 1, n_dominant
             rescaled_values(k) = sorted_values(k) / cumulative_sum
          end do

          !print for diagnosis
          if(i .eq. 150 .and. j .eq. 150) then
              WRITE(*,'(A,I2,A,I2,A)') 'Grid (', i, ',', j, '): Dominant Land Cover Types (>90% of total)'
              DO k = 1, n_dominant
                WRITE(*,'(A,I2,A,F7.4,A,F7.4)') '  Type: ', sorted_indices(k), &
                        '  Original: ', sorted_values(k), '  Rescaled: ', rescaled_values(k)
              END DO
              WRITE(*,'(A,F7.4)') '  Sum of Rescaled Fractions: ', SUM(rescaled_values(1:n_dominant))
              PRINT *, '-------------------------------------------------------------'
          end if

          do k = 1, n_dominant
             SubGrdFracRescaled(i,j,k)  =  rescaled_values(k)
             SubGrdIndexSorted(i,j,k)   =  sorted_indices(k)
             NumberOfTiles(i,j)         =  n_dominant ! it is <= NTiles_user
             NTilesMax                  =  min(max(NTilesMax,n_dominant), NTiles_user) ! maximum value across domain
          enddo
       end do
    end do

    end associate

  end subroutine NoahmpMosaicSortTileCat

end module NoahmpMosaicSortTileCatMod
