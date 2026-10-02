module NoahmpInitMainMod

!!!  Module to initialize Noah-MP 2-D variables

  use Machine
  use NoahmpIOVarType, only : NoahmpIO_type
  use NoahmpSnowInitMod
 
  implicit none
  
contains

  subroutine NoahmpInitMain(NoahmpIO)

! ------------------------ Code history -------------------------------------
! Original Noah-MP subroutine: NOAHMP_INIT
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! Sep 13, 2026: NoahmpIO%xx change to 1-D vector for MPAS, Cenlin He (NCAR)
! ---------------------------------------------------------------------------

    implicit none 
   
    type(NoahmpIO_type), intent(inout) :: NoahmpIO

    ! local variables
    integer                                     :: its,ite
    integer                                     :: I,errflag,NS,IZ
    logical                                     :: urbanpt_flag
    real(kind=kind_noahmp)                      :: BEXP, SMCMAX, PSISAT, FK
    real(kind=kind_noahmp), parameter           :: BLIM  = 5.5
    real(kind=kind_noahmp), parameter           :: HLICE = 3.335E5
    real(kind=kind_noahmp), parameter           :: GRAV0 = 9.81
    real(kind=kind_noahmp), parameter           :: T0    = 273.15
! --------------------------------------------------------------------------- 

    ! initialize
    its = NoahmpIO%its
    ite = NoahmpIO%ite

    ! only initialize for non-restart case
    if ( .not. NoahmpIO%restart_flag ) then

       ! initialize physical snow height SNOWH
       if ( .not. NoahmpIO%FNDSNOWH ) then
          ! If no SNOWH do the following
          print*, 'SNOW HEIGHT NOT FOUND - VALUE DEFINED IN LSMINIT'
          do I = its, ite
             NoahmpIO%SNOWH(I) = NoahmpIO%SNOW(I) * 0.005  ! SNOW in mm and SNOWH in m
          enddo
       endif
   
       ! Check if snow/snowh are consistent and cap SWE at 2000mm
       ! the Noah-MP code does it internally but if we don't do it here, problems ensue
       do I = its, ite
          if ( NoahmpIO%SNOW(I)  < 0.0 ) NoahmpIO%SNOW(I)  = 0.0 
          if ( NoahmpIO%SNOWH(I) < 0.0 ) NoahmpIO%SNOWH(I) = 0.0
          if ( (NoahmpIO%SNOW(I) > 0.0) .and. (NoahmpIO%SNOWH(I) == 0.0) ) &
             NoahmpIO%SNOWH(I) = NoahmpIO%SNOW(I) * 0.005
          if ( (NoahmpIO%SNOWH(I) > 0.0) .and. (NoahmpIO%SNOW(I) == 0.0) ) &
             NoahmpIO%SNOW(I)  = NoahmpIO%SNOWH(I) / 0.005
          if ( NoahmpIO%SNOW(I) > 2000.0 ) then
             NoahmpIO%SNOWH(I) = NoahmpIO%SNOWH(I) * 2000.0 / NoahmpIO%SNOW(I)      ! SNOW in mm and SNOWH in m
             NoahmpIO%SNOW (I) = 2000.0                                             ! cap SNOW at 2000, maintain density
          endif
       enddo

       ! Given the soil layer thicknesses (in DZS), initialize the soil layer
       ! depths from the surface.
       NoahmpIO%ZSOIL(1) = -NoahmpIO%DZS(1)          ! negative
       do NS = 2, NoahmpIO%NSOIL
          NoahmpIO%ZSOIL(NS) = NoahmpIO%ZSOIL(NS-1) - NoahmpIO%DZS(NS)
       enddo

       ! check soil type
       errflag = 0
       do I = its, ite
          if ( NoahmpIO%ISLTYP(I) < 1 ) then
             errflag = 1
             write(*,*) "lsminit: out of range ISLTYP ",I, NoahmpIO%ISLTYP(I)
             stop
          endif
       enddo

       ! initialize soil liquid water content SH2O
       do I = its, ite
          if ( (NoahmpIO%IVGTYP(I) == NoahmpIO%ISICE_TABLE) .and. &
               (NoahmpIO%XICE(I) <= 0.0) ) then
             do NS = 1, NoahmpIO%NSOIL
                NoahmpIO%SMOIS(I,NS) = 1.0  ! glacier starts all frozen
                NoahmpIO%SH2O(I,NS)  = 0.0
                NoahmpIO%TSLB(I,NS)  = min(NoahmpIO%TSLB(I,NS), 263.15) ! set glacier temp to at most -10C
             enddo
             ! NoahmpIO%TMN(I) = min(NoahmpIO%TMN(I), 263.15)           ! set deep temp to at most -10C
             NoahmpIO%SNOW(I)  = max(NoahmpIO%SNOW(I), 10.0)            ! set SWE to at least 10mm
             NoahmpIO%SNOWH(I) = NoahmpIO%SNOW(I) * 0.005               ! SNOW in mm and SNOWH in m
          else
             BEXP   = NoahmpIO%BEXP_TABLE  (NoahmpIO%ISLTYP(I))
             SMCMAX = NoahmpIO%SMCMAX_TABLE(NoahmpIO%ISLTYP(I))
             PSISAT = NoahmpIO%PSISAT_TABLE(NoahmpIO%ISLTYP(I))
             do NS = 1, NoahmpIO%NSOIL
                if ( NoahmpIO%SMOIS(I,NS) > SMCMAX ) NoahmpIO%SMOIS(I,NS) = SMCMAX
             enddo
             if ( (BEXP > 0.0) .and. (SMCMAX > 0.0) .and. (PSISAT > 0.0) ) then
                do NS = 1, NoahmpIO%NSOIL
                   if ( NoahmpIO%TSLB(I,NS) < T0 ) then
                      FK = (((HLICE / (GRAV0*(-PSISAT))) * &
                            ((NoahmpIO%TSLB(I,NS)-T0) / NoahmpIO%TSLB(I,NS)))**(-1/BEXP))*SMCMAX
                      FK = max(FK, 0.02)
                      NoahmpIO%SH2O(I,NS) = min(FK, NoahmpIO%SMOIS(I,NS))
                   else
                      NoahmpIO%SH2O(I,NS) = NoahmpIO%SMOIS(I,NS)
                   endif
                enddo
             else
                do NS = 1, NoahmpIO%NSOIL
                   NoahmpIO%SH2O(I,NS) = NoahmpIO%SMOIS(I,NS)
                enddo
             endif
          endif
       enddo ! I

       ! initilize other quantities
       do I = its, ite
          NoahmpIO%QTDRAIN(I)  = 0.0
          NoahmpIO%TVXY(I)     = NoahmpIO%TSK(I)
          NoahmpIO%TGXY(I)     = NoahmpIO%TSK(I)
          if ( (NoahmpIO%SNOW(I) > 0.0) .and. (NoahmpIO%TSK(I) > T0) ) NoahmpIO%TVXY(I) = T0
          if ( (NoahmpIO%SNOW(I) > 0.0) .and. (NoahmpIO%TSK(I) > T0) ) NoahmpIO%TGXY(I) = T0
          NoahmpIO%CANWAT(I)   = 0.0
          NoahmpIO%CANLIQXY(I) = NoahmpIO%CANWAT(I)
          NoahmpIO%CANICEXY(I) = 0.0
          NoahmpIO%EAHXY(I)    = 2000.0
          NoahmpIO%TAHXY(I)    = NoahmpIO%TSK(I)
          NoahmpIO%T2MVXY(I)   = NoahmpIO%TSK(I)
          NoahmpIO%T2MBXY(I)   = NoahmpIO%TSK(I)
          NoahmpIO%T2MXY(I)    = NoahmpIO%TSK(I)
          if ( (NoahmpIO%SNOW(I) > 0.0) .and. (NoahmpIO%TSK(I) > T0) ) NoahmpIO%TAHXY(I)  = T0
          if ( (NoahmpIO%SNOW(I) > 0.0) .and. (NoahmpIO%TSK(I) > T0) ) NoahmpIO%T2MVXY(I) = T0
          if ( (NoahmpIO%SNOW(I) > 0.0) .and. (NoahmpIO%TSK(I) > T0) ) NoahmpIO%T2MBXY(I) = T0
          if ( (NoahmpIO%SNOW(I) > 0.0) .and. (NoahmpIO%TSK(I) > T0) ) NoahmpIO%T2MXY(I)  = T0
          NoahmpIO%CMXY(I)     = 0.0
          NoahmpIO%CHXY(I)     = 0.0
          NoahmpIO%FWETXY(I)   = 0.0
          NoahmpIO%SNEQVOXY(I) = 0.0
          NoahmpIO%ALBOLDXY(I) = 0.65
          NoahmpIO%QSNOWXY(I)  = 0.0
          NoahmpIO%QRAINXY(I)  = 0.0
          NoahmpIO%WSLAKEXY(I) = 0.0
          if ( NoahmpIO%IOPT_WETLAND > 0 ) then
             NoahmpIO%FSATXY(I)   = 0.0
             NoahmpIO%WSURFXY(I)  = 0.0
          endif
          if ( NoahmpIO%IOPT_RUNSUB /= 5 ) then 
             NoahmpIO%WAXY(I)   = 4900.0 
             NoahmpIO%WTXY(I)   = NoahmpIO%WAXY(I) 
             NoahmpIO%ZWTXY(I)  = (25.0 - NoahmpIO%ZSOIL(NoahmpIO%NSOIL)) - NoahmpIO%WAXY(I)/1000/0.2
          else
             NoahmpIO%WAXY(I)   = 0.0
             NoahmpIO%WTXY(I)   = 0.0
             ! for MPAS, it is defined in mpas_atmphys_lsm_noahmpinit.F
             !NoahmpIO%AREAXY(I) = (max(10.0,NoahmpIO%DX(I)) * max(10.0,NoahmpIO%DY(I))) / &
             !                       (NoahmpIO%MSFTX(I) * NoahmpIO%MSFTY(I))
          endif

          urbanpt_flag = .false.
          if ( (NoahmpIO%IVGTYP(I) == NoahmpIO%ISURBAN_TABLE) .or. &
               (NoahmpIO%IVGTYP(I) > NoahmpIO%URBTYPE_beg) ) then
             urbanpt_flag = .true.
          endif

          if ( (NoahmpIO%IVGTYP(I) == NoahmpIO%ISBARREN_TABLE) .or. &
               (NoahmpIO%IVGTYP(I) == NoahmpIO%ISICE_TABLE) .or. &
               ((NoahmpIO%SF_URBAN_PHYSICS == 0) .and. (urbanpt_flag .eqv. .true.)) .or. &
               (NoahmpIO%IVGTYP(I) == NoahmpIO%ISWATER_TABLE) ) then
             NoahmpIO%LAI(I)      = 0.0
             NoahmpIO%XSAIXY(I)   = 0.0
             NoahmpIO%LFMASSXY(I) = 0.0
             NoahmpIO%STMASSXY(I) = 0.0
             NoahmpIO%RTMASSXY(I) = 0.0
             NoahmpIO%WOODXY(I)   = 0.0
             NoahmpIO%STBLCPXY(I) = 0.0
             NoahmpIO%FASTCPXY(I) = 0.0
             NoahmpIO%GRAINXY(I)  = 1.0e-10
             NoahmpIO%GDDXY(I)    = 0
             NoahmpIO%CROPCAT(I)  = 0
          else
             if ( (NoahmpIO%LAI(I) > 100) .or. (NoahmpIO%LAI(I) < 0) ) &
             NoahmpIO%LAI(I)      = 0.0
             NoahmpIO%LAI(I)      = max(NoahmpIO%LAI(I), 0.05)          ! at least start with 0.05 for arbitrary initialization (v3.7)
             NoahmpIO%XSAIXY(I)   = max(0.1*NoahmpIO%LAI(I), 0.05)      ! MB: arbitrarily initialize SAI using input LAI (v3.7)
             if ( urbanpt_flag .eqv. .true. ) then
                NoahmpIO%LFMASSXY(I) = NoahmpIO%LAI(I) * 1000.0 / &
                                       max(NoahmpIO%SLA_TABLE(NoahmpIO%NATURAL_TABLE),1.0) ! use LAI to initialize (v3.7)
             else
                NoahmpIO%LFMASSXY(I) = NoahmpIO%LAI(I) * 1000.0 / &
                                       max(NoahmpIO%SLA_TABLE(NoahmpIO%IVGTYP(I)),1.0)     ! use LAI to initialize (v3.7)
             endif
             NoahmpIO%STMASSXY(I) = NoahmpIO%XSAIXY(I) * 1000.0 / 3.0    ! use SAI to initialize (v3.7)
             NoahmpIO%RTMASSXY(I) = 500.0                                ! these are all arbitrary and probably should be
             NoahmpIO%WOODXY(I)   = 500.0                                ! in the table or read from initialization
             NoahmpIO%STBLCPXY(I) = 1000.0
             NoahmpIO%FASTCPXY(I) = 1000.0
             NoahmpIO%GRAINXY(I)  = 1.0e-10
             NoahmpIO%GDDXY(I)    = 0    

             ! Initialize crop for crop model
             if ( NoahmpIO%IOPT_CROP == 1 ) then
                NoahmpIO%CROPCAT(I) = NoahmpIO%default_crop_table
                if ( NoahmpIO%CROPTYPE(I,5) >= 0.5 ) then
                   NoahmpIO%RTMASSXY(I) = 0.0
                   NoahmpIO%WOODXY  (I) = 0.0
                   if ( (NoahmpIO%CROPTYPE(I,1) > NoahmpIO%CROPTYPE(I,2)) .and. &
                        (NoahmpIO%CROPTYPE(I,1) > NoahmpIO%CROPTYPE(I,3)) .and. &
                        (NoahmpIO%CROPTYPE(I,1) > NoahmpIO%CROPTYPE(I,4)) ) then      ! choose corn
                      NoahmpIO%CROPCAT(I)  = 1
                      NoahmpIO%LFMASSXY(I) = NoahmpIO%LAI(I) / 0.015                  ! Initialize lfmass Zhe Zhang 2020-07-13
                      NoahmpIO%STMASSXY(I) = NoahmpIO%XSAIXY(I) / 0.003
                   elseif ( (NoahmpIO%CROPTYPE(I,2) > NoahmpIO%CROPTYPE(I,1)) .and. &
                            (NoahmpIO%CROPTYPE(I,2) > NoahmpIO%CROPTYPE(I,3)) .and. &
                            (NoahmpIO%CROPTYPE(I,2) > NoahmpIO%CROPTYPE(I,4)) ) then  ! choose soybean
                      NoahmpIO%CROPCAT(I)  = 2
                      NoahmpIO%LFMASSXY(I) = NoahmpIO%LAI(I) / 0.030                  ! Initialize lfmass Zhe Zhang 2020-07-13
                      NoahmpIO%STMASSXY(I) = NoahmpIO%XSAIXY(I) / 0.003
                   else
                      NoahmpIO%CROPCAT(I)  = NoahmpIO%default_crop_table
                      NoahmpIO%LFMASSXY(I) = NoahmpIO%LAI(I) / 0.035
                      NoahmpIO%STMASSXY(I) = NoahmpIO%XSAIXY(I) / 0.003
                   endif
                endif
             endif

             ! Noah-MP irrigation scheme
             if ( (NoahmpIO%IOPT_IRR >= 1) .and. (NoahmpIO%IOPT_IRR <= 3) ) then
                if ( (NoahmpIO%IOPT_IRRM == 0) .or. (NoahmpIO%IOPT_IRRM ==1) ) then       ! sprinkler
                   NoahmpIO%IRNUMSI(I) = 0
                   NoahmpIO%IRWATSI(I) = 0.0
                   NoahmpIO%IRELOSS(I) = 0.0
                   NoahmpIO%IRRSPLH(I) = 0.0    
                elseif ( (NoahmpIO%IOPT_IRRM == 0) .or. (NoahmpIO%IOPT_IRRM == 2) ) then  ! micro or drip
                   NoahmpIO%IRNUMMI(I) = 0
                   NoahmpIO%IRWATMI(I) = 0.0
                   NoahmpIO%IRMIVOL(I) = 0.0
                elseif ( (NoahmpIO%IOPT_IRRM == 0) .or. (NoahmpIO%IOPT_IRRM == 3) ) then  ! flood 
                   NoahmpIO%IRNUMFI(I) = 0
                   NoahmpIO%IRWATFI(I) = 0.0
                   NoahmpIO%IRFIVOL(I) = 0.0
                endif
             endif
          endif
             
          ! initialize soil albedo
          NoahmpIO%ALBSOILDIRXY(I,:) = 0.0
          NoahmpIO%ALBSOILDIFXY(I,:) = 0.0

       enddo ! I
       
       ! Initialize Noah-MP Snow
       call NoahmpSnowinitMain(NoahmpIO)
 
       ! initialize arrays for groundwater dynamics iopt_runsub=5 
       if ( NoahmpIO%IOPT_RUNSUB == 5 ) then
          NoahmpIO%STEPWTD = nint(NoahmpIO%WTDDT * 60.0 / NoahmpIO%DTBL)
          NoahmpIO%STEPWTD = max(NoahmpIO%STEPWTD,1)
       endif

    endif ! NoahmpIO%restart_flag

    if ( NoahmpIO%IOPT_ALB == 3 ) then ! initialize SNICAR aerosol content in snow
       do I = its, ite
          do IZ = -NoahmpIO%NSNOW+1, 0
             if ( (NoahmpIO%SNLIQXY(I,IZ)+NoahmpIO%SNICEXY(I,IZ)) > 0.0 ) then
                NoahmpIO%MassConcBCPHIXY(I,IZ) = NoahmpIO%BCPHIXY(I,IZ) / (NoahmpIO%SNLIQXY(I,IZ) + NoahmpIO%SNICEXY(I,IZ))
                NoahmpIO%MassConcBCPHOXY(I,IZ) = NoahmpIO%BCPHOXY(I,IZ) / (NoahmpIO%SNLIQXY(I,IZ) + NoahmpIO%SNICEXY(I,IZ))
                NoahmpIO%MassConcOCPHIXY(I,IZ) = NoahmpIO%OCPHIXY(I,IZ) / (NoahmpIO%SNLIQXY(I,IZ) + NoahmpIO%SNICEXY(I,IZ))
                NoahmpIO%MassConcOCPHOXY(I,IZ) = NoahmpIO%OCPHOXY(I,IZ) / (NoahmpIO%SNLIQXY(I,IZ) + NoahmpIO%SNICEXY(I,IZ))
                NoahmpIO%MassConcDUST1XY(I,IZ) = NoahmpIO%DUST1XY(I,IZ) / (NoahmpIO%SNLIQXY(I,IZ) + NoahmpIO%SNICEXY(I,IZ))
                NoahmpIO%MassConcDUST2XY(I,IZ) = NoahmpIO%DUST2XY(I,IZ) / (NoahmpIO%SNLIQXY(I,IZ) + NoahmpIO%SNICEXY(I,IZ))
                NoahmpIO%MassConcDUST3XY(I,IZ) = NoahmpIO%DUST3XY(I,IZ) / (NoahmpIO%SNLIQXY(I,IZ) + NoahmpIO%SNICEXY(I,IZ))
                NoahmpIO%MassConcDUST4XY(I,IZ) = NoahmpIO%DUST4XY(I,IZ) / (NoahmpIO%SNLIQXY(I,IZ) + NoahmpIO%SNICEXY(I,IZ))
                NoahmpIO%MassConcDUST5XY(I,IZ) = NoahmpIO%DUST5XY(I,IZ) / (NoahmpIO%SNLIQXY(I,IZ) + NoahmpIO%SNICEXY(I,IZ))
             else
                NoahmpIO%MassConcBCPHIXY(I,IZ) = 0.0
                NoahmpIO%MassConcBCPHOXY(I,IZ) = 0.0
                NoahmpIO%MassConcOCPHIXY(I,IZ) = 0.0
                NoahmpIO%MassConcOCPHOXY(I,IZ) = 0.0
                NoahmpIO%MassConcDUST1XY(I,IZ) = 0.0
                NoahmpIO%MassConcDUST2XY(I,IZ) = 0.0
                NoahmpIO%MassConcDUST3XY(I,IZ) = 0.0
                NoahmpIO%MassConcDUST4XY(I,IZ) = 0.0
                NoahmpIO%MassConcDUST5XY(I,IZ) = 0.0
             endif
          enddo
       enddo
    endif

  end subroutine NoahmpInitMain    

end module NoahmpInitMainMod
