module GroundWaterMmfMod

!!! Module to calculate lateral groundwater flow and the flux between groundwater and rivers
!!! plus the routine to update soil moisture and water table due to those two fluxes
!!! according to the Miguez-Macho & Fan groundwater scheme (Miguez-Macho et al., JGR 2007).
!!! Module written by Gonzalo Miguez-Macho , U. de Santiago de Compostela, Galicia, Spain
!!! November 2012 

! ------------------------ Code history -----------------------------------
! Original Noah-MP subroutine: module_sf_groundwater.F
! Original code: Miguez-Macho&Fan (Miguez-Macho et al 2007, Fan et al 2007)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! Note: this MMF scheme needs further refactoring
! -------------------------------------------------------------------------

  use NoahmpIOVarType
  use NoahmpVarType
  use Machine

   implicit none

contains

  subroutine WTABLE_mmf_noahmp (NoahmpIO)

! ----------------------------------------------------------------------
  implicit none
! ----------------------------------------------------------------------
! IN only

  type(NoahmpIO_type), intent(inout) :: NoahmpIO

!LOCAL  
  
  integer                                                       :: I,J,K,N  
  real(kind=kind_noahmp),  dimension(0:NoahmpIO%NSOIL)          :: ZSOIL                               !depth of soil layer-bottom [m]
  real(kind=kind_noahmp),  dimension(1:NoahmpIO%NSOIL)          :: SMCEQ                               !equilibrium soil water  content [m3/m3]
  real(kind=kind_noahmp),  dimension(1:NoahmpIO%NSOIL)          :: SMC,SH2O
  real(kind=kind_noahmp)                                        :: DELTAT,RCOND,TOTWATER,PSI,        &
                                                                   WFLUXDEEP,WCNDDEEP,DDZ,SMCWTDMID, &
                                                                   WPLUS,WMINUS
  
  real(kind=kind_noahmp)                                        :: BEXP,DKSAT,PSISAT,SMCMAX,SMCWLT
  real(kind=kind_noahmp),  dimension(NoahmpIO%ims:NoahmpIO%ime, NoahmpIO%jms:NoahmpIO%jme) :: WaterTableAvg                       ! grid average WTPD for laterflow subroutine
! --------------------------------------------------------------------------------    
  associate(                                                   &
            ids                => NoahmpIO%ids                ,&
            ide                => NoahmpIO%ide                ,&
            jds                => NoahmpIO%jds                ,&
            jde                => NoahmpIO%jde                ,&  
            kds                => NoahmpIO%kds                ,&
            kde                => NoahmpIO%kde                ,&
            ims                => NoahmpIO%ims                ,&
            ime                => NoahmpIO%ime                ,&  
            jms                => NoahmpIO%jms                ,&
            jme                => NoahmpIO%jme                ,&  
            kms                => NoahmpIO%kms                ,&
            kme                => NoahmpIO%kme                ,&
            its                => NoahmpIO%its                ,&
            ite                => NoahmpIO%ite                ,&  
            jts                => NoahmpIO%jts                ,&
            jte                => NoahmpIO%jte                ,&  
            kts                => NoahmpIO%kts                ,&
            kte                => NoahmpIO%kte                ,& 
            NSOIL              => NoahmpIO%NSOIL              ,&                  
            ISLTYP             => NoahmpIO%ISLTYP             ,&    
            SMOISEQ            => NoahmpIO%SMOISEQ            ,&  ! mosaic 
            DZS                => NoahmpIO%DZS                ,&       
            WTDDT              => NoahmpIO%WTDDT              ,&  
            FDEPTH             => NoahmpIO%FDEPTHXY           ,&     
            AREA               => NoahmpIO%AREAXY             ,&      
            TOPO               => NoahmpIO%TERRAIN            ,&      
            ISURBAN            => NoahmpIO%ISURBAN            ,&   
            IVGTYP             => NoahmpIO%IVGTYP             ,&   
            RIVERCOND          => NoahmpIO%RIVERCONDXY        ,&   
            RIVERBED           => NoahmpIO%RIVERBEDXY         ,&  
            EQWTD              => NoahmpIO%EQZWT              ,&    
            PEXP               => NoahmpIO%PEXPXY             ,&  
            SMOIS              => NoahmpIO%SMOIS              ,&  ! mosaic    
            SH2OXY             => NoahmpIO%SH2O               ,&  ! mosaic  
            SMCWTD             => NoahmpIO%SMCWTDXY           ,&  ! mosaic 
            WTD                => NoahmpIO%ZWTXY              ,&  ! mosaic 
            QLAT               => NoahmpIO%QLATXY             ,&  
            QRF                => NoahmpIO%QRFXY              ,&  
            DEEPRECH           => NoahmpIO%DEEPRECHXY         ,&   
            QSPRING            => NoahmpIO%QSPRINGXY          ,&    
            QSLAT              => NoahmpIO%QSLATXY            ,&    
            QRFS               => NoahmpIO%QRFSXY             ,&  
            QSPRINGS           => NoahmpIO%QSPRINGSXY         ,&   
            RECH               => NoahmpIO%RECHXY             ,&   
            LANDMASK           => NoahmpIO%LANDMASK           ,&
            NumberOfTiles      => NoahmpIO%NumberOfTiles      ,&
            SubGrdFracRescaled => NoahmpIO%SubGrdFracRescaled  &
           )
! -------------------------------------------------------------------------------- 


    WaterTableAvg = 0.0
    DELTAT = WTDDT * 60. !timestep in seconds for this calculation

    ZSOIL(0) = 0.
    ZSOIL(1) = -DZS(1)
    DO K = 2, NSOIL
       ZSOIL(K)         = -DZS(K) + ZSOIL(K-1)
    END DO

! calculate grid average WTD from subgrid WTD
    IF (NoahmpIO%IOPT_MOSAIC .NE. 0) THEN
      DO J=jts,jte
         DO I=its,ite
            IF(LANDMASK(I,J).GT.0)THEN 
              DO N = 1, NumberOfTiles(I,J)
                 WaterTableAvg (I,J) =  WaterTableAvg (I,J) + &
                                        WTD (I,J,N) * SubGrdFracRescaled (I,J,N)
              ENDDO 
            ENDIF
         ENDDO 
      ENDDO   
    ELSE
      WaterTableAvg = WTD(:,:,1)
    ENDIF
!Calculate lateral flow


  QLAT = 0.
  CALL LATERALFLOW (NoahmpIO, ISLTYP,WaterTableAvg         ,&
                   QLAT,FDEPTH,TOPO,LANDMASK,DELTAT,AREA   ,&
                         ids,ide,jds,jde,kds,kde           ,&
                         ims,ime,jms,jme,kms,kme           ,&
                         its,ite,jts,jte,kts,kte            )


!compute flux from grounwater to rivers in the cell

    DO J=jts,jte
       DO I=its,ite
          IF(LANDMASK(I,J).GT.0)THEN
             IF(WaterTableAvg(I,J) .GT. RIVERBED(I,J) .AND.  EQWTD(I,J) .GT. RIVERBED(I,J)) THEN
               RCOND = RIVERCOND(I,J) * EXP(PEXP(I,J)*(WaterTableAvg(I,J)-EQWTD(I,J)))
             ELSE    
               RCOND = RIVERCOND(I,J)       
             ENDIF
             QRF(I,J) = RCOND * (WaterTableAvg(I,J)-RIVERBED(I,J)) * DELTAT/AREA(I,J)
!for now, dont allow it to go from river to groundwater
             QRF(I,J) = MAX(QRF(I,J),0.)
          ELSE
             QRF(I,J) = 0.
          ENDIF
       ENDDO
    ENDDO

    DO J=jts,jte
       DO I=its,ite
          IF(LANDMASK(I,J).GT.0)THEN 
             BEXP   = NoahmpIO%BEXP_TABLE   (ISLTYP(I,J))  ! mosaic added for soil vars
             DKSAT  = NoahmpIO%DKSAT_TABLE  (ISLTYP(I,J))
             PSISAT = -1.0*NoahmpIO%PSISAT_TABLE (ISLTYP(I,J))
             SMCMAX = NoahmpIO%SMCMAX_TABLE (ISLTYP(I,J))
             SMCWLT = NoahmpIO%SMCWLT_TABLE (ISLTYP(I,J))

            IF(IVGTYP(I,J)==NoahmpIO%ISURBAN)THEN
               SMCMAX = 0.45
               SMCWLT = 0.40
            ENDIF

            DO N = 1, NumberOfTiles(I,J)
!for deep water table calculate recharge
               IF(WTD(I,J,N) < ZSOIL(NSOIL)-DZS(NSOIL))THEN
!assume all liquid if the wtd is deep
                  DDZ             = ZSOIL(NSOIL)-WTD(I,J,N)
                  SMCWTDMID       = 0.5 * (SMCWTD(I,J,N) + SMCMAX )
                  PSI             = PSISAT * ( SMCMAX / SMCWTD(I,J,N) ) ** BEXP
                  WCNDDEEP        = DKSAT * ( SMCWTDMID / SMCMAX ) ** (2.0*BEXP + 3.0)
                  WFLUXDEEP       =  - DELTAT * WCNDDEEP * ( (PSISAT-PSI) / DDZ - 1.)
!update deep soil moisture
                  SMCWTD(I,J,N)   = SMCWTD(I,J,N)  + (DEEPRECH(I,J,N) -  WFLUXDEEP)  / DDZ
                  WPLUS           = MAX((SMCWTD(I,J,N)-SMCMAX), 0.0) * DDZ
                  WMINUS          = MAX((1.E-4-SMCWTD(I,J,N)), 0.0) * DDZ
                  SMCWTD(I,J,N)   = MAX( MIN(SMCWTD(I,J,N),SMCMAX) , 1.E-4)
                  WFLUXDEEP       = WFLUXDEEP + WPLUS - WMINUS
                  DEEPRECH(I,J,N) = WFLUXDEEP
               ENDIF


!Total water flux to or from groundwater in the cell
               IF (NoahmpIO%IOPT_MOSAIC .NE. 0) THEN
                  TOTWATER = SubGrdFracRescaled (I,J,N) * (QLAT(I,J) - QRF(I,J)) + DEEPRECH(I,J,N) 
               ELSE 
                  TOTWATER = QLAT(I,J) - QRF(I,J) + DEEPRECH(I,J,N)                             !If mosaic not ON, N=1
               ENDIF

               SMC(1:NSOIL)   = SMOIS(I,1:NSOIL,J,N)
               SH2O(1:NSOIL)  = SH2OXY(I,1:NSOIL,J,N)
               SMCEQ(1:NSOIL) = SMOISEQ(I,1:NSOIL,J,N)

!Update the water table depth and soil moisture
               CALL UPDATEWTD ( NSOIL, DZS , ZSOIL, SMCEQ, SMCMAX, SMCWLT, PSISAT, BEXP, I, J, &!in
                                TOTWATER, WTD(I,J,N), SMC, SH2O, SMCWTD(I,J,N),                &!inout
                                QSPRING(I,J) ) !out

!now update soil moisture
               SMOIS(I,1:NSOIL,J,N)  = SMC(1:NSOIL)
               SH2OXY(I,1:NSOIL,J,N) = SH2O(1:NSOIL)
            ENDDO ! N
          ENDIF
       ENDDO       ! I
    ENDDO          ! J

!accumulate fluxes for output

    DO J=jts,jte
       DO I=its,ite
         IF(LANDMASK(I,J).GT.0)THEN
           QSLAT(I,J) = QSLAT(I,J) + QLAT(I,J)*1.E3
           QRFS(I,J)  = QRFS(I,J) + QRF(I,J)*1.E3
           DO N = 1, NumberOfTiles(I,J)
              QSPRINGS(I,J) = QSPRINGS(I,J) + QSPRING(I,J)*1.E3
              RECH(I,J,N)   = RECH(I,J,N)   + DEEPRECH(I,J,N)*1.E3
              !zero out DEEPRECH
              DEEPRECH(I,J,N) =0.
           ENDDO
         ENDIF
       ENDDO
    ENDDO

    end associate
  end subroutine WTABLE_mmf_noahmp


! ==================================================================================================
! ----------------------------------------------------------------------
  subroutine LATERALFLOW  (NoahmpIO, ISLTYP,WTD,QLAT,FDEPTH,TOPO,LANDMASK,DELTAT,AREA &
                           ,ids,ide,jds,jde,kds,kde                                   &
                           ,ims,ime,jms,jme,kms,kme                                   &
                           ,its,ite,jts,jte,kts,kte                                   )
! ----------------------------------------------------------------------
!  USE NOAHMP_TABLES, ONLY : DKSAT_TABLE

#ifdef MPP_LAND
     use module_mpp_land
#endif
! ----------------------------------------------------------------------
  IMPLICIT NONE
! ----------------------------------------------------------------------
! input

  type(NoahmpIO_type), intent(in)    :: NoahmpIO

  INTEGER,  INTENT(IN   )   ::     ids,ide, jds,jde, kds,kde,  &
       &                           ims,ime, jms,jme, kms,kme,  &
       &                           its,ite, jts,jte, kts,kte
  REAL(kind=kind_noahmp)                                  , INTENT(IN) :: DELTAT
  INTEGER, DIMENSION( ims:ime, jms:jme ), INTENT(IN) :: ISLTYP, LANDMASK
  REAL(kind=kind_noahmp),    DIMENSION( ims:ime, jms:jme ), INTENT(IN) :: FDEPTH,WTD,TOPO,AREA

!output
  REAL(kind=kind_noahmp), DIMENSION( ims:ime , jms:jme ), INTENT(OUT) :: QLAT

!local
  INTEGER                              :: I, J, itsh, iteh, jtsh, jteh, nx, ny
  REAL(kind=kind_noahmp)               :: Q, KLAT

#ifdef MPP_LAND 
  ! halo'ed arrays
!gmm this should not be like this if memory dimensions are already different than tile dimensions
  REAL(kind=kind_noahmp),    DIMENSION(ims-1:ime+1, jms-1:jme+1) :: KCELL, HEAD
#else
  REAL(kind=kind_noahmp),    DIMENSION(ims:ime, jms:jme) :: KCELL, HEAD
#endif

  REAL(kind=kind_noahmp), DIMENSION(19)      :: KLATFACTOR
  DATA KLATFACTOR /2.,3.,4.,10.,10.,12.,14.,20.,24.,28.,40.,48.,2.,0.,10.,0.,20.,2.,2./

  REAL(kind=kind_noahmp),    PARAMETER :: PI = 3.14159265 
  REAL(kind=kind_noahmp),    PARAMETER :: FANGLE = 0.22754493   ! = 0.5*sqrt(0.5*tan(pi/8))


#ifdef MPP_LAND
    itsh=ims
    iteh=ime
    jtsh=jms
    jteh=jme
#else
    itsh=its
    iteh=ite
    jtsh=jts
    jteh=jte
#endif

    DO J=jtsh,jteh
       DO I=itsh,iteh
           IF(FDEPTH(I,J).GT.0.)THEN
                 KLAT = NoahmpIO%DKSAT_TABLE(ISLTYP(I,J)) * KLATFACTOR(ISLTYP(I,J))
                 IF(WTD(I,J) < -1.5)THEN
                     KCELL(I,J) = FDEPTH(I,J) * KLAT * EXP( (WTD(I,J) + 1.5) / FDEPTH(I,J) )
                 ELSE
                     KCELL(I,J) = KLAT * ( WTD(I,J) + 1.5 + FDEPTH(I,J) )
                 ENDIF
           ELSE
                 KCELL(i,J) = 0.
           ENDIF

           HEAD(I,J) = TOPO(I,J) + WTD(I,J)
       ENDDO
    ENDDO

#ifdef MPP_LAND

  nx = ((ime-ims) + 1) + 2      ! include halos
  ny = ((jme-jms) + 1) + 2      ! include halos
  
!copy neighbor's values for haloed variables
!gmm if memory dimensions already included the halo from the start, only wtd would need to be communicated
!gmm before the kcell and head calculation,
!gmm since the communication of the static arrays could be done at initial time

!first do left and right communication
  call mpp_land_lr_com(kcell, nx, ny, 99)
  call mpp_land_lr_com(head, nx, ny, 99)

  call mpp_land_sync() !need this to make sure that the corners or haloes are passed

!once all know their left and right haloes, do the up and down
  call mpp_land_ub_com(kcell, nx, ny, 99)
  call mpp_land_ub_com(head, nx, ny, 99)

#endif

#ifdef MPP_LAND
    itsh=max(its,2)
    iteh=min(ite,global_nx-1)
    jtsh=max(jts,2)
    jteh=min(jte,global_ny-1)
#else
    itsh=max(its,ids+1)
    iteh=min(ite,ide-1)
    jtsh=max(jts,jds+1)
    jteh=min(jte,jde-1)
#endif

    DO J=jtsh,jteh
       DO I=itsh,iteh
          IF( LANDMASK(I,J).GT.0   )THEN
                 Q=0.
                             
                 Q  = Q + (KCELL(I-1,J+1)+KCELL(I,J)) &
                        * (HEAD(I-1,J+1)-HEAD(I,J))/SQRT(2.)
                             
                 Q  = Q +  (KCELL(I-1,J)+KCELL(I,J)) &
                        *  (HEAD(I-1,J)-HEAD(I,J))

                 Q  = Q +  (KCELL(I-1,J-1)+KCELL(I,J)) &
                        * (HEAD(I-1,J-1)-HEAD(I,J))/SQRT(2.)

                 Q  = Q +  (KCELL(I,J+1)+KCELL(I,J)) &
                        * (HEAD(I,J+1)-HEAD(I,J))

                 Q  = Q +  (KCELL(I,J-1)+KCELL(I,J)) &
                        * (HEAD(I,J-1)-HEAD(I,J))

                 Q  = Q +  (KCELL(I+1,J+1)+KCELL(I,J)) &
                        * (HEAD(I+1,J+1)-HEAD(I,J))/SQRT(2.)
  
                 Q  = Q +  (KCELL(I+1,J)+KCELL(I,J)) &
                        * (HEAD(I+1,J)-HEAD(I,J))

                 Q  = Q +  (KCELL(I+1,J-1)+KCELL(I,J)) &
                        * (HEAD(I+1,J-1)-HEAD(I,J))/SQRT(2.)

                 ! Here, Q is in m3/s. To convert to m, divide it by area of the grid cell.
                 QLAT(I,J) = FANGLE* Q * DELTAT / AREA(I,J)
          ENDIF
       ENDDO
    ENDDO

  end subroutine LATERALFLOW


! ==================================================================================================
! ----------------------------------------------------------------------
  subroutine UPDATEWTD  (NSOIL,  DZS,  ZSOIL ,SMCEQ                ,& !in
                         SMCMAX, SMCWLT, PSISAT, BEXP ,ILOC ,JLOC  ,& !in
                         TOTWATER, WTD ,SMC, SH2O ,SMCWTD          ,& !inout
                         QSPRING                                 )  !out
! ----------------------------------------------------------------------
  IMPLICIT NONE
! ----------------------------------------------------------------------
! input
  INTEGER,                         INTENT(IN) :: NSOIL !no. of soil layers
  INTEGER,                         INTENT(IN) :: ILOC, JLOC
  REAL(kind=kind_noahmp),                         INTENT(IN)    :: SMCMAX
  REAL(kind=kind_noahmp),                         INTENT(IN)    :: SMCWLT
  REAL(kind=kind_noahmp),                         INTENT(IN)    :: PSISAT
  REAL(kind=kind_noahmp),                         INTENT(IN)    :: BEXP
  REAL(kind=kind_noahmp),  DIMENSION(       0:NSOIL), INTENT(IN) :: ZSOIL !depth of soil layer-bottom [m]
  REAL(kind=kind_noahmp),  DIMENSION(       1:NSOIL), INTENT(IN) :: SMCEQ  !equilibrium soil water  content [m3/m3]
  REAL(kind=kind_noahmp),  DIMENSION(       1:NSOIL), INTENT(IN) :: DZS ! soil layer thickness [m]
! input-output
  REAL(kind=kind_noahmp)                           , INTENT(INOUT) :: TOTWATER
  REAL(kind=kind_noahmp)                           , INTENT(INOUT) :: WTD
  REAL(kind=kind_noahmp)                           , INTENT(INOUT) :: SMCWTD
  REAL(kind=kind_noahmp), DIMENSION(       1:NSOIL), INTENT(INOUT) :: SMC
  REAL(kind=kind_noahmp), DIMENSION(       1:NSOIL), INTENT(INOUT) :: SH2O
! output
  REAL(kind=kind_noahmp)                           , INTENT(OUT) :: QSPRING
!local
  INTEGER                                     :: K
  INTEGER                                     :: K1
  INTEGER                                     :: IWTD
  INTEGER                                     :: KWTD
  REAL(kind=kind_noahmp)                                        :: MAXWATUP, MAXWATDW ,WTDOLD
  REAL(kind=kind_noahmp)                                        :: WGPMID
  REAL(kind=kind_noahmp)                                        :: SYIELDDW
  REAL(kind=kind_noahmp)                                        :: DZUP
  REAL(kind=kind_noahmp)                                        :: SMCEQDEEP
  REAL(kind=kind_noahmp), DIMENSION(       1:NSOIL)             :: SICE
! -------------------------------------------------------------



  QSPRING=0.

  SICE = SMC - SH2O

iwtd=1

!case 1: totwater > 0 (water table going up):
IF(totwater.gt.0.)then


         if(wtd.ge.zsoil(nsoil))then

            do k=nsoil-1,1,-1
              if(wtd.lt.zsoil(k))exit
            enddo
            iwtd=k
            kwtd=iwtd+1

!max water that fits in the layer
            maxwatup=dzs(kwtd)*(smcmax-smc(kwtd))

            if(totwater.le.maxwatup)then
               smc(kwtd) = smc(kwtd) + totwater / dzs(kwtd)
               smc(kwtd) = min(smc(kwtd),smcmax)
               if(smc(kwtd).gt.smceq(kwtd))wtd = min ( ( smc(kwtd)*dzs(kwtd) &
                 - smceq(kwtd)*zsoil(iwtd) + smcmax*zsoil(kwtd) ) / &
                     ( smcmax-smceq(kwtd) ) , zsoil(iwtd) )
               totwater=0.
            else   !water enough to saturate the layer
              smc(kwtd) = smcmax
              totwater=totwater-maxwatup
              k1=iwtd
              do k=k1,0,-1
                 wtd = zsoil(k)
                 iwtd=k-1
                 if(k.eq.0)exit
                 maxwatup=dzs(k)*(smcmax-smc(k))
                 if(totwater.le.maxwatup)then
                   smc(k) = smc(k) + totwater / dzs(k)
                   smc(k) = min(smc(k),smcmax)
                   if(smc(k).gt.smceq(k))wtd = min ( ( smc(k)*dzs(k) &
                     - smceq(k)*zsoil(iwtd) + smcmax*zsoil(k) ) / &
                     ( smcmax-smceq(k) ) , zsoil(iwtd) )
                   totwater=0.
                   exit
                 else
                    smc(k) = smcmax
                    totwater=totwater-maxwatup
                 endif

              enddo

            endif

         elseif(wtd.ge.zsoil(nsoil)-dzs(nsoil))then ! wtd below bottom of soil model

            !gmmequilibrium soil moisture content
               smceqdeep = smcmax * ( psisat / &
                           (psisat - dzs(nsoil)) ) ** (1./bexp)
!               smceqdeep = max(smceqdeep,smcwlt)
               smceqdeep = max(smceqdeep,1.E-4)

            maxwatup=(smcmax-smcwtd)*dzs(nsoil)

            if(totwater.le.maxwatup)then
                smcwtd = smcwtd + totwater / dzs(nsoil)
                smcwtd = min(smcwtd,smcmax)
                if(smcwtd.gt.smceqdeep)wtd = min( ( smcwtd*dzs(nsoil) &
                 - smceqdeep*zsoil(nsoil) + smcmax*(zsoil(nsoil)-dzs(nsoil)) ) / &
                     ( smcmax-smceqdeep ) , zsoil(nsoil) )
                totwater=0.
            else
                smcwtd=smcmax
                totwater=totwater-maxwatup
                do k=nsoil,0,-1
                    wtd=zsoil(k)
                    iwtd=k-1
                    if(k.eq.0)exit
                    maxwatup=dzs(k)*(smcmax-smc(k))
                    if(totwater.le.maxwatup)then
                     smc(k) = min(smc(k) + totwater / dzs(k),smcmax)
                     if(smc(k).gt.smceq(k))wtd = min ( ( smc(k)*dzs(k) &
                        - smceq(k)*zsoil(iwtd) + smcmax*zsoil(k) ) / &
                           ( smcmax-smceq(k) ) , zsoil(iwtd) )
                     totwater=0.
                     exit
                    else
                     smc(k) = smcmax
                     totwater=totwater-maxwatup
                    endif
                enddo
             endif

!deep water table
       else

            maxwatup=(smcmax-smcwtd)*(zsoil(nsoil)-dzs(nsoil)-wtd)
            if(totwater.le.maxwatup)then
               wtd = wtd + totwater/(smcmax-smcwtd)
               totwater=0.
            else
               totwater=totwater-maxwatup
               wtd=zsoil(nsoil)-dzs(nsoil)
               maxwatup=(smcmax-smcwtd)*dzs(nsoil)
              if(totwater.le.maxwatup)then

            !gmmequilibrium soil moisture content
               smceqdeep = smcmax * ( psisat / &
                           (psisat - dzs(nsoil)) ) ** (1./bexp)
!               smceqdeep = max(smceqdeep,smcwlt)
               smceqdeep = max(smceqdeep,1.E-4)

                smcwtd = smcwtd + totwater / dzs(nsoil)
                smcwtd = min(smcwtd,smcmax)
                wtd = ( smcwtd*dzs(nsoil) &
                 - smceqdeep*zsoil(nsoil) + smcmax*(zsoil(nsoil)-dzs(nsoil)) ) / &
                     ( smcmax-smceqdeep )
                totwater=0.
              else
                smcwtd=smcmax
                totwater=totwater-maxwatup
                do k=nsoil,0,-1
                    wtd=zsoil(k)
                    iwtd=k-1
                    if(k.eq.0)exit
                    maxwatup=dzs(k)*(smcmax-smc(k))

                    if(totwater.le.maxwatup)then
                     smc(k) = smc(k) + totwater / dzs(k)
                     smc(k) = min(smc(k),smcmax)
                     if(smc(k).gt.smceq(k))wtd = ( smc(k)*dzs(k) &
                        - smceq(k)*zsoil(iwtd) + smcmax*zsoil(k) ) / &
                           ( smcmax-smceq(k) )
                     totwater=0.
                     exit
                    else
                     smc(k) = smcmax
                     totwater=totwater-maxwatup
                    endif
                   enddo
               endif
             endif
         endif

!water springing at the surface
        qspring=totwater

!case 2: totwater < 0 (water table going down):
ELSEIF(totwater.lt.0.)then


         if(wtd.ge.zsoil(nsoil))then !wtd in the resolved layers

            do k=nsoil-1,1,-1
               if(wtd.lt.zsoil(k))exit
            enddo
            iwtd=k

               k1=iwtd+1
               do kwtd=k1,nsoil

!max water that the layer can yield
                  maxwatdw=dzs(kwtd)*(smc(kwtd)-max(smceq(kwtd),sice(kwtd)))

                  if(-totwater.le.maxwatdw)then
                        smc(kwtd) = smc(kwtd) + totwater / dzs(kwtd)
                        if(smc(kwtd).gt.smceq(kwtd))then
                              wtd = ( smc(kwtd)*dzs(kwtd) &
                                 - smceq(kwtd)*zsoil(iwtd) + smcmax*zsoil(kwtd) ) / &
                                 ( smcmax-smceq(kwtd) )
                         else
                              wtd=zsoil(kwtd)
                              iwtd=iwtd+1
                         endif
                         totwater=0.
                         exit
                   else
                         wtd = zsoil(kwtd)
                         iwtd=iwtd+1
                         if(maxwatdw.ge.0.)then
                            smc(kwtd) = smc(kwtd) + maxwatdw / dzs(kwtd)
                            totwater = totwater + maxwatdw
                         endif
                   endif

                enddo

               if(iwtd.eq.nsoil.and.totwater.lt.0.)then
            !gmmequilibrium soil moisture content
               smceqdeep = smcmax * ( psisat / &
                           (psisat - dzs(nsoil)) ) ** (1./bexp)
!               smceqdeep = max(smceqdeep,smcwlt)
               smceqdeep = max(smceqdeep,1.E-4)

                  maxwatdw=dzs(nsoil)*(smcwtd-smceqdeep)

                  if(-totwater.le.maxwatdw)then

                       smcwtd = smcwtd + totwater / dzs(nsoil)
                       wtd = max( ( smcwtd*dzs(nsoil) &
                           - smceqdeep*zsoil(nsoil) + smcmax*(zsoil(nsoil)-dzs(nsoil)) ) / &
                            ( smcmax-smceqdeep ) , zsoil(nsoil)-dzs(nsoil) )

                  else

                       wtd=zsoil(nsoil)-dzs(nsoil)
                       smcwtd = smcwtd + totwater / dzs(nsoil)
!and now even further down
                       dzup=(smceqdeep-smcwtd)*dzs(nsoil)/(smcmax-smceqdeep)
                       wtd=wtd-dzup
                       smcwtd=smceqdeep

                  endif

                endif



        elseif(wtd.ge.zsoil(nsoil)-dzs(nsoil))then

!if wtd was already below the bottom of the resolved soil crust
            !gmmequilibrium soil moisture content
               smceqdeep = smcmax * ( psisat / &
                           (psisat - dzs(nsoil)) ) ** (1./bexp)
!               smceqdeep = max(smceqdeep,smcwlt)
               smceqdeep = max(smceqdeep,1.E-4)

            maxwatdw=dzs(nsoil)*(smcwtd-smceqdeep)

            if(-totwater.le.maxwatdw)then

               smcwtd = smcwtd + totwater / dzs(nsoil)
               wtd = max( ( smcwtd*dzs(nsoil) &
                    - smceqdeep*zsoil(nsoil) + smcmax*(zsoil(nsoil)-dzs(nsoil)) ) / &
                    ( smcmax-smceqdeep ) , zsoil(nsoil)-dzs(nsoil) )

            else

               wtd=zsoil(nsoil)-dzs(nsoil)
               smcwtd = smcwtd + totwater / dzs(nsoil)
!and now even further down
               dzup=(smceqdeep-smcwtd)*dzs(nsoil)/(smcmax-smceqdeep)
               wtd=wtd-dzup
               smcwtd=smceqdeep

             endif

         else
!gmmequilibrium soil moisture content
               wgpmid = smcmax * ( psisat / &
                    (psisat - (zsoil(nsoil)-wtd)) ) ** (1./bexp)
!               wgpmid=max(wgpmid,smcwlt)
               wgpmid=max(wgpmid,1.E-4)
               syielddw=smcmax-wgpmid
               wtdold=wtd
               wtd = wtdold + totwater/syielddw
!update wtdwgp
               smcwtd = (smcwtd*(zsoil(nsoil)-wtdold)+wgpmid*(wtdold-wtd) ) / (zsoil(nsoil)-wtd)

          endif

          qspring=0.

ENDIF

         SH2O = SMC - SICE


  end subroutine UPDATEWTD

! ----------------------------------------------------------------------

END MODULE GroundWaterMmfMod
