!>  \file noahmp_tables.f90
!!  This file contains Fortran versions of the data tables included with NoahMP in mptable.tbl, soilparm.tbl, and genparm.tbl. 

!> \ingroup NoahMP_LSM
!! \brief Data from MPTABLE.TBL, SOILPARM.TBL, GENPARM.TBL for NoahMP
!!
!! Note that a subset of the data in the *.TBL files is represented in this file. For example,
!! only the data in the noah_mp_modis_parameters section of MPTABLE.TBL and the STAS section of
!! SOILPARM.TBL are included in this module.
module noahmp_tables
use machine ,   only : kind_phys
    implicit none

    integer, private, parameter :: mvt   = 30 ! use 30 instead of 27
    integer, private, parameter :: mband = 2
    integer, private, parameter :: msc   = 20
    integer, private, parameter :: max_soiltyp = 30
    integer, private, parameter :: ncrop = 5
    integer, private, parameter :: nstage = 8
    integer, private, parameter :: num_slope = 9

! mptable.tbl vegetation parameters

    integer :: isurban_table
    integer :: iswater_table 
    integer :: isbarren_table
    integer :: isice_table 
    integer :: iscrop_table 
    integer :: eblforest_table
    integer :: natural_table
    integer :: lcz_1_table
    integer :: lcz_2_table
    integer :: lcz_3_table
    integer :: lcz_4_table
    integer :: lcz_5_table
    integer :: lcz_6_table
    integer :: lcz_7_table
    integer :: lcz_8_table
    integer :: lcz_9_table
    integer :: lcz_10_table
    integer :: lcz_11_table
    real (kind=kind_phys) :: ch2op_table(mvt)       !< maximum intercepted h2o per unit lai+sai (mm)
    real (kind=kind_phys) :: dleaf_table(mvt)       !< characteristic leaf dimension (m)
    real (kind=kind_phys) :: z0mvt_table(mvt)       !< momentum roughness length (m)
    real (kind=kind_phys) :: hvt_table(mvt)         !< top of canopy (m)
    real (kind=kind_phys) :: hvb_table(mvt)         !< bottom of canopy (m)
    real (kind=kind_phys) :: z0mhvt_table(mvt)      !< ratio of z0m to hvt
    real (kind=kind_phys) :: den_table(mvt)         !< tree density (no. of trunks per m2)
    real (kind=kind_phys) :: rc_table(mvt)          !< tree crown radius (m)
    real (kind=kind_phys) :: mfsno_table(mvt)       !< snowmelt curve parameter ()
    real (kind=kind_phys) :: scffac_table(mvt)      !< snow cover factor (m)
    real (kind=kind_phys) :: cbiom_table(mvt)       !< canopy biomass heat capacity parameter (m)
    real (kind=kind_phys) :: saim_table(mvt,12)     !< monthly stem area index, one-sided
    real (kind=kind_phys) :: laim_table(mvt,12)     !< monthly leaf area index, one-sided
    real (kind=kind_phys) :: sla_table(mvt)         !< single-side leaf area per kg [m2/kg]
    real (kind=kind_phys) :: dilefc_table(mvt)      !< coeficient for leaf stress death [1/s]
    real (kind=kind_phys) :: dilefw_table(mvt)      !< coeficient for leaf stress death [1/s]
    real (kind=kind_phys) :: fragr_table(mvt)       !< fraction of growth respiration  !original was 0.3 
    real (kind=kind_phys) :: ltovrc_table(mvt)      !< leaf turnover [1/s]

    real (kind=kind_phys) :: c3psn_table(mvt)       !< photosynthetic pathway: 0. = c4, 1. = c3
    real (kind=kind_phys) :: kc25_table(mvt)        !< co2 michaelis-menten constant at 25c (pa)
    real (kind=kind_phys) :: akc_table(mvt)         !< q10 for kc25
    real (kind=kind_phys) :: ko25_table(mvt)        !< o2 michaelis-menten constant at 25c (pa)
    real (kind=kind_phys) :: ako_table(mvt)         !< q10 for ko25
    real (kind=kind_phys) :: vcmx25_table(mvt)      !< maximum rate of carboxylation at 25c (umol co2/m**2/s)
    real (kind=kind_phys) :: avcmx_table(mvt)       !< q10 for vcmx25
    real (kind=kind_phys) :: bp_table(mvt)          !< minimum leaf conductance (umol/m**2/s)
    real (kind=kind_phys) :: mp_table(mvt)          !< slope of conductance-to-photosynthesis relationship
    real (kind=kind_phys) :: qe25_table(mvt)        !< quantum efficiency at 25c (umol co2 / umo photon)
    real (kind=kind_phys) :: aqe_table(mvt)         !< q10 for qe25
    real (kind=kind_phys) :: rmf25_table(mvt)       !< leaf maintenance respiration at 25c (umol co2/m**2/s)
    real (kind=kind_phys) :: rms25_table(mvt)       !< stem maintenance respiration at 25c (umol co2/kg bio/s)
    real (kind=kind_phys) :: rmr25_table(mvt)       !< root maintenance respiration at 25c (umol co2/kg bio/s)
    real (kind=kind_phys) :: arm_table(mvt)         !< q10 for maintenance respiration
    real (kind=kind_phys) :: folnmx_table(mvt)      !< foliage nitrogen concentration when f(n)=1 (%)
    real (kind=kind_phys) :: tmin_table(mvt)        !< minimum temperature for photosynthesis (k)

    real (kind=kind_phys) :: xl_table(mvt)          !< leaf/stem orientation index
    real (kind=kind_phys) :: rhol_table(mvt,mband)  !< leaf reflectance: 1=vis, 2=nir
    real (kind=kind_phys) :: rhos_table(mvt,mband)  !< stem reflectance: 1=vis, 2=nir
    real (kind=kind_phys) :: taul_table(mvt,mband)  !< leaf transmittance: 1=vis, 2=nir
    real (kind=kind_phys) :: taus_table(mvt,mband)  !< stem transmittance: 1=vis, 2=nir

    real (kind=kind_phys) :: mrp_table(mvt)         !< microbial respiration parameter (umol co2 /kg c/ s)
    real (kind=kind_phys) :: cwpvt_table(mvt)       !< empirical canopy wind parameter

    real (kind=kind_phys) :: wrrat_table(mvt)       !< wood to non-wood ratio
    real (kind=kind_phys) :: wdpool_table(mvt)      !< wood pool (switch 1 or 0) depending on woody or not [-]
    real (kind=kind_phys) :: tdlef_table(mvt)       !< characteristic t for leaf freezing [k]

    real (kind=kind_phys) :: nroot_table(mvt)       !< number of soil layers with root present
    real (kind=kind_phys) :: rgl_table(mvt)         !< parameter used in radiation stress function
    real (kind=kind_phys) :: rs_table(mvt)          !< minimum stomatal resistance [s m-1]
    real (kind=kind_phys) :: hs_table(mvt)          !< parameter used in vapor pressure deficit function
    real (kind=kind_phys) :: topt_table(mvt)        !< optimum transpiration air temperature [k]
    real (kind=kind_phys) :: rsmax_table(mvt)       !< maximal stomatal resistance [s m-1]

! soilparm.tbl parameters

    integer :: slcats

    real (kind=kind_phys) :: bexp_table(max_soiltyp)   
    real (kind=kind_phys) :: smcdry_table(max_soiltyp)  
    real (kind=kind_phys) :: f1_table(max_soiltyp)     
    real (kind=kind_phys) :: smcmax_table(max_soiltyp)
    real (kind=kind_phys) :: smcref_table(max_soiltyp)  
    real (kind=kind_phys) :: psisat_table(max_soiltyp) 
    real (kind=kind_phys) :: dksat_table(max_soiltyp) 
    real (kind=kind_phys) :: dwsat_table(max_soiltyp)
    real (kind=kind_phys) :: smcwlt_table(max_soiltyp)   
    real (kind=kind_phys) :: quartz_table(max_soiltyp)  
    real (kind=kind_phys) :: bvic_table(max_soiltyp)        !vic model infiltration parameter (-) for opt_run=6
    real (kind=kind_phys) :: axaj_table(max_soiltyp)        !Xinanjiang: Tension water distribution inflection parameter [-] for opt_run=7
    real (kind=kind_phys) :: bxaj_table(max_soiltyp)        !Xinanjiang: Tension water distribution shape parameter [-] for opt_run=7
    real (kind=kind_phys) :: xxaj_table(max_soiltyp)        !Xinanjiang: Free water distribution shape parameter [-] for opt_run=7
    real (kind=kind_phys) :: bdvic_table(max_soiltyp)       !VIC model infiltration parameter (-)
    real (kind=kind_phys) :: gdvic_table(max_soiltyp)       !mean capilary drive (m)
    real (kind=kind_phys) :: bbvic_table(max_soiltyp)       !heterogeniety parameter for DVIC infiltration [-]

! genparm.tbl parameters

    real (kind=kind_phys) :: slope_table(num_slope)                     !< slope factor for soil drainage
    
    real (kind=kind_phys) :: csoil_table       !< soil heat capacity [j m-3 k-1]
    real (kind=kind_phys) :: refdk_table       !< parameter in the surface runoff parameterization
    real (kind=kind_phys) :: refkdt_table      !< parameter in the surface runoff parameterization
    real (kind=kind_phys) :: frzk_table        !< frozen ground parameter
    real (kind=kind_phys) :: zbot_table        !< depth [m] of lower boundary soil temperature
    real (kind=kind_phys) :: czil_table        !< parameter used in the calculation of the roughness length for heat

! mptable.tbl radiation parameters

    real (kind=kind_phys) :: albsat_table(msc,mband)   !< saturated soil albedos: 1=vis, 2=nir
    real (kind=kind_phys) :: albdry_table(msc,mband)   !< dry soil albedos: 1=vis, 2=nir
    real (kind=kind_phys) :: albice_table(mband)       !< albedo land ice: 1=vis, 2=nir
    real (kind=kind_phys) :: alblak_table(mband)       !< albedo frozen lakes: 1=vis, 2=nir
    real (kind=kind_phys) :: omegas_table(mband)       !< two-stream parameter omega for snow
    real (kind=kind_phys) :: betads_table              !< two-stream parameter betad for snow
    real (kind=kind_phys) :: betais_table              !< two-stream parameter betad for snow
    real (kind=kind_phys) :: eg_table(2)               !< emissivity

! mptable.tbl global parameters
 
    real (kind=kind_phys) :: co2_table                 !< co2 partial pressure
    real (kind=kind_phys) :: o2_table                  !< o2 partial pressure
    real (kind=kind_phys) :: timean_table              !< gridcell mean topgraphic index (global mean)
    real (kind=kind_phys) :: fsatmx_table              !< maximum surface saturated fraction (global mean)
    real (kind=kind_phys) :: z0sno_table               !< snow surface roughness length (m) (0.002)
    real (kind=kind_phys) :: ssi_table                 !< liquid water holding capacity for snowpack (m3/m3) (0.03)
    real (kind=kind_phys) :: snow_ret_fac_table        !< snowpack water release timescale factor (1/s)
    real (kind=kind_phys) :: snow_emis_table           !< surface emissivity
    real (kind=kind_phys) :: swemx_table               !< new snow mass to fully cover old snow (mm)
    real (kind=kind_phys) :: tau0_table                !< tau0 from yang97 eqn. 10a
    real (kind=kind_phys) :: grain_growth_table        !< growth from vapor diffusion yang97 eqn. 10b
    real (kind=kind_phys) :: extra_growth_table        !< extra growth near freezing yang97 eqn. 10c
    real (kind=kind_phys) :: dirt_soot_table           !< dirt and soot term yang97 eqn. 10d
    real (kind=kind_phys) :: bats_cosz_table           !< zenith angle snow albedo adjustment; b in yang97 eqn. 15
    real (kind=kind_phys) :: bats_vis_new_table        !< new snow visible albedo
    real (kind=kind_phys) :: bats_nir_new_table        !< new snow nir albedo
    real (kind=kind_phys) :: bats_vis_age_table        !< age factor for diffuse visible snow albedo yang97 eqn. 17
    real (kind=kind_phys) :: bats_nir_age_table        !< age factor for diffuse nir snow albedo yang97 eqn. 18
    real (kind=kind_phys) :: bats_vis_dir_table        !< cosz factor for direct visible snow albedo yang97 eqn. 15
    real (kind=kind_phys) :: bats_nir_dir_table        !< cosz factor for direct nir snow albedo yang97 eqn. 16
    real (kind=kind_phys) :: rsurf_snow_table          !< surface resistance for snow(s/m)
    real (kind=kind_phys) :: rsurf_exp_table           !< exponent in the shape parameter for soil resistance option 1

! mptable.tbl irrigation parameters

    real (kind=kind_phys) :: irr_frac_table              ! irrigation Fraction
 integer :: irr_har_table               ! number of days before harvest date to stop irrigation
    real (kind=kind_phys) :: irr_lai_table               ! Minimum lai to trigger irrigation
    real (kind=kind_phys) :: irr_mad_table               ! management allowable deficit (0-1)
    real (kind=kind_phys) :: filoss_table                ! fraction of flood irrigation loss (0-1)
    real (kind=kind_phys) :: sprir_rate_table            ! mm/h, sprinkler irrigation rate
    real (kind=kind_phys) :: micir_rate_table            ! mm/h, micro irrigation rate
    real (kind=kind_phys) :: firtfac_table               ! flood application rate factor
    real (kind=kind_phys) :: ir_rain_table               ! maximum precipitation to stop irrigation trigger

! mptable.tbl crop parameters

    integer :: default_crop_table          ! Default crop index
    integer :: pltday_table(ncrop)         !< planting date
    integer :: hsday_table(ncrop)          !< harvest date
    real (kind=kind_phys) :: plantpop_table(ncrop)       !< plant density [per ha] - used?
    real (kind=kind_phys) :: irri_table(ncrop)           !< irrigation strategy 0= non-irrigation 1=irrigation (no water-stress)

    real (kind=kind_phys) :: gddtbase_table(ncrop)       !< base temperature for gdd accumulation [c]
    real (kind=kind_phys) :: gddtcut_table(ncrop)        !< upper temperature for gdd accumulation [c]
    real (kind=kind_phys) :: gdds1_table(ncrop)          !< gdd from seeding to emergence
    real (kind=kind_phys) :: gdds2_table(ncrop)          !< gdd from seeding to initial vegetative 
    real (kind=kind_phys) :: gdds3_table(ncrop)          !< gdd from seeding to  post vegetative
    real (kind=kind_phys) :: gdds4_table(ncrop)          !< gdd from seeding to  intial reproductive
    real (kind=kind_phys) :: gdds5_table(ncrop)          !< gdd from seeding to pysical maturity 

    real (kind=kind_phys) :: c3psni_table(ncrop)       !photosynthetic pathway: 0. = c4, 1. = c3 ! Zhe Zhang 2020-07-03
    real (kind=kind_phys) :: kc25i_table(ncrop)        !co2 michaelis-menten constant at 25c (pa)
    real (kind=kind_phys) :: akci_table(ncrop)         !q10 for kc25
    real (kind=kind_phys) :: ko25i_table(ncrop)        !o2 michaelis-menten constant at 25c (pa)
    real (kind=kind_phys) :: akoi_table(ncrop)         !q10 for ko25
    real (kind=kind_phys) :: vcmx25i_table(ncrop)      !maximum rate of carboxylation at 25c (umol co2/m**2/s)
    real (kind=kind_phys) :: avcmxi_table(ncrop)       !q10 for vcmx25
    real (kind=kind_phys) :: bpi_table(ncrop)          !minimum leaf conductance (umol/m**2/s)
    real (kind=kind_phys) :: mpi_table(ncrop)          !slope of conductance-to-photosynthesis relationship
    real (kind=kind_phys) :: qe25i_table(ncrop)        !quantum efficiency at 25c (umol co2 / umol photon)
    real (kind=kind_phys) :: folnmxi_table(ncrop)      !foliage nitrogen concentration when

    integer :: c3c4_table(ncrop)           !< photosynthetic pathway:  1. = c3 2. = c4
    real (kind=kind_phys) :: aref_table(ncrop)           !< reference maximum co2 assimulation rate 
    real (kind=kind_phys) :: psnrf_table(ncrop)          !< co2 assimulation reduction factor(0-1) (caused by non-modeling part,e.g.pest,weeds)
    real (kind=kind_phys) :: i2par_table(ncrop)          !< fraction of incoming solar radiation to photosynthetically active radiation
    real (kind=kind_phys) :: tassim0_table(ncrop)        !< minimum temperature for co2 assimulation [c]
    real (kind=kind_phys) :: tassim1_table(ncrop)        !< co2 assimulation linearly increasing until temperature reaches t1 [c]
    real (kind=kind_phys) :: tassim2_table(ncrop)        !< co2 assmilation rate remain at aref until temperature reaches t2 [c]
    real (kind=kind_phys) :: k_table(ncrop)              !< light extinction coefficient
    real (kind=kind_phys) :: epsi_table(ncrop)           !< initial light use efficiency

    real (kind=kind_phys) :: q10mr_table(ncrop)          !< q10 for maintainance respiration
    real (kind=kind_phys) :: foln_mx_table(ncrop)        !< foliage nitrogen concentration when f(n)=1 (%)
    real (kind=kind_phys) :: lefreez_table(ncrop)        !< characteristic t for leaf freezing [k]

    real (kind=kind_phys) :: dile_fc_table(ncrop,nstage) !< coeficient for temperature leaf stress death [1/s]
    real (kind=kind_phys) :: dile_fw_table(ncrop,nstage) !< coeficient for water leaf stress death [1/s]
    real (kind=kind_phys) :: fra_gr_table(ncrop)         !< fraction of growth respiration

    real (kind=kind_phys) :: lf_ovrc_table(ncrop,nstage) !< fraction of leaf turnover  [1/s]
    real (kind=kind_phys) :: st_ovrc_table(ncrop,nstage) !< fraction of stem turnover  [1/s]
    real (kind=kind_phys) :: rt_ovrc_table(ncrop,nstage) !< fraction of root tunrover  [1/s]
    real (kind=kind_phys) :: lfmr25_table(ncrop)         !<  leaf maintenance respiration at 25c [umol co2/m**2  /s]
    real (kind=kind_phys) :: stmr25_table(ncrop)         !<  stem maintenance respiration at 25c [umol co2/kg bio/s]
    real (kind=kind_phys) :: rtmr25_table(ncrop)         !<  root maintenance respiration at 25c [umol co2/kg bio/s]
    real (kind=kind_phys) :: grainmr25_table(ncrop)      !< grain maintenance respiration at 25c [umol co2/kg bio/s]

    real (kind=kind_phys) :: lfpt_table(ncrop,nstage)    !< fraction of carbohydrate flux to leaf
    real (kind=kind_phys) :: stpt_table(ncrop,nstage)    !< fraction of carbohydrate flux to stem
    real (kind=kind_phys) :: rtpt_table(ncrop,nstage)    !< fraction of carbohydrate flux to root
    real (kind=kind_phys) :: grainpt_table(ncrop,nstage) !< fraction of carbohydrate flux to grain
    real (kind=kind_phys) :: lfct_table(ncrop,nstage)    ! fraction of carbohydrate translocation from leaf to grain ! Zhe Zhang 2020-07-13
    real (kind=kind_phys) :: stct_table(ncrop,nstage)    !  stem to grain
    real (kind=kind_phys) :: rtct_table(ncrop,nstage)    !  root to grain
    real (kind=kind_phys) :: bio2lai_table(ncrop)        !< leaf are per living leaf biomass [m^2/kg]

! tile drainage parameters
    real (kind=kind_phys)    :: tdsmc_fac_table(max_soiltyp)
    real (kind=kind_phys)    :: td_dc_table(max_soiltyp)
    integer :: td_depth_table(max_soiltyp)
    integer :: drain_layer_opt_table
    real (kind=kind_phys)    :: td_dcoef_table(max_soiltyp)
    real (kind=kind_phys)    :: td_d_table(max_soiltyp)
    real (kind=kind_phys)    :: td_adepth_table(max_soiltyp)
    real (kind=kind_phys)    :: td_radi_table(max_soiltyp)
    real (kind=kind_phys)    :: td_spac_table(max_soiltyp)
    real (kind=kind_phys)    :: td_ddrain_table(max_soiltyp)
    real (kind=kind_phys)    :: klat_fac_table(max_soiltyp)

! mptable.tbl optional parameters

 !------------------------------------------------------------------------------
 ! Saxton and Rawls 2006 Pedo-transfer function coefficients
 !------------------------------------------------------------------------------

    real (kind=kind_phys) ::  sr2006_theta_1500t_a    !< sand coefficient
    real (kind=kind_phys) ::  sr2006_theta_1500t_b    !< clay coefficient
    real (kind=kind_phys) ::  sr2006_theta_1500t_c    !< orgm coefficient
    real (kind=kind_phys) ::  sr2006_theta_1500t_d    !< sand*orgm coefficient
    real (kind=kind_phys) ::  sr2006_theta_1500t_e    !< clay*orgm coefficient
    real (kind=kind_phys) ::  sr2006_theta_1500t_f    !< sand*clay coefficient
    real (kind=kind_phys) ::  sr2006_theta_1500t_g    !< constant adjustment

    real (kind=kind_phys) ::  sr2006_theta_1500_a     !< theta_1500t coefficient
    real (kind=kind_phys) ::  sr2006_theta_1500_b     !< constant adjustment

    real (kind=kind_phys) ::  sr2006_theta_33t_a      !< sand coefficient
    real (kind=kind_phys) ::  sr2006_theta_33t_b      !< clay coefficient
    real (kind=kind_phys) ::  sr2006_theta_33t_c      !< orgm coefficient
    real (kind=kind_phys) ::  sr2006_theta_33t_d      !< sand*orgm coefficient
    real (kind=kind_phys) ::  sr2006_theta_33t_e      !< clay*orgm coefficient
    real (kind=kind_phys) ::  sr2006_theta_33t_f      !< sand*clay coefficient
    real (kind=kind_phys) ::  sr2006_theta_33t_g      !< constant adjustment

    real (kind=kind_phys) ::  sr2006_theta_33_a       !< theta_33t*theta_33t coefficient
    real (kind=kind_phys) ::  sr2006_theta_33_b       !< theta_33t coefficient
    real (kind=kind_phys) ::  sr2006_theta_33_c       !< constant adjustment

    real (kind=kind_phys) ::  sr2006_theta_s33t_a     !< sand coefficient
    real (kind=kind_phys) ::  sr2006_theta_s33t_b     !< clay coefficient
    real (kind=kind_phys) ::  sr2006_theta_s33t_c     !< orgm coefficient
    real (kind=kind_phys) ::  sr2006_theta_s33t_d     !< sand*orgm coefficient
    real (kind=kind_phys) ::  sr2006_theta_s33t_e     !< clay*orgm coefficient
    real (kind=kind_phys) ::  sr2006_theta_s33t_f     !< sand*clay coefficient
    real (kind=kind_phys) ::  sr2006_theta_s33t_g     !< constant adjustment

    real (kind=kind_phys) ::  sr2006_theta_s33_a      !< theta_s33t coefficient
    real (kind=kind_phys) ::  sr2006_theta_s33_b      !< constant adjustment

    real (kind=kind_phys) ::  sr2006_psi_et_a         !< sand coefficient
    real (kind=kind_phys) ::  sr2006_psi_et_b         !< clay coefficient
    real (kind=kind_phys) ::  sr2006_psi_et_c         !< theta_s33 coefficient
    real (kind=kind_phys) ::  sr2006_psi_et_d         !< sand*theta_s33 coefficient
    real (kind=kind_phys) ::  sr2006_psi_et_e         !< clay*theta_s33 coefficient
    real (kind=kind_phys) ::  sr2006_psi_et_f         !< sand*clay coefficient
    real (kind=kind_phys) ::  sr2006_psi_et_g         !< constant adjustment

    real (kind=kind_phys) ::  sr2006_psi_e_a          !< psi_et*psi_et coefficient
    real (kind=kind_phys) ::  sr2006_psi_e_b          !< psi_et coefficient
    real (kind=kind_phys) ::  sr2006_psi_e_c          !< constant adjustment

    real (kind=kind_phys) ::  sr2006_smcmax_a         !< sand adjustment
    real (kind=kind_phys) ::  sr2006_smcmax_b         !< constant adjustment

contains

  subroutine read_mp_table_parameters(errmsg, errflg)
    implicit none

    character(len=*),     intent(out) :: errmsg
    integer,              intent(out) :: errflg

    ! vegetation parameters
    character(len=256)                     :: dataset_identifier
    character(len=256)                     :: veg_dataset_description
    logical                                :: file_named
    integer                                :: ierr, ik, im
    integer                                :: nveg, isurban, iswater, isbarren, isice, iscrop, eblforest, natural
    integer                                :: lcz_1, lcz_2, lcz_3, lcz_4, lcz_5, lcz_6, lcz_7, lcz_8, lcz_9, lcz_10, lcz_11
    real (kind=kind_phys), dimension(mvt) ::  sai_jan, sai_feb, sai_mar, sai_apr, sai_may, sai_jun, sai_jul, sai_aug,        &
                                              sai_sep, sai_oct, sai_nov, sai_dec, lai_jan, lai_feb, lai_mar, lai_apr,        &
                                              lai_may, lai_jun, lai_jul, lai_aug, lai_sep, lai_oct, lai_nov, lai_dec,        &
                                              rhol_vis, rhol_nir, rhos_vis, rhos_nir, taul_vis, taul_nir, taus_vis, taus_nir,&
                                              ch2op, dleaf, z0mvt, hvt, hvb, z0mhvt,                                         &
                                              den, rc, mfsno, scffac, cbiom, xl, cwpvt, c3psn, kc25,                         &
                                              akc, ko25, ako, avcmx, aqe, ltovrc, dilefc, dilefw, rmf25, sla, fragr, tmin,   &
                                              vcmx25, tdlef, bp, mp, qe25, rms25, rmr25, arm, folnmx, wdpool, wrrat, mrp,    &
                                              nroot, rgl, rs, hs, topt, rsmax, rtovrc, rswoodc, bf, wstrc, laimin,           &
                                              xsamin, eps1, eps2, eps3, eps4, eps5
    namelist / noahmp_usgs_veg_categories /   veg_dataset_description, nveg
    namelist / noahmp_usgs_parameters     /   isurban, iswater, isbarren, isice, iscrop, eblforest, natural,                 &
                                              lcz_1, lcz_2, lcz_3, lcz_4, lcz_5, lcz_6, lcz_7, lcz_8, lcz_9, lcz_10, lcz_11, &
                                              ch2op, dleaf, z0mvt, hvt, hvb, z0mhvt,                                         &
                                              den, rc, mfsno, scffac, cbiom, xl, cwpvt, c3psn, kc25,                         &
                                              akc, ko25, ako, avcmx, aqe, ltovrc, dilefc, dilefw, rmf25, sla, fragr, tmin,   &
                                              vcmx25, tdlef, bp, mp, qe25, rms25, rmr25, arm, folnmx, wdpool, wrrat, mrp,    &
                                              nroot, rgl, rs, hs, topt, rsmax, rtovrc, rswoodc, bf, wstrc, laimin,           &
                                              xsamin, sai_jan, sai_feb, sai_mar, sai_apr, sai_may,                           &
                                              sai_jun, sai_jul, sai_aug, sai_sep, sai_oct, sai_nov, sai_dec, lai_jan,        &
                                              lai_feb, lai_mar, lai_apr, lai_may, lai_jun, lai_jul, lai_aug, lai_sep,        &
                                              lai_oct, lai_nov, lai_dec, rhol_vis, rhol_nir, rhos_vis, rhos_nir, taul_vis,   &
                                              taul_nir, taus_vis, taus_nir, eps1, eps2, eps3, eps4, eps5
    namelist / noahmp_modis_veg_categories /  veg_dataset_description, nveg
    namelist / noahmp_modis_parameters     /  isurban, iswater, isbarren, isice, iscrop, eblforest, natural,                 &
                                              lcz_1, lcz_2, lcz_3, lcz_4, lcz_5, lcz_6, lcz_7, lcz_8, lcz_9, lcz_10, lcz_11, &
                                              ch2op, dleaf, z0mvt, hvt, hvb, z0mhvt,                                         &
                                              den, rc, mfsno, scffac, cbiom, xl, cwpvt, c3psn, kc25,                         &
                                              akc, ko25, ako, avcmx, aqe, ltovrc, dilefc, dilefw, rmf25, sla, fragr, tmin,   &
                                              vcmx25, tdlef, bp, mp, qe25, rms25, rmr25, arm, folnmx, wdpool, wrrat, mrp,    &
                                              nroot, rgl, rs, hs, topt, rsmax, rtovrc, rswoodc, bf, wstrc, laimin,           &
                                              xsamin, sai_jan, sai_feb, sai_mar, sai_apr, sai_may,                           &
                                              sai_jun, sai_jul, sai_aug, sai_sep, sai_oct, sai_nov, sai_dec, lai_jan,        &
                                              lai_feb, lai_mar, lai_apr, lai_may, lai_jun, lai_jul, lai_aug, lai_sep,        &
                                              lai_oct, lai_nov, lai_dec, rhol_vis, rhol_nir, rhos_vis, rhos_nir, taul_vis,   &
                                              taul_nir, taus_vis, taus_nir, eps1, eps2, eps3, eps4, eps5
    ! soil parameters
    character(len=256)                             :: message
    character(len=10)                              :: sltype
    integer                                        :: slcats
    real (kind=kind_phys), dimension(max_soiltyp)  :: bb, drysmc, maxsmc, refsmc, satpsi, satdk, satdw, wltsmc, qtz,    &
                                                      bvic, axaj, bxaj, xxaj, bdvic, bbvic, gdvic, hc
    namelist / noahmp_stas_soil_categories /          sltype, slcats
    namelist / noahmp_soil_stas_parameters /          bb, drysmc, maxsmc, refsmc, satpsi, satdk, satdw, wltsmc, qtz,    &
                                                      bvic, axaj, bxaj, xxaj, bdvic, bbvic, gdvic
    namelist / noahmp_soil_stas_ruc_parameters /      bb, drysmc, hc, maxsmc, refsmc, satpsi, satdk, satdw, wltsmc, qtz,    &
                                                      bvic, axaj, bxaj, xxaj, bdvic, bbvic, gdvic

    ! general parameters
    real (kind=kind_phys)                       :: csoil_data, refdk_data, refkdt_data, frzk_data, zbot_data, czil_data
    real (kind=kind_phys), dimension(num_slope) :: slope_data
    namelist / noahmp_general_parameters /          slope_data, csoil_data, refdk_data, refkdt_data, frzk_data, zbot_data,   &
                                                    czil_data

    ! radiation parameters
    real (kind=kind_phys)                   :: betads, betais, eice
    real (kind=kind_phys), dimension(mband) :: albice, alblak, omegas 
    real (kind=kind_phys), dimension(2)     :: eg
    real (kind=kind_phys), dimension(msc)   :: albsat_vis, albsat_nir, albdry_vis, albdry_nir
    namelist / noahmp_rad_parameters /          albsat_vis, albsat_nir, albdry_vis, albdry_nir, albice, alblak, omegas,      &
                                                betads, betais, eg, eice

    ! global parameters
    real (kind=kind_phys)                    :: co2, o2, timean, fsatmx, z0sno, ssi, snow_ret_fac ,snow_emis, swemx, tau0,   &
                                                grain_growth, extra_growth, dirt_soot, bats_cosz, bats_vis_new,              &
                                                bats_nir_new, bats_vis_age, bats_nir_age, bats_vis_dir, bats_nir_dir,        &
                                                rsurf_snow, rsurf_exp, c2_snowcompact, c3_snowcompact, c4_snowcompact,       &
                                                c5_snowcompact, dm_snowcompact, eta0_snowcompact, snliqmaxfrac, swemaxgla,   &
                                                wslmax, rous, cmic, snowden_max, class_alb_ref, class_sno_age, class_alb_new,&
                                                psiwlt, z0soil, z0lake
    namelist / noahmp_global_parameters /       co2, o2, timean, fsatmx, z0sno, ssi, snow_ret_fac ,snow_emis, swemx, tau0,   &
                                                grain_growth, extra_growth, dirt_soot, bats_cosz, bats_vis_new,              &
                                                bats_nir_new, bats_vis_age, bats_nir_age, bats_vis_dir, bats_nir_dir,        &
                                                rsurf_snow, rsurf_exp, c2_snowcompact, c3_snowcompact, c4_snowcompact,       &
                                                c5_snowcompact, dm_snowcompact, eta0_snowcompact, snliqmaxfrac, swemaxgla,   &
                                                wslmax, rous, cmic, snowden_max, class_alb_ref, class_sno_age, class_alb_new,&
                                                psiwlt, z0soil, z0lake

    ! irrigation parameters
    integer                                  :: irr_har
    real (kind=kind_phys)                    :: irr_frac, irr_lai, irr_mad, filoss, sprir_rate, micir_rate, firtfac, ir_rain
    namelist / noahmp_irrigation_parameters /   irr_frac, irr_har, irr_lai, irr_mad, filoss, sprir_rate, micir_rate, firtfac,&
                                                ir_rain

    ! crop parameters
    integer                                  :: default_crop
    integer               , dimension(ncrop) :: pltday, hsday
    real (kind=kind_phys), dimension(ncrop)  :: plantpop, irri, gddtbase, gddtcut, gdds1, gdds2, gdds3, gdds4, gdds5, c3psni,&
                                                kc25i, akci, ko25i, akoi, avcmxi, vcmx25i, bpi, mpi, folnmxi, qe25i, aref,   &
                                                psnrf, i2par, tassim0, tassim1, tassim2, k, epsi, q10mr, lefreez,            &
                                                dile_fc_s1, dile_fc_s2, dile_fc_s3, dile_fc_s4, dile_fc_s5, dile_fc_s6,      &
                                                dile_fc_s7, dile_fc_s8, dile_fw_s1, dile_fw_s2, dile_fw_s3, dile_fw_s4,      &
                                                dile_fw_s5, dile_fw_s6, dile_fw_s7, dile_fw_s8, fra_gr, lf_ovrc_s1,          &
                                                lf_ovrc_s2, lf_ovrc_s3, lf_ovrc_s4, lf_ovrc_s5, lf_ovrc_s6, lf_ovrc_s7,      &
                                                lf_ovrc_s8, st_ovrc_s1, st_ovrc_s2, st_ovrc_s3, st_ovrc_s4, st_ovrc_s5,      &
                                                st_ovrc_s6, st_ovrc_s7, st_ovrc_s8, rt_ovrc_s1, rt_ovrc_s2, rt_ovrc_s3,      &
                                                rt_ovrc_s4, rt_ovrc_s5, rt_ovrc_s6, rt_ovrc_s7, rt_ovrc_s8, lfmr25, stmr25,  &
                                                rtmr25, grainmr25, lfpt_s1, lfpt_s2, lfpt_s3, lfpt_s4, lfpt_s5, lfpt_s6,     &
                                                lfpt_s7, lfpt_s8, stpt_s1, stpt_s2, stpt_s3, stpt_s4, stpt_s5, stpt_s6,      &
                                                stpt_s7, stpt_s8, rtpt_s1, rtpt_s2, rtpt_s3, rtpt_s4, rtpt_s5, rtpt_s6,      &
                                                rtpt_s7, rtpt_s8, grainpt_s1, grainpt_s2, grainpt_s3, grainpt_s4, grainpt_s5,&
                                                grainpt_s6, grainpt_s7, grainpt_s8, lfct_s1, lfct_s2, lfct_s3, lfct_s4,      &
                                                lfct_s5, lfct_s6, lfct_s7, lfct_s8, stct_s1, stct_s2, stct_s3, stct_s4,      &
                                                stct_s5, stct_s6, stct_s7, stct_s8, rtct_s1, rtct_s2, rtct_s3, rtct_s4,      &
                                                rtct_s5, rtct_s6, rtct_s7, rtct_s8, bio2lai
    namelist / noahmp_crop_parameters /         default_crop, pltday, hsday, plantpop, irri, gddtbase, gddtcut, gdds1, gdds2,&
                                                gdds3, gdds4, gdds5, c3psni, kc25i, akci, ko25i, akoi, avcmxi, vcmx25i, bpi, &
                                                mpi, folnmxi, qe25i, aref, psnrf, i2par, tassim0, tassim1, tassim2, k,       &
                                                epsi,q10mr, lefreez, dile_fc_s1, dile_fc_s2, dile_fc_s3, dile_fc_s4,         &
                                                dile_fc_s5, dile_fc_s6, dile_fc_s7, dile_fc_s8, dile_fw_s1, dile_fw_s2,      &
                                                dile_fw_s3, dile_fw_s4, dile_fw_s5, dile_fw_s6, dile_fw_s7, dile_fw_s8,      &
                                                fra_gr, lf_ovrc_s1, lf_ovrc_s2, lf_ovrc_s3, lf_ovrc_s4, lf_ovrc_s5,          &
                                                lf_ovrc_s6, lf_ovrc_s7, lf_ovrc_s8, st_ovrc_s1, st_ovrc_s2, st_ovrc_s3,      &
                                                st_ovrc_s4, st_ovrc_s5, st_ovrc_s6, st_ovrc_s7, st_ovrc_s8, rt_ovrc_s1,      &
                                                rt_ovrc_s2, rt_ovrc_s3, rt_ovrc_s4, rt_ovrc_s5, rt_ovrc_s6, rt_ovrc_s7,      &
                                                rt_ovrc_s8, lfmr25, stmr25, rtmr25, grainmr25, lfpt_s1, lfpt_s2, lfpt_s3,    &
                                                lfpt_s4, lfpt_s5, lfpt_s6, lfpt_s7, lfpt_s8, stpt_s1, stpt_s2, stpt_s3,      &
                                                stpt_s4, stpt_s5, stpt_s6, stpt_s7, stpt_s8, rtpt_s1, rtpt_s2, rtpt_s3,      &
                                                rtpt_s4, rtpt_s5, rtpt_s6, rtpt_s7, rtpt_s8, grainpt_s1, grainpt_s2,         &
                                                grainpt_s3, grainpt_s4, grainpt_s5, grainpt_s6, grainpt_s7, grainpt_s8,      &
                                                lfct_s1, lfct_s2, lfct_s3, lfct_s4, lfct_s5, lfct_s6, lfct_s7, lfct_s8,      &
                                                stct_s1, stct_s2, stct_s3, stct_s4, stct_s5, stct_s6, stct_s7, stct_s8,      &
                                                rtct_s1, rtct_s2, rtct_s3, rtct_s4, rtct_s5, rtct_s6, rtct_s7, rtct_s8,      &
                                                bio2lai

    ! tile drainage parameters
    integer                                        :: nsoiltype, drain_layer_opt
    integer               , dimension(max_soiltyp) :: td_depth
    real (kind=kind_phys), dimension(max_soiltyp)                   :: tdsmc_fac, td_dc, td_dcoef, td_d, td_adepth, td_radi, td_spac,         &
                                                      td_ddrain, klat_fac
    namelist / noahmp_tiledrain_parameters /          nsoiltype, drain_layer_opt, tdsmc_fac, td_depth, td_dc, td_dcoef, td_d,&
                                                      td_adepth, td_radi, td_spac, td_ddrain, klat_fac

    ! optional parameters
    real (kind=kind_phys)                          :: sr2006_theta_1500t_a, sr2006_theta_1500t_b, sr2006_theta_1500t_c,      &
                                                      sr2006_theta_1500t_d, sr2006_theta_1500t_e, sr2006_theta_1500t_f,      &
                                                      sr2006_theta_1500t_g, sr2006_theta_1500_a , sr2006_theta_1500_b,       &
                                                      sr2006_theta_33t_a, sr2006_theta_33t_b, sr2006_theta_33t_c,            &
                                                      sr2006_theta_33t_d, sr2006_theta_33t_e, sr2006_theta_33t_f,            &
                                                      sr2006_theta_33t_g, sr2006_theta_33_a, sr2006_theta_33_b,              &
                                                      sr2006_theta_33_c, sr2006_theta_s33t_a, sr2006_theta_s33t_b,           &
                                                      sr2006_theta_s33t_c, sr2006_theta_s33t_d, sr2006_theta_s33t_e,         &
                                                      sr2006_theta_s33t_f, sr2006_theta_s33t_g, sr2006_theta_s33_a,          &
                                                      sr2006_theta_s33_b, sr2006_psi_et_a, sr2006_psi_et_b, sr2006_psi_et_c, &
                                                      sr2006_psi_et_d, sr2006_psi_et_e, sr2006_psi_et_f, sr2006_psi_et_g,    &
                                                      sr2006_psi_e_a, sr2006_psi_e_b, sr2006_psi_e_c, sr2006_smcmax_a,       &
                                                      sr2006_smcmax_b
    namelist / noahmp_optional_parameters /           sr2006_theta_1500t_a, sr2006_theta_1500t_b, sr2006_theta_1500t_c,      &
                                                      sr2006_theta_1500t_d, sr2006_theta_1500t_e, sr2006_theta_1500t_f,      &
                                                      sr2006_theta_1500t_g, sr2006_theta_1500_a, sr2006_theta_1500_b,        &
                                                      sr2006_theta_33t_a, sr2006_theta_33t_b, sr2006_theta_33t_c,            &
                                                      sr2006_theta_33t_d, sr2006_theta_33t_e, sr2006_theta_33t_f,            &
                                                      sr2006_theta_33t_g, sr2006_theta_33_a, sr2006_theta_33_b,              &
                                                      sr2006_theta_33_c, sr2006_theta_s33t_a, sr2006_theta_s33t_b,           &
                                                      sr2006_theta_s33t_c, sr2006_theta_s33t_d, sr2006_theta_s33t_e,         &
                                                      sr2006_theta_s33t_f, sr2006_theta_s33t_g, sr2006_theta_s33_a,          &
                                                      sr2006_theta_s33_b, sr2006_psi_et_a, sr2006_psi_et_b, sr2006_psi_et_c, &
                                                      sr2006_psi_et_d, sr2006_psi_et_e, sr2006_psi_et_f, sr2006_psi_et_g,    &
                                                      sr2006_psi_e_a, sr2006_psi_e_b, sr2006_psi_e_c, sr2006_smcmax_a,       &
                                                      sr2006_smcmax_b

    errmsg = ''
    errflg = 0

    ! initialize our variables to bad values, so that if the namelist read fails, we come to a screeching halt as soon as we try to use anything.
    ! vegetation parameters
    isurban_table      = -99999
    iswater_table      = -99999
    isbarren_table     = -99999
    isice_table        = -99999
    iscrop_table       = -99999
    eblforest_table    = -99999
    natural_table      = -99999
    lcz_1_table   = -99999
    lcz_2_table   = -99999
    lcz_3_table   = -99999
    lcz_4_table   = -99999
    lcz_5_table   = -99999
    lcz_6_table   = -99999
    lcz_7_table   = -99999
    lcz_8_table   = -99999
    lcz_9_table   = -99999
    lcz_10_table   = -99999
    lcz_11_table   = -99999
    ch2op_table  = -1.0e36
    dleaf_table  = -1.0e36
    z0mvt_table  = -1.0e36
    hvt_table    = -1.0e36
    hvb_table    = -1.0e36
    z0mhvt_table = -1.0e36
    den_table    = -1.0e36
    rc_table     = -1.0e36
    mfsno_table  = -1.0e36
    scffac_table = -1.0e36
    cbiom_table  = -1.0e36
    rhol_table   = -1.0e36
    rhos_table   = -1.0e36
    taul_table   = -1.0e36
    taus_table   = -1.0e36
    xl_table     = -1.0e36
    cwpvt_table  = -1.0e36
    c3psn_table  = -1.0e36
    kc25_table   = -1.0e36
    akc_table    = -1.0e36
    ko25_table   = -1.0e36
    ako_table    = -1.0e36
    avcmx_table  = -1.0e36
    aqe_table    = -1.0e36
    ltovrc_table = -1.0e36
    dilefc_table = -1.0e36
    dilefw_table = -1.0e36
    rmf25_table  = -1.0e36
    sla_table    = -1.0e36
    fragr_table  = -1.0e36
    tmin_table   = -1.0e36
    vcmx25_table = -1.0e36
    tdlef_table  = -1.0e36
    bp_table     = -1.0e36
    mp_table     = -1.0e36
    qe25_table   = -1.0e36
    rms25_table  = -1.0e36
    rmr25_table  = -1.0e36
    arm_table    = -1.0e36
    folnmx_table = -1.0e36
    wdpool_table = -1.0e36
    wrrat_table  = -1.0e36
    mrp_table    = -1.0e36
    saim_table   = -1.0e36
    laim_table   = -1.0e36
    nroot_table  = -1.0e36
    rgl_table    = -1.0e36
    rs_table     = -1.0e36
    hs_table     = -1.0e36
    topt_table   = -1.0e36
    rsmax_table  = -1.0e36
    ! not used in the current ufs version
!   rtovrc_table       = -1.0e36
!   rswoodc_table      = -1.0e36
!   bf_table           = -1.0e36
!   wstrc_table        = -1.0e36
!   laimin_table       = -1.0e36
!   xsamin_table       = -1.0e36

    ! soil parameters

       bexp_table = -1.0e36
     smcdry_table = -1.0e36
         f1_table = -1.0e36
     smcmax_table = -1.0e36
     smcref_table = -1.0e36
     psisat_table = -1.0e36
      dksat_table = -1.0e36
      dwsat_table = -1.0e36
     smcwlt_table = -1.0e36
     quartz_table = -1.0e36
       bvic_table = -1.0e36
       axaj_table = -1.0e36
       bxaj_table = -1.0e36
       xxaj_table = -1.0e36
      bdvic_table = -1.0e36
      gdvic_table = -1.0e36
      bbvic_table = -1.0e36

    ! general parameters
      slope_table = -1.0e36
      csoil_table = -1.0e36
      refdk_table = -1.0e36
     refkdt_table = -1.0e36
       frzk_table = -1.0e36
       zbot_table = -1.0e36
       czil_table = -1.0e36

    ! radiation parameters
    albsat_table     = -1.0e36
    albdry_table     = -1.0e36
    albice_table     = -1.0e36
    alblak_table     = -1.0e36
    omegas_table     = -1.0e36
    betads_table     = -1.0e36
    betais_table     = -1.0e36
    eg_table         = -1.0e36
!   eice_table       = -1.0e36

    ! global parameters
       co2_table     = -1.0e36
        o2_table     = -1.0e36
    timean_table     = -1.0e36
    fsatmx_table     = -1.0e36
     z0sno_table     = -1.0e36
       ssi_table     = -1.0e36
snow_ret_fac_table   = -1.0e36
   snow_emis_table   = -1.0e36
       swemx_table   = -1.0e36
        tau0_table   = -1.0e36
grain_growth_table   = -1.0e36
extra_growth_table   = -1.0e36
   dirt_soot_table   = -1.0e36
   bats_cosz_table   = -1.0e36
bats_vis_new_table   = -1.0e36
bats_nir_new_table   = -1.0e36
bats_vis_age_table   = -1.0e36
bats_nir_age_table   = -1.0e36
bats_vis_dir_table   = -1.0e36
bats_nir_dir_table   = -1.0e36
rsurf_snow_table     = -1.0e36
 rsurf_exp_table     = -1.0e36

!   c2_snowcompact_table   = -1.0e36
!   c3_snowcompact_table   = -1.0e36
!   c4_snowcompact_table   = -1.0e36
!   c5_snowcompact_table   = -1.0e36
!   dm_snowcompact_table   = -1.0e36
!   eta0_snowcompact_table = -1.0e36
!   snliqmaxfrac_table     = -1.0e36
!   swemaxgla_table        = -1.0e36
!   wslmax_table           = -1.0e36
!   rous_table             = -1.0e36
!   cmic_table             = -1.0e36
!   snowden_max_table      = -1.0e36
!   class_alb_ref_table    = -1.0e36
!   class_sno_age_table    = -1.0e36
!   class_alb_new_table    = -1.0e36
!   psiwlt_table           = -1.0e36
!   z0soil_table           = -1.0e36
!   z0lake_table           = -1.0e36

    ! irrigation parameters
    irr_har_table    =  -99999        ! number of days before harvest date to stop irrigation 
    irr_frac_table   = -1.0e36    ! irrigation fraction
    irr_lai_table    = -1.0e36    ! minimum lai to trigger irrigation
    irr_mad_table    = -1.0e36    ! management allowable deficit (0-1)
    filoss_table     = -1.0e36    ! fraction of flood irrigation loss (0-1) 
    sprir_rate_table = -1.0e36    ! mm/h, sprinkler irrigation rate
    micir_rate_table = -1.0e36    ! mm/h, micro irrigation rate
    firtfac_table    = -1.0e36    ! flood application rate factor
    ir_rain_table    = -1.0e36    ! maximum precipitation to stop irrigation trigger

    ! crop parameters
 default_crop_table     = -99999
       pltday_table     = -99999
        hsday_table     = -99999
     plantpop_table     = -1.0e36
         irri_table     = -1.0e36
     gddtbase_table     = -1.0e36
      gddtcut_table     = -1.0e36
        gdds1_table     = -1.0e36
        gdds2_table     = -1.0e36
        gdds3_table     = -1.0e36
        gdds4_table     = -1.0e36
        gdds5_table     = -1.0e36
       c3psni_table     = -1.0e36 ! parameter from psn copied from stomata ! zhe zhang 2020-07-13
        kc25i_table     = -1.0e36
         akci_table     = -1.0e36
        ko25i_table     = -1.0e36
         akoi_table     = -1.0e36
       avcmxi_table     = -1.0e36
      vcmx25i_table     = -1.0e36
          bpi_table     = -1.0e36
          mpi_table     = -1.0e36
      folnmxi_table     = -1.0e36
        qe25i_table     = -1.0e36 ! ends here
!???         c3c4_table     = -99999
         aref_table     = -1.0e36
        psnrf_table     = -1.0e36
        i2par_table     = -1.0e36
      tassim0_table     = -1.0e36
      tassim1_table     = -1.0e36
      tassim2_table     = -1.0e36
            k_table     = -1.0e36
         epsi_table     = -1.0e36
        q10mr_table     = -1.0e36
      foln_mx_table     = -1.0e36
      lefreez_table     = -1.0e36
      dile_fc_table     = -1.0e36
      dile_fw_table     = -1.0e36
       fra_gr_table     = -1.0e36
      lf_ovrc_table     = -1.0e36
      st_ovrc_table     = -1.0e36
      rt_ovrc_table     = -1.0e36
       lfmr25_table     = -1.0e36
       stmr25_table     = -1.0e36
       rtmr25_table     = -1.0e36
    grainmr25_table     = -1.0e36
         lfpt_table     = -1.0e36
         stpt_table     = -1.0e36
         rtpt_table     = -1.0e36
      grainpt_table     = -1.0e36
         lfct_table     = -1.0e36 ! convert start
         stct_table     = -1.0e36
         rtct_table     = -1.0e36 ! convert end
      bio2lai_table     = -1.0e36

    ! tile drainage parameters

    drain_layer_opt_table  = -99999
    td_depth_table         = -99999
    tdsmc_fac_table        = -1.0e36
    td_dc_table            = -1.0e36
    td_dcoef_table         = -1.0e36
    td_d_table             = -1.0e36
    td_adepth_table        = -1.0e36
    td_radi_table          = -1.0e36
    td_spac_table          = -1.0e36
    td_ddrain_table        = -1.0e36
    klat_fac_table         = -1.0e36

    ! optional parameters
!   sr2006_theta_1500t_a_table = -1.0e36
!   sr2006_theta_1500t_b_table = -1.0e36
!   sr2006_theta_1500t_c_table = -1.0e36
!   sr2006_theta_1500t_d_table = -1.0e36
!   sr2006_theta_1500t_e_table = -1.0e36
!   sr2006_theta_1500t_f_table = -1.0e36
!   sr2006_theta_1500t_g_table = -1.0e36
!   sr2006_theta_1500_a_table  = -1.0e36
!   sr2006_theta_1500_b_table  = -1.0e36
!   sr2006_theta_33t_a_table   = -1.0e36
!   sr2006_theta_33t_b_table   = -1.0e36
!   sr2006_theta_33t_c_table   = -1.0e36
!   sr2006_theta_33t_d_table   = -1.0e36
!   sr2006_theta_33t_e_table   = -1.0e36
!   sr2006_theta_33t_f_table   = -1.0e36
!   sr2006_theta_33t_g_table   = -1.0e36
!   sr2006_theta_33_a_table    = -1.0e36
!   sr2006_theta_33_b_table    = -1.0e36
!   sr2006_theta_33_c_table    = -1.0e36
!   sr2006_theta_s33t_a_table  = -1.0e36
!   sr2006_theta_s33t_b_table  = -1.0e36
!   sr2006_theta_s33t_c_table  = -1.0e36
!   sr2006_theta_s33t_d_table  = -1.0e36
!   sr2006_theta_s33t_e_table  = -1.0e36
!   sr2006_theta_s33t_f_table  = -1.0e36
!   sr2006_theta_s33t_g_table  = -1.0e36
!   sr2006_theta_s33_a_table   = -1.0e36
!   sr2006_theta_s33_b_table   = -1.0e36
!   sr2006_psi_et_a_table      = -1.0e36
!   sr2006_psi_et_b_table      = -1.0e36
!   sr2006_psi_et_c_table      = -1.0e36
!   sr2006_psi_et_d_table      = -1.0e36
!   sr2006_psi_et_e_table      = -1.0e36
!   sr2006_psi_et_f_table      = -1.0e36
!   sr2006_psi_et_g_table      = -1.0e36
!   sr2006_psi_e_a_table       = -1.0e36
!   sr2006_psi_e_b_table       = -1.0e36
!   sr2006_psi_e_c_table       = -1.0e36
!   sr2006_smcmax_a_table      = -1.0e36
!   sr2006_smcmax_b_table      = -1.0e36

    !---------------------------------------------------------------
    ! transfer values from table to input variables
    !---------------------------------------------------------------

    !---------------- noahmptable.tbl vegetation parameters

    dataset_identifier = "modified_igbp_modis_noah"

    inquire( file='noahmptable.tbl', exist=file_named )
    if ( file_named ) then
       open(15, file="noahmptable.tbl", status='old', form='formatted', action='read', iostat=ierr)
    else
       open(15, status='old', form='formatted', action='read', iostat=ierr)
    end if
    if ( ierr /= 0 ) then
       errmsg = 'warning: cannot find file noahmptable.tbl'
       errflg = 1
       return
!      write(*,'("warning: cannot find file noahmptable.tbl")')
    endif

    if ( trim(dataset_identifier) == "usgs" ) then
       read(15, noahmp_usgs_veg_categories)
       read(15, noahmp_usgs_parameters)
    elseif ( trim(dataset_identifier) == "modified_igbp_modis_noah" ) then
       read(15,noahmp_modis_veg_categories)
       read(15,noahmp_modis_parameters)
    else
       write(*,'("warning: unrecognized dataset_identifier in subroutine readnoahmptable")')
       write(*,'("warning: dataset_identifier = ''", a, "''")') trim(dataset_identifier)
    endif
    close(15)


    ! assign values
    isurban_table         = isurban
    iswater_table         = iswater
    isbarren_table        = isbarren
    isice_table           = isice
    iscrop_table          = iscrop
    eblforest_table       = eblforest
    natural_table         = natural
    lcz_1_table           = lcz_1
    lcz_2_table           = lcz_2
    lcz_3_table           = lcz_3
    lcz_4_table           = lcz_4
    lcz_5_table           = lcz_5
    lcz_6_table           = lcz_6
    lcz_7_table           = lcz_7
    lcz_8_table           = lcz_8
    lcz_9_table           = lcz_9
    lcz_10_table          = lcz_10
    lcz_11_table          = lcz_11
    ch2op_table  (1:nveg) = ch2op  (1:nveg)
    dleaf_table  (1:nveg) = dleaf  (1:nveg)
    z0mvt_table  (1:nveg) = z0mvt  (1:nveg)
    hvt_table    (1:nveg) = hvt    (1:nveg)
    hvb_table    (1:nveg) = hvb    (1:nveg)
    z0mhvt_table (1:nveg) = z0mhvt (1:nveg)
    den_table    (1:nveg) = den    (1:nveg)
    rc_table     (1:nveg) = rc     (1:nveg)
    mfsno_table  (1:nveg) = mfsno  (1:nveg)
    scffac_table (1:nveg) = scffac (1:nveg)
    cbiom_table  (1:nveg) = cbiom  (1:nveg)
    xl_table     (1:nveg) = xl     (1:nveg)
    cwpvt_table  (1:nveg) = cwpvt  (1:nveg)
    c3psn_table  (1:nveg) = c3psn  (1:nveg)
    kc25_table   (1:nveg) = kc25   (1:nveg)
    akc_table    (1:nveg) = akc    (1:nveg)
    ko25_table   (1:nveg) = ko25   (1:nveg)
    ako_table    (1:nveg) = ako    (1:nveg)
    avcmx_table  (1:nveg) = avcmx  (1:nveg)
    aqe_table    (1:nveg) = aqe    (1:nveg)
    ltovrc_table (1:nveg) = ltovrc (1:nveg)
    dilefc_table (1:nveg) = dilefc (1:nveg)
    dilefw_table (1:nveg) = dilefw (1:nveg)
    rmf25_table  (1:nveg) = rmf25  (1:nveg)
    sla_table    (1:nveg) = sla    (1:nveg)
    fragr_table  (1:nveg) = fragr  (1:nveg)
    tmin_table   (1:nveg) = tmin   (1:nveg)
    vcmx25_table (1:nveg) = vcmx25 (1:nveg)
    tdlef_table  (1:nveg) = tdlef  (1:nveg)
    bp_table     (1:nveg) = bp     (1:nveg)
    mp_table     (1:nveg) = mp     (1:nveg)
    qe25_table   (1:nveg) = qe25   (1:nveg)
    rms25_table  (1:nveg) = rms25  (1:nveg)
    rmr25_table  (1:nveg) = rmr25  (1:nveg)
    arm_table    (1:nveg) = arm    (1:nveg)
    folnmx_table (1:nveg) = folnmx (1:nveg)
    wdpool_table (1:nveg) = wdpool (1:nveg)
    wrrat_table  (1:nveg) = wrrat  (1:nveg)
    mrp_table    (1:nveg) = mrp    (1:nveg)
    nroot_table  (1:nveg) = nroot  (1:nveg)
    rgl_table    (1:nveg) = rgl    (1:nveg)
    rs_table     (1:nveg) = rs     (1:nveg)
    hs_table     (1:nveg) = hs     (1:nveg)
    topt_table   (1:nveg) = topt   (1:nveg)
    rsmax_table  (1:nveg) = rsmax  (1:nveg)
!   rtovrc_table (1:nveg) = rtovrc (1:nveg)
!   rswoodc_table(1:nveg) = rswoodc(1:nveg)
!   bf_table     (1:nveg) = bf     (1:nveg)
!   wstrc_table  (1:nveg) = wstrc  (1:nveg)
!   laimin_table (1:nveg) = laimin (1:nveg)
!   xsamin_table (1:nveg) = xsamin (1:nveg)

    saim_table(1:nveg, 1) = sai_jan(1:nveg)
    saim_table(1:nveg, 2) = sai_feb(1:nveg)
    saim_table(1:nveg, 3) = sai_mar(1:nveg)
    saim_table(1:nveg, 4) = sai_apr(1:nveg)
    saim_table(1:nveg, 5) = sai_may(1:nveg)
    saim_table(1:nveg, 6) = sai_jun(1:nveg)
    saim_table(1:nveg, 7) = sai_jul(1:nveg)
    saim_table(1:nveg, 8) = sai_aug(1:nveg)
    saim_table(1:nveg, 9) = sai_sep(1:nveg)
    saim_table(1:nveg,10) = sai_oct(1:nveg)
    saim_table(1:nveg,11) = sai_nov(1:nveg)
    saim_table(1:nveg,12) = sai_dec(1:nveg)
    laim_table(1:nveg, 1) = lai_jan(1:nveg)
    laim_table(1:nveg, 2) = lai_feb(1:nveg)
    laim_table(1:nveg, 3) = lai_mar(1:nveg)
    laim_table(1:nveg, 4) = lai_apr(1:nveg)
    laim_table(1:nveg, 5) = lai_may(1:nveg)
    laim_table(1:nveg, 6) = lai_jun(1:nveg)
    laim_table(1:nveg, 7) = lai_jul(1:nveg)
    laim_table(1:nveg, 8) = lai_aug(1:nveg)
    laim_table(1:nveg, 9) = lai_sep(1:nveg)
    laim_table(1:nveg,10) = lai_oct(1:nveg)
    laim_table(1:nveg,11) = lai_nov(1:nveg)
    laim_table(1:nveg,12) = lai_dec(1:nveg)
    rhol_table(1:nveg,1)  = rhol_vis(1:nveg) !leaf reflectance: 1=vis, 2=nir
    rhol_table(1:nveg,2)  = rhol_nir(1:nveg) !leaf reflectance: 1=vis, 2=nir
    rhos_table(1:nveg,1)  = rhos_vis(1:nveg) !stem reflectance: 1=vis, 2=nir
    rhos_table(1:nveg,2)  = rhos_nir(1:nveg) !stem reflectance: 1=vis, 2=nir
    taul_table(1:nveg,1)  = taul_vis(1:nveg) !leaf transmittance: 1=vis, 2=nir
    taul_table(1:nveg,2)  = taul_nir(1:nveg) !leaf transmittance: 1=vis, 2=nir
    taus_table(1:nveg,1)  = taus_vis(1:nveg) !stem transmittance: 1=vis, 2=nir
    taus_table(1:nveg,2)  = taus_nir(1:nveg) !stem transmittance: 1=vis, 2=nir

    !---------------- noahmptable.tbl soil parameters
    inquire( file='noahmptable.tbl', exist=file_named )
    if ( file_named ) then
       open(15, file="noahmptable.tbl", status='old', form='formatted', action='read', iostat=ierr)
    else
       open(15, status='old', form='formatted', action='read', iostat=ierr)
    end if
    if ( ierr /= 0 ) then
       errmsg = 'warning: cannot find file noahmptable.tbl'
       errflg = 1
       return
!      write(*,'("warning: cannot find file noahmptable.tbl")')
    endif
    read(15, noahmp_stas_soil_categories)
    if ( trim(sltype) == "stas" ) then
       read(15, noahmp_soil_stas_parameters)
    elseif ( trim(sltype) == "stas_ruc" ) then
       read(15, noahmp_soil_stas_ruc_parameters)
    else
       write(*,'("warning: unrecognized soiltype in subroutine readnoahmptable")')
       write(*,'("warning: dataset_identifier = ''", a, "''")') trim(sltype)
    endif
    close(15)

    ! assign values
!    slcats_table           = slcats
     bexp_table  (1:slcats) = bb    (1:slcats)
     smcdry_table(1:slcats) = drysmc(1:slcats)
     smcmax_table(1:slcats) = maxsmc(1:slcats)
     smcref_table(1:slcats) = refsmc(1:slcats)
     psisat_table(1:slcats) = satpsi(1:slcats)
     dksat_table (1:slcats) = satdk (1:slcats)
     dwsat_table (1:slcats) = satdw (1:slcats)
     smcwlt_table(1:slcats) = wltsmc(1:slcats)
     quartz_table(1:slcats) = qtz   (1:slcats)
     bvic_table  (1:slcats) = bvic  (1:slcats)
     axaj_table  (1:slcats) = axaj  (1:slcats)
     bxaj_table  (1:slcats) = bxaj  (1:slcats)
     xxaj_table  (1:slcats) = xxaj  (1:slcats)
     bdvic_table (1:slcats) = bdvic (1:slcats)
     gdvic_table (1:slcats) = gdvic (1:slcats)
     bbvic_table (1:slcats) = bbvic (1:slcats)

    !---------------- noahmptable.tbl general parameters
    inquire( file='noahmptable.tbl', exist=file_named )
    if ( file_named ) then
       open(15, file="noahmptable.tbl", status='old', form='formatted', action='read', iostat=ierr)
    else
       open(15, status='old', form='formatted', action='read', iostat=ierr)
    end if
    if ( ierr /= 0 ) then
       errmsg = 'warning: cannot find file noahmptable.tbl'
       errflg = 1
       return
!      write(*,'("warning: cannot find file noahmptable.tbl")')
    endif
    read(15, noahmp_general_parameters)
    close(15)

    ! assign values
     slope_table(1:num_slope) = slope_data(1:num_slope)
     csoil_table              = csoil_data
     refdk_table              = refdk_data
     refkdt_table             = refkdt_data
     frzk_table               = frzk_data
     zbot_table               = zbot_data
     czil_table               = czil_data

    !---------------- noahmptable.tbl radiation parameters
    inquire( file='noahmptable.tbl', exist=file_named )
    if ( file_named ) then
      open(15, file="noahmptable.tbl", status='old', form='formatted', action='read', iostat=ierr)
    else
      open(15, status='old', form='formatted', action='read', iostat=ierr)
    end if
    if (ierr /= 0) then
       errmsg = 'warning: cannot find file noahmptable.tbl'
       errflg = 1
       return
!      write(*,'("warning: cannot find file noahmptable.tbl")')
    endif
    read(15,noahmp_rad_parameters)
    close(15)

    ! assign values
     albsat_table(:,1) = albsat_vis ! saturated soil albedos: 1=vis, 2=nir
     albsat_table(:,2) = albsat_nir ! saturated soil albedos: 1=vis, 2=nir
     albdry_table(:,1) = albdry_vis ! dry soil albedos: 1=vis, 2=nir
     albdry_table(:,2) = albdry_nir ! dry soil albedos: 1=vis, 2=nir
     albice_table      = albice
     alblak_table      = alblak
     omegas_table      = omegas
     betads_table      = betads
     betais_table      = betais
     eg_table          = eg
!    eice_table        = eice

    !---------------- noahmptable.tbl global parameters
    inquire( file='noahmptable.tbl', exist=file_named )
    if ( file_named ) then
      open(15, file="noahmptable.tbl", status='old', form='formatted', action='read', iostat=ierr)
    else
      open(15, status='old', form='formatted', action='read', iostat=ierr)
    end if
    if (ierr /= 0) then
       errmsg = 'warning: cannot find file noahmptable.tbl'
       errflg = 1
       return
!      write(*,'("warning: cannot find file noahmptable.tbl")')
    endif
    read(15,noahmp_global_parameters)
    close(15)

    ! assign values
     co2_table              = co2
     o2_table               = o2
     timean_table           = timean
     fsatmx_table           = fsatmx
     z0sno_table            = z0sno
     ssi_table              = ssi
     snow_ret_fac_table     = snow_ret_fac
     snow_emis_table        = snow_emis
     swemx_table            = swemx
     tau0_table             = tau0
     grain_growth_table     = grain_growth
     extra_growth_table     = extra_growth
     dirt_soot_table        = dirt_soot
     bats_cosz_table        = bats_cosz
     bats_vis_new_table     = bats_vis_new
     bats_nir_new_table     = bats_nir_new
     bats_vis_age_table     = bats_vis_age
     bats_nir_age_table     = bats_nir_age
     bats_vis_dir_table     = bats_vis_dir
     bats_nir_dir_table     = bats_nir_dir
     rsurf_snow_table       = rsurf_snow
     rsurf_exp_table        = rsurf_exp
!    c2_snowcompact_table   = c2_snowcompact
!    c3_snowcompact_table   = c3_snowcompact
!    c4_snowcompact_table   = c4_snowcompact
!    c5_snowcompact_table   = c5_snowcompact
!    dm_snowcompact_table   = dm_snowcompact
!    eta0_snowcompact_table = eta0_snowcompact
!    snliqmaxfrac_table     = snliqmaxfrac
!    swemaxgla_table        = swemaxgla
!    wslmax_table           = wslmax
!    rous_table             = rous
!    cmic_table             = cmic
!    snowden_max_table      = snowden_max
!    class_alb_ref_table    = class_alb_ref
!    class_sno_age_table    = class_sno_age
!    class_alb_new_table    = class_alb_new
!    psiwlt_table           = psiwlt
!    z0soil_table           = z0soil
!    z0lake_table           = z0lake

    !---------------- noahmptable.tbl irrigation parameters
    inquire( file='noahmptable.tbl', exist=file_named )
    if ( file_named ) then
      open(15, file="noahmptable.tbl", status='old', form='formatted', action='read', iostat=ierr)
    else
      open(15, status='old', form='formatted', action='read', iostat=ierr)
    end if
    if (ierr /= 0) then
       errmsg = 'warning: cannot find file noahmptable.tbl'
       errflg = 1
       return
!      write(*,'("warning: cannot find file noahmptable.tbl")')
    endif
    read(15,noahmp_irrigation_parameters)
    close(15)

    ! assign values
     irr_frac_table   = irr_frac
     irr_har_table    = irr_har
     irr_lai_table    = irr_lai
     irr_mad_table    = irr_mad
     filoss_table     = filoss  
     sprir_rate_table = sprir_rate
     micir_rate_table = micir_rate
     firtfac_table    = firtfac
     ir_rain_table    = ir_rain 

    !---------------- noahmptable.tbl crop parameters
    inquire( file='noahmptable.tbl', exist=file_named )
    if ( file_named ) then
      open(15, file="noahmptable.tbl", status='old', form='formatted', action='read', iostat=ierr)
    else
      open(15, status='old', form='formatted', action='read', iostat=ierr)
    end if
    if (ierr /= 0) then
       errmsg = 'warning: cannot find file noahmptable.tbl'
       errflg = 1
       return
!      write(*,'("warning: cannot find file noahmptable.tbl")')
    endif
    read(15,noahmp_crop_parameters)
    close(15)

    ! assign values
     default_crop_table     = default_crop
     pltday_table           = pltday
     hsday_table            = hsday
     plantpop_table         = plantpop
     irri_table             = irri
     gddtbase_table         = gddtbase
     gddtcut_table          = gddtcut
     gdds1_table            = gdds1
     gdds2_table            = gdds2
     gdds3_table            = gdds3
     gdds4_table            = gdds4
     gdds5_table            = gdds5
     c3psni_table (1:5)     = c3psni (1:5)
     kc25i_table  (1:5)     = kc25i  (1:5)
     akci_table   (1:5)     = akci   (1:5)
     ko25i_table  (1:5)     = ko25i  (1:5)
     akoi_table   (1:5)     = akoi   (1:5)
     avcmxi_table (1:5)     = avcmxi (1:5)
     vcmx25i_table(1:5)     = vcmx25i(1:5)
     bpi_table    (1:5)     = bpi    (1:5)
     mpi_table    (1:5)     = mpi    (1:5)
     folnmxi_table(1:5)     = folnmxi(1:5)
     qe25i_table  (1:5)     = qe25i  (1:5)
     aref_table             = aref
     psnrf_table            = psnrf
     i2par_table            = i2par
     tassim0_table          = tassim0
     tassim1_table          = tassim1
     tassim2_table          = tassim2
     k_table                = k
     epsi_table             = epsi
     q10mr_table            = q10mr
     lefreez_table          = lefreez
     fra_gr_table           = fra_gr
     lfmr25_table           = lfmr25
     stmr25_table           = stmr25
     rtmr25_table           = rtmr25
     grainmr25_table        = grainmr25
     bio2lai_table          = bio2lai
     dile_fc_table(:,1)     = dile_fc_s1
     dile_fc_table(:,2)     = dile_fc_s2
     dile_fc_table(:,3)     = dile_fc_s3
     dile_fc_table(:,4)     = dile_fc_s4
     dile_fc_table(:,5)     = dile_fc_s5
     dile_fc_table(:,6)     = dile_fc_s6
     dile_fc_table(:,7)     = dile_fc_s7
     dile_fc_table(:,8)     = dile_fc_s8
     dile_fw_table(:,1)     = dile_fw_s1
     dile_fw_table(:,2)     = dile_fw_s2
     dile_fw_table(:,3)     = dile_fw_s3
     dile_fw_table(:,4)     = dile_fw_s4
     dile_fw_table(:,5)     = dile_fw_s5
     dile_fw_table(:,6)     = dile_fw_s6
     dile_fw_table(:,7)     = dile_fw_s7
     dile_fw_table(:,8)     = dile_fw_s8
     lf_ovrc_table(:,1)     = lf_ovrc_s1
     lf_ovrc_table(:,2)     = lf_ovrc_s2
     lf_ovrc_table(:,3)     = lf_ovrc_s3
     lf_ovrc_table(:,4)     = lf_ovrc_s4
     lf_ovrc_table(:,5)     = lf_ovrc_s5
     lf_ovrc_table(:,6)     = lf_ovrc_s6
     lf_ovrc_table(:,7)     = lf_ovrc_s7
     lf_ovrc_table(:,8)     = lf_ovrc_s8
     st_ovrc_table(:,1)     = st_ovrc_s1
     st_ovrc_table(:,2)     = st_ovrc_s2
     st_ovrc_table(:,3)     = st_ovrc_s3
     st_ovrc_table(:,4)     = st_ovrc_s4
     st_ovrc_table(:,5)     = st_ovrc_s5
     st_ovrc_table(:,6)     = st_ovrc_s6
     st_ovrc_table(:,7)     = st_ovrc_s7
     st_ovrc_table(:,8)     = st_ovrc_s8
     rt_ovrc_table(:,1)     = rt_ovrc_s1
     rt_ovrc_table(:,2)     = rt_ovrc_s2
     rt_ovrc_table(:,3)     = rt_ovrc_s3
     rt_ovrc_table(:,4)     = rt_ovrc_s4
     rt_ovrc_table(:,5)     = rt_ovrc_s5
     rt_ovrc_table(:,6)     = rt_ovrc_s6
     rt_ovrc_table(:,7)     = rt_ovrc_s7
     rt_ovrc_table(:,8)     = rt_ovrc_s8
     lfpt_table   (:,1)     = lfpt_s1
     lfpt_table   (:,2)     = lfpt_s2
     lfpt_table   (:,3)     = lfpt_s3
     lfpt_table   (:,4)     = lfpt_s4
     lfpt_table   (:,5)     = lfpt_s5
     lfpt_table   (:,6)     = lfpt_s6
     lfpt_table   (:,7)     = lfpt_s7
     lfpt_table   (:,8)     = lfpt_s8
     stpt_table   (:,1)     = stpt_s1
     stpt_table   (:,2)     = stpt_s2
     stpt_table   (:,3)     = stpt_s3
     stpt_table   (:,4)     = stpt_s4
     stpt_table   (:,5)     = stpt_s5
     stpt_table   (:,6)     = stpt_s6
     stpt_table   (:,7)     = stpt_s7
     stpt_table   (:,8)     = stpt_s8
     rtpt_table   (:,1)     = rtpt_s1
     rtpt_table   (:,2)     = rtpt_s2
     rtpt_table   (:,3)     = rtpt_s3
     rtpt_table   (:,4)     = rtpt_s4
     rtpt_table   (:,5)     = rtpt_s5
     rtpt_table   (:,6)     = rtpt_s6
     rtpt_table   (:,7)     = rtpt_s7
     rtpt_table   (:,8)     = rtpt_s8
     grainpt_table(:,1)     = grainpt_s1
     grainpt_table(:,2)     = grainpt_s2
     grainpt_table(:,3)     = grainpt_s3
     grainpt_table(:,4)     = grainpt_s4
     grainpt_table(:,5)     = grainpt_s5
     grainpt_table(:,6)     = grainpt_s6
     grainpt_table(:,7)     = grainpt_s7
     grainpt_table(:,8)     = grainpt_s8
     lfct_table   (:,1)     = lfct_s1
     lfct_table   (:,2)     = lfct_s2
     lfct_table   (:,3)     = lfct_s3
     lfct_table   (:,4)     = lfct_s4
     lfct_table   (:,5)     = lfct_s5
     lfct_table   (:,6)     = lfct_s6
     lfct_table   (:,7)     = lfct_s7
     lfct_table   (:,8)     = lfct_s8
     stct_table   (:,1)     = stct_s1
     stct_table   (:,2)     = stct_s2
     stct_table   (:,3)     = stct_s3
     stct_table   (:,4)     = stct_s4
     stct_table   (:,5)     = stct_s5
     stct_table   (:,6)     = stct_s6
     stct_table   (:,7)     = stct_s7
     stct_table   (:,8)     = stct_s8
     rtct_table   (:,1)     = rtct_s1
     rtct_table   (:,2)     = rtct_s2
     rtct_table   (:,3)     = rtct_s3
     rtct_table   (:,4)     = rtct_s4
     rtct_table   (:,5)     = rtct_s5
     rtct_table   (:,6)     = rtct_s6
     rtct_table   (:,7)     = rtct_s7
     rtct_table   (:,8)     = rtct_s8

    !---------------- noahmptable.tbl tile drainage parameters
    inquire( file='noahmptable.tbl', exist=file_named )
    if ( file_named ) then
      open(15, file="noahmptable.tbl", status='old', form='formatted', action='read', iostat=ierr)
    else
      open(15, status='old', form='formatted', action='read', iostat=ierr)
    end if
    if (ierr /= 0) then
       errmsg = 'warning: cannot find file noahmptable.tbl'
       errflg = 1
       return
!      write(*,'("warning: cannot find file noahmptable.tbl")')
    endif
    read(15,noahmp_tiledrain_parameters)
    close(15)

    ! assign values
     drain_layer_opt_table        = drain_layer_opt
     tdsmc_fac_table(1:nsoiltype) = tdsmc_fac(1:nsoiltype)
     td_depth_table (1:nsoiltype) = td_depth (1:nsoiltype)
     td_dc_table    (1:nsoiltype) = td_dc    (1:nsoiltype)
     td_dcoef_table (1:nsoiltype) = td_dcoef (1:nsoiltype)
     td_d_table     (1:nsoiltype) = td_d     (1:nsoiltype)
     td_adepth_table(1:nsoiltype) = td_adepth(1:nsoiltype)
     td_radi_table  (1:nsoiltype) = td_radi  (1:nsoiltype)
     td_spac_table  (1:nsoiltype) = td_spac  (1:nsoiltype)
     td_ddrain_table(1:nsoiltype) = td_ddrain(1:nsoiltype)
     klat_fac_table (1:nsoiltype) = klat_fac (1:nsoiltype)

    !---------------- noahmptable.tbl optional parameters
    inquire( file='noahmptable.tbl', exist=file_named )
    if ( file_named ) then
      open(15, file="noahmptable.tbl", status='old', form='formatted', action='read', iostat=ierr)
    else
      open(15, status='old', form='formatted', action='read', iostat=ierr)
    end if
    if (ierr /= 0) then
       errmsg = 'warning: cannot find file noahmptable.tbl'
       errflg = 1
       return
!      write(*,'("warning: cannot find file noahmptable.tbl")')
    endif
    read(15,noahmp_optional_parameters)
    close(15)

    ! assign values
!    sr2006_theta_1500t_a_table = sr2006_theta_1500t_a
!    sr2006_theta_1500t_b_table = sr2006_theta_1500t_b
!    sr2006_theta_1500t_c_table = sr2006_theta_1500t_c
!    sr2006_theta_1500t_d_table = sr2006_theta_1500t_d
!    sr2006_theta_1500t_e_table = sr2006_theta_1500t_e
!    sr2006_theta_1500t_f_table = sr2006_theta_1500t_f
!    sr2006_theta_1500t_g_table = sr2006_theta_1500t_g
!    sr2006_theta_1500_a_table  = sr2006_theta_1500_a
!    sr2006_theta_1500_b_table  = sr2006_theta_1500_b
!    sr2006_theta_33t_a_table   = sr2006_theta_33t_a
!    sr2006_theta_33t_b_table   = sr2006_theta_33t_b
!    sr2006_theta_33t_c_table   = sr2006_theta_33t_c
!    sr2006_theta_33t_d_table   = sr2006_theta_33t_d
!    sr2006_theta_33t_e_table   = sr2006_theta_33t_e
!    sr2006_theta_33t_f_table   = sr2006_theta_33t_f
!    sr2006_theta_33t_g_table   = sr2006_theta_33t_g
!    sr2006_theta_33_a_table    = sr2006_theta_33_a
!    sr2006_theta_33_b_table    = sr2006_theta_33_b
!    sr2006_theta_33_c_table    = sr2006_theta_33_c
!    sr2006_theta_s33t_a_table  = sr2006_theta_s33t_a
!    sr2006_theta_s33t_b_table  = sr2006_theta_s33t_b
!    sr2006_theta_s33t_c_table  = sr2006_theta_s33t_c
!    sr2006_theta_s33t_d_table  = sr2006_theta_s33t_d
!    sr2006_theta_s33t_e_table  = sr2006_theta_s33t_e
!    sr2006_theta_s33t_f_table  = sr2006_theta_s33t_f
!    sr2006_theta_s33t_g_table  = sr2006_theta_s33t_g
!    sr2006_theta_s33_a_table   = sr2006_theta_s33_a
!    sr2006_theta_s33_b_table   = sr2006_theta_s33_b
!    sr2006_psi_et_a_table      = sr2006_psi_et_a
!    sr2006_psi_et_b_table      = sr2006_psi_et_b
!    sr2006_psi_et_c_table      = sr2006_psi_et_c
!    sr2006_psi_et_d_table      = sr2006_psi_et_d
!    sr2006_psi_et_e_table      = sr2006_psi_et_e
!    sr2006_psi_et_f_table      = sr2006_psi_et_f
!    sr2006_psi_et_g_table      = sr2006_psi_et_g
!    sr2006_psi_e_a_table       = sr2006_psi_e_a
!    sr2006_psi_e_b_table       = sr2006_psi_e_b
!    sr2006_psi_e_c_table       = sr2006_psi_e_c
!    sr2006_smcmax_a_table      = sr2006_smcmax_a
!    sr2006_smcmax_b_table      = sr2006_smcmax_b

  end subroutine read_mp_table_parameters

end module noahmp_tables

