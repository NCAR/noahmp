.. _DriverAndCoupling:

********************
Driver and Coupling 
********************

This section describes the NoahMP land model National Unified Operation Prediction Capability (NUOPC) "cap", which is a light weight software layer that is required when the model is used as a part of the coupled modeling systems compatible with NOUPC such as `UFS Weather Model <https://ufs-weather-model.readthedocs.io/en/ufs-v2.0.0/>`_ 

=================
Supported Drivers
=================

The driver and coupling layer is found under **drivers/** directory. 

Within the **drivers/nuopc/** directory, the following files are found,

.. list-table:: List of files under **drivers/nuopc/** directory
   :widths: 25 50
   :header-rows: 1

   * - File 
     - Short Description
   * - lnd_comp_domain.F90
     - Includes subroutines/functions that are used to create ESMF Grid/Mesh using mosaic or ESMF Mesh files 
   * - lnd_comp_driver.F90
     - Includes subroutines (init, run and finalize) that are used to drive the NoahMP land model
   * - lnd_comp_import_export.F90
     - Includes subroutines to define import and export states as well as diagnostic code to check the fields in states 
   * - lnd_comp_io.F90
     - Includes subroutines to read in static information and initial conditions from tiled files and writing model output in tiled format. 
   * - lnd_comp_kind.F90
     - Stores parameters that define used kind types etc.
   * - lnd_comp_nuopc.F90
     - Top level file that defines NUOPC cap for NoahMP Land Model
   * - lnd_comp_shr.F90
     - Shared code that aims to read parameters from ESMF config file and manuplating strings
   * - lnd_comp_types.F90
     - Definition of NoahMP Land component specific data types for forcing, initial condition, model etc.

=====================================================
NoahMP Land Model Specific ESMF Configuration Options
=====================================================

The NUOPC "cap" uses set of namelist options provided as ESMF Config file format to define domain and configure the NoahMP Land Model. The list of relevant namelist options, which needs to be defined under `LND_attributes::` group can be seen in the following table.

.. list-table:: Namelist options for NoahMP Land Model
   :widths: 25 50
   :header-rows: 1

   * - Option
     - Short Description
   * - mosaic_file
     - The path and name of the mosaic grid file (i.e. INPUT/C96_mosaic.nc) in Grid Spec format
   * - input_dir
     - The directory that stores initial conditions, static information and grid related files that are pointed by mosaic file
   * - ic_type
     - Indicates the source of the initial conditions. Two options are supported 'custom' (i.e. C96.initial.tile[1-6].nc) and 'sfc' (default, sfc_data.tile[1-6].nc).
   * - layout
     - Defines decompositions in each direction on each tile (i.e. 3:8 for C96). This needs to be consistent with resolution.
   * - num_soil_levels
     - Number of soil levels used by NoahMP Land Model (i.e. 4)
   * - forcing_height
     - Height of the atmospheric forcing in meters (i.e. 10)
   * - soil_level_thickness (Note: soil thickness are hard-coded in the NoahMP at this point)
     - Thickness of the soil levels (i.e. 0.10:0.30:0.60:1.00). The number of soil level thickness needs to be consistent with `num_soil_levels`.
   * - soil_level_nodes
     - Depths of the node points for each soil level (i.e. 0.05:0.25:0.70:1.50). The number of node points needs to be consistent with `num_soil_levels`.
   * - dynamic_vegetation_option (idveg)
     - Options for dynamic vegetation (i.e. 4)
   * - canopy_stomatal_resistance_option (iopt_crs)
     - Canopy stomatal resistance (1->ball-berry; 2->jarvis)
   * - soil_wetness_option (iopt_btr)
     - Options for soil moisture factor for stomatal resistance (1->noah; 2->clm; 3->ssib)
   * - runoff_option (iopt_run)
     - Options for runoff and groundwater (1->simgm; 2->simtop; 3->schaake96; 4->bats)
   * - surface_exchange_option (iopt_sfc)
     - Options for surface layer drag coeff (ch and cm) (1->m-o; 2->chen97)
   * - supercooled_soilwater_option (iopt_frz)
     - Options for supercooled liquid water (1-> ny06; 2->koren99)
   * - frozen_soil_adjust_option (iopt_inf)
     - Options for frozen soil permeability (1-> ny06; 2->koren99)
   * - radiative_transfer_option (iopt_rad)
     - Options for radiation transfer (1->gap=f(3d,cosz); 2->gap=0; 3->gap=1-fveg)
   * - snow_albedo_option (iopt_alb)
     - Options for snow surface albedo (1->bats; 2->class)
   * - precip_partition_option (iopt_snf)
     - Options for rainfall & snowfall (1-jordan91; 2->bats; 3->noah)
   * - soil_temp_lower_bdy_option (iopt_tbot)
     - Options for lower boundary of soil temperature (1->zero-flux; 2->noah)
   * - soil_temp_time_scheme_option (iopt_stc)
     - Options for snow/soil temperature time scheme (only layer 1, 1->semi-implicit; 2->full implicit/original noah)
   * - surface_evap_resistance_option (iopt_rsf, Note: this is not used currently and fixed to 4 in `noahmpdrv.F90`)
     - Options for surface resistance (1->sakaguchi/zeng; 2->seller; 3->mod sellers; 4->1+snow)
   * - glacier_option (iopt_gla, Note: this is not used currently and fixed to 2 in `noahmpdrv.F90`)
     - Options for glacier treatment (1->phase change; 2->simple)
   * - output_freq
     - Option for output frequency in seconds (i.e. 21600 for 6-hourly output, 1 output every coupling time-step and can be used for debugging)
   * - restart_freq
     - Option for restart frequency in seconds. If it is not provided, it will be same with `output_freq`
   * - restart_file
     - Option for specifying the restart file (`.tile#.nc` will be added to the given file name). If it is not provided the model specifies the file name internally using current model time and coupling time step. 
   * - do_mynnedmf
     - Option for MYNN-EDMF (default value is `.false.`)
   * - do_mynnsfclay
     - Option for MYNN surface layer scheme (default value is `.false.`)
   * - soil_type_category (isot)
     - Option for soil type (0->Zobler - 9 category; 1->STATSGO - 19 category; 2->STAS-RUC - 19 category)
   * - veg_type_category (ivegsrc)
     - Option for source of vegetation data (0->USGS; 1->IGBP; 2->UMD; 3->NLCD40; 4->USGS-RUC; 5->MODI-RUC)
   * - initial_emiss
     - Option for initial surface lw emissivity in fraction (default value is 0.95)
   * - initial_albedo
     - Option for initial mean surface albedo (value is default 0.2)
   * - calc_snet
     - Option for calculating net shortwave radiation using downward component and surface albedo

.. note::
   ``:`` symbol is used as a seperator for namelist options with multiple values such as `layout`, `soil_level_thickness`.

===========================
Underlying Model Interfaces
===========================

---------------
Domain Creation
---------------

The current version of the NUOPC "cap" is able to create ESMF grid by reading mosaic grid file. Then, NOUPC "cap" converts created ESMF multi-tile grid to ESMF Mesh to standardize interface and allow running same NUOPC "cap" also reading domain information in ESMF Mesh format (not implemented yet). In this design, only ESMF Mesh is exposed to coupler or mediator via defining fields in import and export states using ESMF Mesh. The land fraction information (`land_frac`) is provided by reading `*oro_data.tile*` files and it is also used to define the land-sea mask (if land_frac is greater than 0 then it is assumed as land, otherwise it is water). The orography data is also defined using `orog_raw` variable in the same files.

--------------
Initialization
--------------

During the `InitializeAdvertise` phase, call is made to `advertise_fields()` to setup import and export states.

---
Run
---

During the `ModelAdvance` phase, the `cap` updates the import state and calls NoahMP driver routine (`drv_run`, which is found in `drivers/nuopc/lnd_comp_driver.F90`) to run the model and updates the export state with the information calculated by model. The `drv_run` call mainly read in static information as well as initial conditions when it is first called and interpolate monthly data provided by the static information such as fractional coverage of green vegetation and surface albedo to the date of the simulation. Then calculates solar zenith angle based on the time information extracted from `cap` and calls `noahmpdrv_run` subroutine provided by the NoahMP. This phase also responsible to write NoahMP model output in tiled format by taking advantage of ESMF I/O multi-tile support. 

.. note::
   : the restart capability is only tested with DATM+NOAHMP configuration.

--------
Finalize
--------

This phase is not implemented yet.

------------------------------
Model Fields Used for Coupling
------------------------------


.. list-table:: Import fields
   :widths: 25 10 10 25 25
   :header-rows: 1

   * - Standard Name
     - Units
     - Model Variable
     - Description
     - Notes
   * - inst_height_lowest (`Sa_z`)
     - m
     - noahmp%forc%hgt
     - bottom layer height
     - namelist option `forcing_height`
   * - inst_temp_height_lowest (`Sa_tbot`)
     - K
     - noahmp%forc%t1
     - bottom layer temperature
     - 
   * - inst_temp_height_lowest_from_phys (`Sa_ta`)
     - K
     - noahmp%forc%t1
     - bottom layer temperature
     - used if coupled with active atmosphere
   * - inst_temp_height_surface (`Sa_tskn`)
     - K
     - noahmp%forc%tskin
     - surface skin temperature
     - 
   * - inst_pres_height_lowest (`Sa_pbot`)
     - Pa
     - noahmp%forc%pbot
     - pressure at lowest model layer
     -
   * - inst_pres_height_lowest_from_phys (`Sa_prsl`)
     - Pa
     - noahmp%forc%pbot
     - pressure at lowest model layer
     - used if coupled with active atmosphere
   * - inst_pres_height_surface (`Sa_pslv`)
     - Pa
     - noahmp%forc%ps
     - surface pressure
     -
   * - inst_spec_humid_height_lowest (`Sa_shum`)
     - kg kg-1 
     - noahmp%forc%q1
     - bottom layer specific humidity
     -
   * - inst_spec_humid_height_lowest_from_phys (`Sa_qa`)
     - kg kg-1 
     - noahmp%forc%q1
     - bottom layer specific humidity
     - used if coupled with active atmosphere
   * - inst_zonal_wind_height_lowest (`Sa_u`)
     - m s-1 
     - noahmp%forc%u1
     - bottom layer zonal wind
     -
   * - inst_merid_wind_height_lowest (`Sa_v`)
     - m s-1 
     - noahmp%forc%v1
     - bottom layer meridional wind
     -
   * - inst_zonal_wind_height_lowest_from_phys (`Sa_ua`)
     - m s-1 
     - noahmp%forc%u1
     - bottom layer zonal wind
     - used if coupled with active atmosphere
   * - inst_merid_wind_height_lowest_from_phys (`Sa_va`)
     - m s-1 
     - noahmp%forc%v1
     - bottom layer meridional wind
     - used if coupled with  active atmosphere
   * - inst_exner_function_height_lowest (`Sa_exner`)
     - 1 
     - noahmp%forc%prslk1
     - dimensionless exner function at surface adjacent layer
     -
   * - surface_friction_velocity (`Sa_ustar`)
     - m s-1 
     - noahmp%forc%ustar1
     - surface friction velocity
     -
   * - mean_down_sw_flx (`Faxa_swdn`)
     - W m-2 
     - noahmp%forc%dswsfc
     - mean downward SW heat flux
     -
   * - mean_down_lw_flx (`Faxa_lwdn`)
     - W m-2 
     - noahmp%forc%dlwflx
     - mean downward LW heat flux
     -
   * - mean_net_sw_flx (`Faxa_swnet`)
     - W m-2 
     - noahmp%forc%dlwflx
     - net SW radiation 
     - if it is not available, it will be calculated by using `mean_down_sw_flx` and surface albedo (see `calc_snet` option)
   * - mean_prec_rate_conv (`Faxa_rainc`)
     - kg m-2 s-1
     - noahmp%forc%tprcpc
     - convective precipitation
     - provided when coupled with data atmosphere
   * - Faxa_rainl
     - kg m-2 s-1
     - noahmp%forc%tprcpl
     - large-scale precipitation
     - provided when coupled with data atmosphere
   * - mean_prec_rate (`Faxa_rain`)
     - kg m-2 s-1
     - noahmp%forc%tprcp
     - total precipitation
     - total precipitation can be calculated by its components (tprcpc and tprcpl) or provided by directly from active atmosphere 
   * - Faxa_snowc
     - kg m-2 s-1
     - noahmp%forc%snowc
     - convective part of snow precipitation
     - 
   * - Faxa_snowl
     - kg m-2 s-1
     - noahmp%forc%snowl
     - large-scale part of snow precipitation
     - 
   * - mean_fprec_rate (`Faxa_snow`)
     - kg m-2 s-1
     - noahmp%forc%snow
     - total snow precipitation
     -
   * - vfrac
     - 1
     - noahmp%forc%vegfrac
     - areal fractional cover of green vegetation
     -
   * - zorl
     - cm
     - noahmp%forc%zorl
     - surface roughness
     -

.. list-table:: Export fields
   :widths: 25 10 10 25 25
   :header-rows: 1

   * - Standard Name
     - Units
     - Model Variable
     - Description
     - Notes
   * - Sl_lfrin
     - 0-1
     - noahmp%domain%frac
     - land fraction
     - required by mediator
   * - Sl_sfrac
     - 0-1
     - noahmp%model%sncovr1
     - mean snow area fraction
     -
   * - Fall_lat
     - kg kg-1 m s-1
     - noahmp%model%evap
     - mean latent heat flux
     -
   * - Fall_sen
     - kg kg-1 m s-1
     - noahmp%model%hflx
     - mean sensible heat flux
     -
   * - Fall_evap
     - W m-2
     - noahmp%model%ep
     - mean potential latent heat flux
     -
   * - Sl_tref
     - K
     - noahmp%model%t2mmp
     - instantenous temperature at 2 meters
     -
   * - Sl_qref
     - kg kg-1
     - noahmp%model%q2mp
     - instantenous specific humidity at 2 meters
     -
   * - Sl_q
     - kg kg-1
     - noahmp%model%qsurf
     - instantenous specific humidity (at lowest model layer)
     -
   * - Fall_gflx
     - W m-2
     - noahmp%model%gflux
     - mean upward heat flux (ground)
     -
   * - Fall_roff
     - kg m-2 s-1
     - noahmp%model%runoff
     - mean runoff rate (surface)
     -
   * - Fall_soff
     - kg m-2 s-1
     - noahmp%model%drain
     - mean runoff rate (sub-surface)
     -
   * - Sl_cmm
     - m s-1
     - noahmp%model%cmm
     - instantenous drag wind speed for momentum
     -
   * - Sl_chh
     - kg m-2 s-1
     - noahmp%model%chh
     - instantenous drag wind speed for heat and moisture
     -
   * - Sl_zvfun
     - 0-1
     - noahmp%model%zvfun
     - instantenous function of roughness length and areal fractional cover of green vegetation
     -
