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
     - Subroutines that are used to create ESMF Grid/Mesh using mosaic or ESMF Mesh files 
   * - lnd_comp_driver.F90
     - Subroutines (init, run and finalize) that are used to drive the NoahMP land model
   * - lnd_comp_import_export.F90
     - Subroutines to define import and export states as well as diagnostic code to check the fields in states 
   * - lnd_comp_io.F90
     - Subroutines to read in static information and initial conditions from tiled files and writing model output in tiled format
   * - lnd_comp_kind.F90
     - Parameters that define used kind types
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
     - The path and name of the mosaic grid file (i.e. INPUT/C96_mosaic.nc)
   * - input_dir
     - The directory that stores initial conditions, static information and grid related files 
   * - ic_type
     - Indicates the source of the initial conditions. Two options are supported 'custom' (i.e. C96.initial.tile[1-6].nc) and 'sfc' (default, sfc_data.tile[1-6].nc).
   * - layout
     - Defines decompositions in each direction on each tile (i.e. 3:8 for C96). This needs to be consistent with resolution.
   * - num_soil_levels
     - Number of soil levels used by NoahMP Land Model (i.e. 4)
   * - forcing_height
     - Height of the atmospheric forcing in meters (i.e. 10)
   * - soil_level_thickness
     - Thickness of the soil levels (i.e. 0.10:0.30:0.60:1.00). The number of soil level thickness needs to be consistent with `num_soil_levels`.
   * - soil_level_nodes
     - Depths of the node points for each soil level (i.e. 0.05:0.25:0.70:1.50). The number of node points needs to be consistent with `num_soil_levels`.
   * - dynamic_vegetation_option
     - Options for dynamic vegetation (i.e. 4)
   * - canopy_stomatal_resistance_option
     - Canopy stomatal resistance (1->ball-berry; 2->jarvis)
   * - soil_wetness_option
     - Options for soil moisture factor for stomatal resistance (1->noah; 2->clm; 3->ssib)
   * - runoff_option
     - Options for runoff and groundwater (1->simgm; 2->simtop; 3->schaake96; 4->bats)
   * - surface_exchange_option
     - Options for surface layer drag coeff (ch and cm) (1->m-o; 2->chen97)
   * - supercooled_soilwater_option
     - Options for supercooled liquid water (1-> ny06; 2->koren99)
   * - frozen_soil_adjust_option
     - Options for frozen soil permeability (1-> ny06; 2->koren99)
   * - radiative_transfer_option
     - Options for radiation transfer (1->gap=f(3d,cosz); 2->gap=0; 3->gap=1-fveg)
   * - snow_albedo_option
     - Options for snow surface albedo (1->bats; 2->class)
   * - precip_partition_option
     - Options for rainfall & snowfall (1-jordan91; 2->bats; 3->noah)
   * - soil_temp_lower_bdy_option
     - Options for lower boundary of soil temperature (1->zero-flux; 2->noah)
   * - soil_temp_time_scheme_option
     - Options for snow/soil temperature time scheme (only layer 1, 1->semi-implicit; 2->full implicit/original noah)
   * - surface_evap_resistance_option
     - Options for surface resistance (1->sakaguchi/zeng; 2->seller; 3->mod sellers; 4->1+snow)
   * - glacier_option
     - Options for glacier treatment (1->phase change; 2->simple)
   * - output_freq
     - Options for output frequency in hours

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

During the `InitializeAdvertise` phase, call is made to `fms_init()` to use `Flexible Modeling System (FMS) <https://www.gfdl.noaa.gov/fms/>`_ for reading and writing cubed-sphere tiled output. In this case, The MPI communicator is pulled in through the ESMF VM object and used by FMS. This phase also calls `advertise_fields()` to setup import and export states.

---
Run
---

During the `ModelAdvance` phase, the `cap` updates the import state and calls NoahMP driver routine (`drv_run`, which is found in `drivers/nuopc/lnd_comp_driver.F90`) to run the model and updates the export state with the information calculated by model. The `drv_run` call mainly read in static information as well as initial conditions when it is first called and interpolate monthly data provided by the static information such as fractional coverage of green vegetation and surface albedo to the date of the simulation. Then calculates solar zenith angle based on the time information extracted from `cap` and calls `noahmpdrv_run` subroutine provided by the NoahMP. This phase also responsible to write NoahMP model output in tiled format by taking advantage of FMS and ESMF routines. 

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
   * - inst_height_lowest
     - m
     - noahmp%forc%hgt
     - bottom layer height
     - namelist option `forcing_height`
   * - inst_temp_height_lowest
     - K
     - noahmp%forc%t1
     - bottom layer temperature
     -
   * - inst_pres_height_lowest
     - Pa
     - noahmp%forc%ps
     - pressure at lowest model layer
     -
   * - inst_spec_humid_height_lowest
     - kg kg-1 
     - noahmp%forc%q1
     - bottom layer specific humidity
     -
   * - inst_zonal_wind_height_lowest
     - m s-1 
     - noahmp%forc%u1
     - bottom layer zonal wind
     -
   * - inst_merid_wind_height_lowest
     - m s-1 
     - noahmp%forc%v1
     - bottom layer meridional wind
     -
   * - mean_down_sw_flx
     - W m-2 
     - noahmp%forc%dswsfc
     - mean downward SW heat flux
     -
   * - mean_down_lw_flx
     - W m-2 
     - noahmp%forc%dlwflx
     - mean downward LW heat flux
     -
   * - Faxa_rainc
     - kg m-2 s-1
     - noahmp%forc%tprcpc
     - convective precipitation
     - provided when coupled with CDEPS data atmosphere
   * - Faxa_rainl
     - kg m-2 s-1
     - noahmp%forc%tprcpl
     - large-scale precipitation
     - provided when coupled with CDEPS data atmosphere
   * - mean_prec_rate
     - kg m-2 s-1
     - noahmp%forc%tprcp
     - total precipitation
     - total precipitation can be calculated by its components (tprcpc and tprcpl)

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
   * - Sl_t
     - K
     - noahmp%model%t2mmp
     - land surface temperature
     - 
