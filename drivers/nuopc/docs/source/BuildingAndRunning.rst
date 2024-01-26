.. _BuildingAndRunning:

**************************************************************
Building and Running the UFS Weather Model Land Configurations
**************************************************************

==================================
Downloading the Weather Model Code
==================================

To clone the ufs-weather-model repository that is capable to run land as a seperate component, execute the following commands:

.. code-block:: console

  git clone --recursive https://github.com/ufs-community/ufs-weather-model.git

==========================
Building the Weather Model
==========================

The UFS Weather Model uses the CMake build system. The build system is updated to include new land model related applications: `LND`, `ATML` and `S2SWAL`.

.. list-table:: List of applications that supports external alnd component
   :widths: 25 50 50
   :header-rows: 1

   * - Application Name
     - Enabled Components
     - Short Description 
   * - LND
     - CDEPS, NOAHMP, CMEPS
     - Land model forced by GSWP3 and ERA5 data atmosphere
   * - ATML
     - FV3ATM, NOAHMP and CMEPS
     - Land model forced by active atmosphere (FV3ATM) 
   * - S2SWAL
     - FV3ATM, MOM6, CICE6, WW3, NOAHMP, FMS, CMEPS
     - All components used by S2S plus NOAHMP. There is no RT to test this configuration.

To compile the model with land component model support, following command can be used for NCAR's Cheyenne. Note that this is using `contol_p8` as an example for CCPP suites and uses GNU compiler. The platform definition can be changed to build model in other platforms such as MSU's Orion.

.. code-block:: console

  cd tests
  ./compile.sh "cheyenne.gnu" "-DAPP=ATML -DCCPP_SUITES=FV3_GFS_v16,FV3_GFS_v15_thompson_mynn,FV3_GFS_v17_p8,FV3_GFS_v17_p8_rrtmgp,FV3_GFS_v15_thompson_mynn_lam3km" noahmp NO NO

========================================
Running NoahMP Specific Regression Tests
========================================

The new regression tests that inclulde NoahMP land component:

.. list-table:: List of regression tests 
   :widths: 25 50
   :header-rows: 1

   * - Test Name
     - Short Description
   * - datm_cdeps_lnd_gswp3
     - NoahMP forced by the CDEPS "data atmosphere" using Global Soil Wetness Project v3 forcings. 24 hour forecast with 1 hour coupling interval.
   * - datm_cdeps_lnd_era5
     - NoahMP forced by the CDEPS "data atmosphere" using ECMWF's ERA5 ranalysis forcings. 6 hour forecast with 1 hour coupling interval.
   * - datm_cdeps_lnd_era5_rst
     - Restart reproducibility test (compare results with datm_cdeps_lnd_era5)
   * - control_p8_atmlnd_sbs
     - Side-by-side test that forces external land component with active atmosphere in `contol_p8` configuration. This is mainly used to compare land output coming from CCPP/Physics with external NoahMP. In this configuration there is no feedback to active atmosphere.
   * - control_p8_atmlnd
     - Fully coupled atmosphere-land configuration
   * - control_restart_p8_atmlnd
     - Restart reproducibility test for fully coupled atmosphere-land configuration (compare results with control_p8_atmlnd).

Newly introduced RTs can be run with following command,

.. code-block:: console

  cd tests
  ./rt.sh -k -n datm_cdeps_lnd_gswp3

.. code-block:: console

  cd tests
  ./rt.sh -k -n datm_cdeps_lnd_era5_rst

.. code-block:: console

  cd tests
  ./rt.sh -k -n control_p8_atmlnd_sbs

The ``-k`` basically keeps the run directory for future reference. The user could use this run directory to modify out-of-box configuration for different purposes.
