.. _BuildingAndRunning:

**************************************************************
Building and Running the UFS Weather Model Land Configurations
**************************************************************

==================================
Downloading the Weather Model Code
==================================

To clone the ufs-weather-model repository that is capable to run land as a seperate component, execute the following commands:

.. code-block:: console

  git clone -b feature/lnd_noahmp --recursive https://github.com/uturuncoglu/ufs-weather-model.git ufs-weather-model

==========================
Building the Weather Model
==========================

The UFS Weather Model uses the CMake build system. The build system is updated to include new land model related applications: LND and LND-ALL.

.. list-table:: List of regression tests 
   :widths: 25 50
   :header-rows: 1

   * - Application Name
     - Enabled Components
   * - LND
     - CMEPS, CDEPS, NOAHMP and FMS
   * - LND-ALL
     - All components used by S2S (w/o wave) as well as NOAHMP

To compile the model with land model support, following command can be used for NCAR's Cheyenne. The platform definition can be changed to build model in other platforms such as MSU's Orion but initial implementation is only tested on NCAR's Cheyenne at this point.

.. code-block:: console

  cd tests
  ./compile.sh "cheyenne.intel" "-DAPP=LND-ALL -DCCPP_SUITES=FV3_GFS_v16_coupled_nsstNoahmpUGWPv1,FV3_GFS_v16_coupled_p7_rrtmgp" noahmp NO NO

========================================
Running NoahMP Specific Regression Tests
========================================

Two new regression test are included to test the external NoahMP land component:

.. list-table:: List of regression tests 
   :widths: 25 50
   :header-rows: 1

   * - Test Name
     - Short Description
   * - datm_cdeps_lnd_gswp3
     - NoahMP forced by the CDEPS "data atmosphere" using Global Soil Wetness Project v3 forcings. 24 hour forecast with 1 hour coupling interval.
   * - datm_cdeps_lnd_gswp3_rst
     - Restart reproducibility test (compare results with datm_cdeps_lnd_gswp3)

Since the initial implementation is still in a development fork, the required input files and baseline directory to run the new regression tests are not the part of the common UFS Weather Model input directory. As a result, the user need to point temporary input and baseline directory to run newly defined RTs under NCAR's Cheyenne platform. This can be done easily by modifying ``tests/rt.sh`` script as follows,

.. code-block:: console

  # Find Cheyenne section and change the DISKNM variable from
  DISKNM=/glade/p/ral/jntp/GMTB/ufs-weather-model/RT
  # to
  DISKNM=/glade/work/turuncu/NOAHMP/RT

This will allow the RTs find the reqired input and baseline files. Then, RTs can be run with following command,

.. code-block:: console

  cd tests
  ./rt.sh -k -n datm_cdeps_lnd_gswp3

The ``-k`` basically keeps the run directory for future reference. The user could use this run directory to modify out-of-box configuration for different purposes.
