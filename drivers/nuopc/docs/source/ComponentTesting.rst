.. _ComponentTesting:

*********************************************************************
Building, Running and Testing Land Model Outside of UFS Weather Model 
*********************************************************************

============================
Creating working environment
============================

This section mainly aims to give brief information about creating work environment to test, build, run and develop NoahMP land component outside of the UFS Weather Model. For this purpose, `Docker <https://www.docker.com>`_ container can be used. 

To create Docker container that will be the foundation of development and testing environment:

.. code-block:: console

  docker run -it ubuntu:latest

  apt-get update
  apt-get -y install unzip file gringo
  apt-get -y install build-essential binutils-dev gfortran
  apt-get -y install python3-dev python3-yaml
  apt-get -y install cmake wget ca-certificates vim git

.. note::
   Using container to create test and development environment is optional. The following steps can be used also in other platforms or/and OS in case of lack of depencenies in there.

Then, `spack <https://spack.io>`_ package manager can be used to install dependencies such as ESMF. In this case, the content of the `spack.yaml` can be structured as following. Also note that the installadtion directory can be changed through the `config` section of the YAML file. 

.. code-block:: yaml 

  spack:
    concretizer:
      unify: true
    specs:
    - esmf@8.5.0b17+external-parallelio
    - parallelio@2.5.10+pnetcdf~shared
    view: /root/.spack-ci/view
    config:
      source_cache: ~/.spack-ci/source_cache
      misc_cache: ~/.spack-ci/misc_cache
      test_cache: ~/.spack-ci/test_cache
      install_tree:
        root: ~/.spack-ci/opt

To install dependencies:

.. code-block:: console

  mkdir $HOME/test
  cd $HOMEtest

  git clone https://github.com/NOAA-EMC/spack.git
  . spack/share/spack/setup-env.sh
  spack -e . concretize -f
  spack -e . install -j3

.. note::
  In this case, the `spack.yaml` is assumed in the test/ directory. In case of having it an other directory, user just need to point the directory after `-e` option (it is set to `.` in this example).

.. note::
  If there is a message like "Fetch anyway?", then just hit "Y".

====================================================
Building Data Atmosphere Forced NoahMP Configuration 
====================================================

Once the working environment is ready, the components of the simple configuration (`datm+lnd`) can be build. In this case, first model components are need to be be build individually.

**NoahMP**

.. code-block:: console

  export PATH=$HOME/.spack-ci/view/bin:$PATH
  export ESMFMKFILE=$HOME/.spack-ci/view/lib/esmf.mk
  export NetCDF_ROOT=`nc-config --prefix`
  export FC=gfortran
  cd $HOME/test
  git clone https://github.com/NOAA-EMC/noahmp
  cd  noahmp
  mkdir build
  cd build
  cmake -DCMAKE_INSTALL_PREFIX=../install -DOPENMP=ON ../
  make
  make install

.. note::
  This will install component specific files (`libnoahmp.a` and `*.cmake`) under installation directory.

**CDEPS**

.. code-block:: console

  export PATH=$HOME/.spack-ci/view/bin:$PATH
  export ESMFMKFILE=$HOME/.spack-ci/view/lib/esmf.mk
  export FC=gfortran
  cd $HOME/test
  git clone -b hotfix/std_build https://github.com/uturuncoglu/CDEPS.git cdeps
  cd cdeps
  mkdir build
  cd build
  cmake -DCMAKE_INSTALL_PREFIX=../ \
    -DPIO_C_LIBRARY=$HOME/.spack-ci/view/lib \
    -DPIO_C_INCLUDE_DIR=$HOME/.spack-ci/view/include \
    -DPIO_Fortran_LIBRARY=$HOME/.spack-ci/view/lib \
    -DPIO_Fortran_INCLUDE_DIR=$HOME/.spack-ci/view/include \
    -DCMAKE_Fortran_FLAGS="-ffree-line-length-none -fallow-argument-mismatch -fallow-invalid-boz" \
    -DDISABLE_FoX=ON ../
  make
  make install

.. note::
  Again, this will install component specific files under installation directory.

Then, component libraries can be used to create executable through the use of `ESMX (Earth System Model eXecutable) <https://github.com/esmf-org/esmf/tree/develop/src/addon/ESMX>`_ layer provided by `ESMF <http://earthsystemmodeling.org/docs/nightly/develop/ESMF_refdoc/>`_. 

This requires creating a simple YAML file (`esmxBuild.yaml`) that points component specific files with following content:

.. code-block:: yaml 

  components:
    datm:
      cmake_config: $HOME/test/cdeps/install/lib/cmake/datm-esmx.cmake
      fort_module: cdeps_datm_comp
    noahmp:
      cmake_config: $HOME/test/noahmp/lib/cmake/noahmp-esmx.cmake
      fort_module: lnd_comp_nuopc

.. note::
  File `esmxBuild.yaml` needs to be placed in the $HOME/test/app directory.


The application can be build with following commands:

.. code-block:: console

  export PATH=$HOME/.spack-ci/view/bin:$PATH
  export ESMFMKFILE=$HOME/.spack-ci/view/lib/esmf.mk
  export ESMF_ESMXDIR=$HOME/.spack-ci/view/include/ESMX
  cd $HOME/test
  mkdir app
  cd app
  cmake -H$ESMF_ESMXDIR -Bbuild
  cd build
  make

.. note::
  This will create `esmx` executable under build/ directory.

===================================================
Running Data Atmosphere Forced NoahMP Configuration
===================================================

To run `esmx` executable, the run directory, which includes input files for both CDEPS and NoahMP components, is required. The more information about running UFS Weather Model datm+lnd configuration can be found in :ref:`BuildingAndRunning`.

===========================
Automated Component Testing
===========================

The NoahMP model uses GitHub Actions (GHA), a GitHub-hosted continuous integration service, to perform CI/CD testing. The GHA is triggered in following cases,

#. When a developer makes a pull request (PR) to the NoahMP repository

#. When a developer makes a direct push to default branch (i.e. develop)

#. Twice a week on Monday and Friday (This is required to prevent auto-removing cache entries after 7-days)

#. Repository admin could also trigger GHA manually

The GHA-related ``yaml`` script is located in the ``.github/workflows/`` directory. ``datm_noahmp.yaml`` is the main workflow file that aim to run `datm+lnd` configuration. 

* `.github/workflows/tests/` directory includes YAML file that will be used to create required configuration files for ESMX driver and retrive required input files.
* `.github/workflows/tests/test_datm_lnd/` directory includes YAML files that will be used to create required configuration files for CDEPS and NoahMP and retrive required input files.
* `.github/workflows/data/` directory includes additional the input files (initial conditions) that are not found on the web to retrive.

The action uses composite action for isolated component testing that can be found in `esmf-org nuopc-comp-testing repository <https://github.com/esmf-org/nuopc-comp-testing>`_.
