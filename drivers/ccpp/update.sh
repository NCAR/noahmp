#!/bin/bash

cp -f ../../../../FV3/ccpp/physics/physics/tools/funcphys.f90 .
md5sum ../../../../FV3/ccpp/physics/physics/tools/funcphys.f90 funcphys.f90
cp -f ../../../../FV3/ccpp/physics/physics/hooks/machine.F .
md5sum ../../../../FV3/ccpp/physics/physics/hooks/machine.F machine.F
cp -f ../../../../FV3/ccpp/physics/physics/hooks/physcons.F90 .
md5sum ../../../../FV3/ccpp/physics/physics/hooks/physcons.F90 physcons.F90
cp -f ../../../../FV3/ccpp/physics/physics/SFC_Models/Land/Noah/namelist_soilveg.f .
md5sum ../../../../FV3/ccpp/physics/physics/SFC_Models/Land/Noah/namelist_soilveg.f namelist_soilveg.f
cp -f ../../../../FV3/ccpp/physics/physics/SFC_Models/Land/Noah/set_soilveg.f .
md5sum ../../../../FV3/ccpp/physics/physics/SFC_Models/Land/Noah/set_soilveg.f set_soilveg.f
cp -f ../../../../FV3/ccpp/physics/physics/SFC_Layer/UFS/sfc_diff.f .
md5sum ../../../../FV3/ccpp/physics/physics/SFC_Layer/UFS/sfc_diff.f sfc_diff.f
cp -f ../../../../FV3/ccpp/physics/physics/SFC_Models/Land/Noahmp/noahmp_tables.f90 .
md5sum ../../../../FV3/ccpp/physics/physics/SFC_Models/Land/Noahmp/noahmp_tables.f90 noahmp_tables.f90
cp -f ../../../../FV3/ccpp/physics/physics/SFC_Models/Land/Noahmp/noahmpdrv.F90 .
md5sum ../../../../FV3/ccpp/physics/physics/SFC_Models/Land/Noahmp/noahmpdrv.F90 noahmpdrv.F90
cp -f ../../../../FV3/ccpp/physics/physics/SFC_Models/Land/Noahmp/noahmpdrv.meta . 
md5sum ../../../../FV3/ccpp/physics/physics/SFC_Models/Land/Noahmp/noahmpdrv.meta noahmpdrv.meta
cp -f ../../../../FV3/ccpp/physics/physics/SFC_Models/Land/Noahmp/module_sf_noahmp_glacier.F90 ../../src/.
md5sum ../../../../FV3/ccpp/physics/physics/SFC_Models/Land/Noahmp/module_sf_noahmp_glacier.F90 ../../src/module_sf_noahmp_glacier.F90
cp -f ../../../../FV3/ccpp/physics/physics/SFC_Models/Land/Noahmp/module_sf_noahmplsm.F90 ../../src/.
md5sum ../../../../FV3/ccpp/physics/physics/SFC_Models/Land/Noahmp/module_sf_noahmplsm.F90 ../../src/module_sf_noahmplsm.F90
