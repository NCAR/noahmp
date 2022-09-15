#!/bin/bash

lst=`ls -al *.[f-F]* | awk '{print $9}'`
for i in $lst
do
  echo $i
  cp -f ../../../../FV3/ccpp/physics/physics/$i .
  md5sum ../../../../FV3/ccpp/physics/physics/$i $i
done

echo "module_sf_noahmp_glacier.f90"
cp -f ../../../../FV3/ccpp/physics/physics/module_sf_noahmp_glacier.f90 ../../src/.
md5sum ../../../../FV3/ccpp/physics/physics/module_sf_noahmp_glacier.f90 ../../src/module_sf_noahmp_glacier.f90
echo "module_sf_noahmplsm.f90"
cp -f ../../../../FV3/ccpp/physics/physics/module_sf_noahmplsm.f90 ../../src/.
md5sum ../../../../FV3/ccpp/physics/physics/module_sf_noahmplsm.f90 ../../src/module_sf_noahmplsm.f90
