cmake_minimum_required(VERSION 3.10)
project(NOAHMP LANGUAGES Fortran)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
enable_language(Fortran)
add_compile_options(-cpp)

#------------------------------------------------------------------------------
# Add source files
# NUOPC cap
list(APPEND _noahmp_cap_files drivers/nuopc/lnd_comp_kind.F90
                              drivers/nuopc/lnd_comp_types.F90
                              drivers/nuopc/lnd_comp_shr.F90
                              drivers/nuopc/lnd_comp_io.F90
                              drivers/nuopc/lnd_comp_domain.F90
                              drivers/nuopc/lnd_comp_import_export.F90
                              drivers/nuopc/lnd_comp_nuopc.F90
                              drivers/nuopc/lnd_comp_driver.F90)

# CCPP interface 
list(APPEND _noahmp_ccpp_files drivers/ccpp/noahmpdrv.F90
                               drivers/ccpp/sfc_diff.f
                               drivers/ccpp/machine.F
                               drivers/ccpp/noahmp_tables.f90
                               drivers/ccpp/namelist_soilveg.f
                               drivers/ccpp/set_soilveg.f
                               drivers/ccpp/funcphys.f90
                               drivers/ccpp/physcons.F90)

# NoahMP
list(APPEND _noahmp_files src/module_sf_noahmplsm.F90
                          src/module_sf_noahmp_glacier.F90)

#------------------------------------------------------------------------------
# Build options 
set(OPENMP OFF CACHE BOOL "Enable OpenMP threading")

#------------------------------------------------------------------------------
# Find dependencies
find_package(MPI REQUIRED)
if(OPENMP)
  find_package(OpenMP REQUIRED)
endif()
find_package(NetCDF REQUIRED)
if (TARGET esmf)
  message("Target esmf is already found. Skip find_package ...")
else()
  set(USE_ESMF_STATIC_LIBS TRUE)
  find_package(ESMF REQUIRED)
  if (ESMF_FOUND)
    set(CMAKE_Fortran_COMPILER "${ESMF_F90COMPILER}")
    set(CMAKE_Fortran_FLAGS "${ESMF_F90COMPILEOPTS} ${ESMF_F90COMPILEPATHS}")
    set(CMAKE_Fortran_LINK_FLAGS "${ESMF_F90LINKOPTS} ${ESMF_F90LINKPATHS} ${ESMF_F90LINKRPATHS} ${ESMF_F90ESMFLINKLIBS}")
  endif()
endif()
find_package(FMS REQUIRED COMPONENTS R8)
add_library(fms ALIAS FMS::fms_r8)

#------------------------------------------------------------------------------
# Set CCPP flags for C/C++/Fortran preprocessor
add_definitions(-DCCPP)

#------------------------------------------------------------------------------
# NOAHMP
add_library(noahmp STATIC ${_noahmp_cap_files} ${_noahmp_ccpp_files} ${_noahmp_files})
set_target_properties(noahmp PROPERTIES 
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include
)
target_include_directories(noahmp PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
                                         $<INSTALL_INTERFACE:include>)
target_link_libraries(noahmp PUBLIC esmf fms)

#------------------------------------------------------------------------------
# CMake export for ESMX driver
file(WRITE  ${CMAKE_INSTALL_PREFIX}/lib/cmake/noahmp-esmx.cmake "add_library(FMS::fms_r8 STATIC IMPORTED)\n")
file(APPEND ${CMAKE_INSTALL_PREFIX}/lib/cmake/noahmp-esmx.cmake "set_target_properties(FMS::fms_r8 PROPERTIES\n")
file(APPEND ${CMAKE_INSTALL_PREFIX}/lib/cmake/noahmp-esmx.cmake "  IMPORTED_LOCATION \"${FMS_R8_LIBRARIES}\"\n") 
file(APPEND ${CMAKE_INSTALL_PREFIX}/lib/cmake/noahmp-esmx.cmake "  INTERFACE_INCLUDE_DIRECTORIES \"${FMS_INSTALL_PREFIX}/include\"\n") 
file(APPEND ${CMAKE_INSTALL_PREFIX}/lib/cmake/noahmp-esmx.cmake "  INTERFACE_LINK_DIRECTORIES \"${FMS_INSTALL_PREFIX}/lib\"\n") 
file(APPEND ${CMAKE_INSTALL_PREFIX}/lib/cmake/noahmp-esmx.cmake ")\n\n")

file(APPEND ${CMAKE_INSTALL_PREFIX}/lib/cmake/noahmp-esmx.cmake "add_library(noahmp STATIC IMPORTED)\n")
file(APPEND ${CMAKE_INSTALL_PREFIX}/lib/cmake/noahmp-esmx.cmake "set_target_properties(noahmp PROPERTIES\n")
file(APPEND ${CMAKE_INSTALL_PREFIX}/lib/cmake/noahmp-esmx.cmake "  IMPORTED_LOCATION \"${CMAKE_INSTALL_PREFIX}/lib/libnoahmp.a\"\n")
file(APPEND ${CMAKE_INSTALL_PREFIX}/lib/cmake/noahmp-esmx.cmake "  INTERFACE_INCLUDE_DIRECTORIES \"${CMAKE_INSTALL_PREFIX}/include\"\n")
file(APPEND ${CMAKE_INSTALL_PREFIX}/lib/cmake/noahmp-esmx.cmake "  INTERFACE_LINK_DIRECTORIES  \"${CMAKE_INSTALL_PREFIX}/lib\"\n")
file(APPEND ${CMAKE_INSTALL_PREFIX}/lib/cmake/noahmp-esmx.cmake "  INTERFACE_LINK_LIBRARIES  \"noahmp;FMS::fms_r8\"\n")
file(APPEND ${CMAKE_INSTALL_PREFIX}/lib/cmake/noahmp-esmx.cmake ")\n")
file(APPEND ${CMAKE_INSTALL_PREFIX}/lib/cmake/noahmp-esmx.cmake "target_link_libraries(esmx_driver PUBLIC noahmp)\n")

#------------------------------------------------------------------------------
# Install
install(
  TARGETS noahmp
  EXPORT  noahmp
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include DESTINATION ${CMAKE_INSTALL_PREFIX})

install(EXPORT noahmp
	DESTINATION lib/cmake)
