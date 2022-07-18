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
list(APPEND _noahmp_files src/module_sf_noahmplsm.f90
                          src/module_sf_noahmp_glacier.f90)

#------------------------------------------------------------------------------
# Customization for compilers
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(CMAKE_Fortran_FLAGS "-g -fbacktrace -ffree-line-length-none")
  if(${CMAKE_Fortran_COMPILER_VERSION} VERSION_GREATER_EQUAL 10)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fallow-argument-mismatch -fallow-invalid-boz")
  endif()
  set(CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -fcheck=bounds -ffpe-trap=invalid,zero,overflow,underflow" )
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set(CMAKE_Fortran_FLAGS "-g -traceback")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O -assume realloc_lhs")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -fp-model precise")
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -check -check noarg_temp_created -check nopointer -fpe0 -ftrapuv -init=snan,arrays")
else()
  message(WARNING "Fortran compiler with ID ${CMAKE_Fortran_COMPILER_ID} will be used with CMake default options")
endif()

#------------------------------------------------------------------------------
# Find dependencies
#if (TARGET esmf)
#  message("Target esmf is already found. Skip find_package ...")
#else()
  set(USE_ESMF_STATIC_LIBS TRUE)
  find_package(ESMF REQUIRED)
  if (ESMF_FOUND)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${ESMF_F90COMPILEOPTS} ${ESMF_F90COMPILEPATHS}")
    set(CMAKE_Fortran_LINK_FLAGS "${CMAKE_Fortran_LINK_FLAGS} ${ESMF_F90LINKOPTS} ${ESMF_F90LINKPATHS} ${ESMF_F90LINKRPATHS} ${ESMF_F90ESMFLINKLIBS}")
  endif()
  #endif()
find_package(MPI REQUIRED)
find_package(FMS REQUIRED COMPONENTS R4 R8)
add_library(fms ALIAS FMS::fms_r8)

#------------------------------------------------------------------------------
# Set CCPP flags for C/C++/Fortran preprocessor
add_definitions(-DCCPP)

#------------------------------------------------------------------------------
# NOAHMP
add_library(noahmp STATIC ${_noahmp_cap_files} ${_noahmp_ccpp_files} ${_noahmp_files})
set_target_properties(noahmp PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)
target_include_directories(noahmp PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod>
                                       $<INSTALL_INTERFACE:mod>)
target_link_libraries(noahmp PUBLIC esmf fms)

#------------------------------------------------------------------------------
# Install
install(
 TARGETS noahmp
  EXPORT  noahmp-config
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  COMPONENT Library)

install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod DESTINATION ${CMAKE_INSTALL_PREFIX})

install(EXPORT      noahmp-config
        DESTINATION lib/cmake)
