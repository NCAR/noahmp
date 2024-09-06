cmake_minimum_required(VERSION 3.17)
project(NoahMP LANGUAGES Fortran)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake-modules)

set(NOAHMP_MASTER_PROJECT OFF)

if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    set(NOAHMP_MASTER_PROJECT ON)
endif()

set(CMAKE_PREFIX_PATH ${NETCDF_DIR} ${CMAKE_PREFIX_PATH})
set(NETCDF_F90 "YES")
find_package(NetCDF REQUIRED)

file(GLOB NoahMP_Sources
    ${CMAKE_CURRENT_SOURCE_DIR}/src/*.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/utility/*.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/drivers/hrldas/*.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/hrldas_utility_routines/*.F90
)
add_library(NoahMP ${NoahMP_Sources})

if("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
    target_compile_options(NoahMP PRIVATE "-ffree-line-length-none")
endif()

target_include_directories(NoahMP PRIVATE ${NETCDF_INCLUDES})
target_link_libraries(NoahMP PRIVATE ${NETCDF_LIBRARIES_F90})
target_include_directories(NoahMP INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/drivers/hrldas)
add_library(NoahMP::NoahMP ALIAS NoahMP)
