cmake_minimum_required(VERSION 3.17)
project(NoahMP LANGUAGES Fortran CXX C)
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
    ${CMAKE_CURRENT_SOURCE_DIR}/drivers/hrldas/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/hrldas_utility_routines/*.F90
)
add_library(noahmp ${NoahMP_Sources})
add_library(NoahMP::noahmp ALIAS noahmp)

if("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
    target_compile_options(noahmp PRIVATE "-ffree-line-length-none")
endif()

target_include_directories(noahmp PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/drivers/hrldas>)
target_include_directories(noahmp PRIVATE $<BUILD_INTERFACE:${NETCDF_INCLUDES}>)
target_link_libraries(noahmp PRIVATE ${NETCDF_LIBRARIES_F90})
target_link_libraries(noahmp PUBLIC gfortran)
target_include_directories(noahmp INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/drivers/hrldas>)

# Install the library and export the target
install(TARGETS noahmp
    EXPORT NoahMPTargets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    INCLUDES DESTINATION include
)

install(EXPORT NoahMPTargets
    FILE NoahMPTargets.cmake
    NAMESPACE NoahMP::
    DESTINATION lib/cmake/noahmp
)
