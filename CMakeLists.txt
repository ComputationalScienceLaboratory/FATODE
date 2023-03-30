################################################################################
# Preamble
cmake_minimum_required(VERSION 3.21)

project(
  FATODE
  VERSION 1.2.0
  LANGUAGES Fortran
)


message ( STATUS "CMake build configuration for ${PROJECT_NAME}(${CMAKE_BUILD_TYPE}) ${PROJECT_VERSION}" )

################################################################################
# Projet wide setup options

option(INCLUDE_FORWARD_SOLVERS "Include the forward model integrators in the library" ON)
option(INCLUDE_ADJOINT_SOLVERS "Include the adjoint solvers in the library" ON)
option(INCLUDE_TANGENT_LINEAR_SOLVERS "Include the tangent linear methods in the library" ON)

# TODO: instead of using the included files pulled from the host repositories, add an option to link to the libraries
option(USE_SUITESPARSE "Build with suitespare support" ON)
option(USE_LAPACK "Build with lapack support" OFF)

if (USE_SUITESPARSE AND USE_LAPACK)
    message(FATAL_ERROR "Options USE_SUITESPARSE and USE_LAPACK are mutually exclusive. Choose one of them.")
elseif (NOT USE_SUITESPARSE AND NOT USE_LAPACK)
    message(FATAL_ERROR "You must choose one to use one of USE_SUITESPARE or USE_LAPACK.")
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

set(FATODE_MOD_DIR "${PROJECT_BINARY_DIR}/include")
set(FATODE_LIB_DIR "${PROJECT_BINARY_DIR}/lib")
include(GNUInstallDirs)
set(FATODE_INSTALL_PREFIX "${PROJECT_NAME}-${PROJECT_VERSION}" )
set(FATODE_INSTALL_MOD_DIR "${INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}")

################################################################################
# FATODE targets

add_library(FATODE)
add_library(csl::FATODE ALIAS FATODE)

if(INCLUDE_FORWARD_SOLVERS)
  add_subdirectory(FWD)
  target_link_libraries(FATODE
    PUBLIC
      csl::FWD
  )
endif()

if(INCLUDE_ADJOINT_SOLVERS)
  add_subdirectory(ADJ)
  target_link_libraries(FATODE
    PUBLIC
      csl::ADJ
  )
endif()

if(INCLUDE_TANGENT_LINEAR_SOLVERS)
  add_subdirectory(TLM)
  target_link_libraries(FATODE
    PUBLIC
      csl::TLM
  )
endif()

set_target_properties(FATODE PROPERTIES
  Fortran_MODULE_DIRECTORY ${FATODE_MOD_DIR}
)

target_include_directories(FATODE
  PUBLIC  
    $<BUILD_INTERFACE:${FATODE_MOD_DIR}>
    $<INSTALL_INTERFACE:${FATODE_INSTALL_MOD_DIR}>
)

set_target_properties(FATODE PROPERTIES 
  ARCHIVE_OUTPUT_DIRECTORY ${FATODE_LIB_DIR}
  VERSION ${PROJECT_VERSION}
  SOVERSION ${PROJECT_VERSION_MAJOR}
)

