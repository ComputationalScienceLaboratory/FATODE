################################################################################
# Preamble
cmake_minimum_required(VERSION 3.21)

project(
  FATODE
  VERSION 1.2.0
  LANGUAGES Fortran C
)

message(STATUS "CMake build configuration for ${PROJECT_NAME}(${CMAKE_BUILD_TYPE}) ${PROJECT_VERSION}")

set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${PROJECT_SOURCE_DIR}/cmake")

################################################################################
# Projet wide setup options

option(INCLUDE_FORWARD_SOLVERS "Include the forward model integrators in the library" ON)
option(INCLUDE_ADJOINT_SOLVERS "Include the adjoint solvers in the library" ON)
option(INCLUDE_TANGENT_LINEAR_SOLVERS "Include the tangent linear methods in the library" ON)

option(ENABLE_TESTS "Build tests" ON)

# TODO: instead of using the included files pulled from the host repositories, add an option to link to the libraries
option(USE_SUITESPARSE_UMF "Build with SuiteSpare UMF support" OFF)
option(USE_SUPER_LU "Build with SuperLU SUPPORT" OFF)
option(USE_FULL_ALGEBRA "Set this if you want to do all of the algebra. Otherwise set the SuperLU or UMF options to use a sparse matrix approach." OFF)

include(exclusive_options)
check_mutually_exclusive_options(USE_SUPER_LU USE_SUITESPARSE_UMF USE_FULL_ALGEBRA)

find_package(LAPACK REQUIRED)

if(USE_SUITESPARSE_UMF)
  set(CMAKE_MODULE_PATH "/opt/homebrew/lib/cmake/SuiteSparse;${CMAKE_MODULE_PATH}")
  set(CMAKE_MODULE_PATH "/usr/local/lib64/cmake/SuiteSparse;${CMAKE_MODULE_PATH}")
  find_package(UMFPACK
    MODULE
    REQUIRED
  )
endif()

if(USE_SUPER_LU)
  find_package(superlu REQUIRED)
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

set(FATODE_MOD_DIR "${PROJECT_BINARY_DIR}/include")
set(FATODE_LIB_DIR "${PROJECT_BINARY_DIR}/lib")
include(GNUInstallDirs)

set(INSTALL_PREFIX "${PROJECT_NAME}-${PROJECT_VERSION}")
set(INSTALL_MOD_DIR "${INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}")

################################################################################
# FATODE targets

add_library(FATODE)
add_library(csl::FATODE ALIAS FATODE)

if(INCLUDE_FORWARD_SOLVERS)
  add_subdirectory(FWD)
  target_sources(FATODE PUBLIC $<TARGET_OBJECTS:FWD_OBJECTS>)
  target_link_libraries(FATODE
    PUBLIC
      csl::FWD
  )
endif()

if(INCLUDE_ADJOINT_SOLVERS)
  add_subdirectory(ADJ)
  target_sources(FATODE PUBLIC $<TARGET_OBJECTS:ADJ_OBJECTS>)
  target_link_libraries(FATODE
    PUBLIC
      csl::ADJ
  )
endif()

if(INCLUDE_TANGENT_LINEAR_SOLVERS)
  add_subdirectory(TLM)
  target_sources(FATODE PUBLIC $<TARGET_OBJECTS:TLM_OBJECTS>)
  target_link_libraries(FATODE
    PUBLIC
      csl::TLM
  )
endif()

add_subdirectory(LSS_LIBS)

target_link_libraries(FATODE
  PUBLIC
    LAPACK::LAPACK
)

if(USE_SUITESPARSE_UMF)
  target_link_libraries(FATODE
    PUBLIC
      ${UMFPACK_LIBRARY}
  )
  target_include_directories(FATODE
    PUBLIC
      ${UMFPACK_INCLUDE_DIR}
  )
endif()

if(USE_SUPER_LU)
  target_link_libraries(FATODE
    PUBLIC
      superlu::superlu
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

################################################################################
# FATODE tests

if(PROJECT_IS_TOP_LEVEL AND ENABLE_TESTS)
  enable_testing()
  add_subdirectory(test)
endif()

################################################################################
# Packaging

# only include packaging if we are the top level project being built
if(PROJECT_IS_TOP_LEVEL)
  add_subdirectory(packaging)
endif()