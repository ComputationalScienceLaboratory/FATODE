add_library(FWD)
add_library(csl::FWD ALIAS FWD)

set(FWD_MOD_DIR "${FATODE_MOD_DIR}/FWD")

set_target_properties(FWD PROPERTIES
  Fortran_MODULE_DIRECTORY ${FWD_MOD_DIR}
)

target_include_directories(FWD
  PUBLIC  
    $<BUILD_INTERFACE:${FWD_MOD_DIR}>
    $<INSTALL_INTERFACE:${FATODE_INSTALL_MOD_DIR}>
)

set_target_properties(FWD PROPERTIES 
  ARCHIVE_OUTPUT_DIRECTORY ${FATODE_LIB_DIR}
  VERSION ${PROJECT_VERSION}
  SOVERSION ${PROJECT_VERSION_MAJOR}
)

if (USE_SUITESPARSE_LU)
  target_compile_definitions(FWD PUBLIC $<$<BOOL:${USE_SUITESPARSE_LU}>:SPARSE_LU>)
endif()
if (USE_SUITESPARSE_UMF)
  target_compile_definitions(FWD PUBLIC $<$<BOOL:${USE_SUITESPARSE_UMF}>:SPARSE_UMF>)
endif()
if (USE_LAPACK)
  target_compile_definitions(FWD PUBLIC $<$<BOOL:${USE_LAPACK}>:FULL_ALGEBRA>)
endif()

add_subdirectory(ERK)
add_subdirectory(RK)
add_subdirectory(ROS)
add_subdirectory(SDIRK)