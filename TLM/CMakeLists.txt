add_library(TLM)
add_library(csl::TLM ALIAS TLM)

set(TLM_MOD_DIR "${FATODE_MOD_DIR}/TLM")

set_target_properties(TLM PROPERTIES
  Fortran_MODULE_DIRECTORY ${TLM_MOD_DIR}
)

target_include_directories(TLM
  PUBLIC  
    $<BUILD_INTERFACE:${TLM_MOD_DIR}>
    $<INSTALL_INTERFACE:${FATODE_INSTALL_MOD_DIR}>
)

set_target_properties(TLM PROPERTIES 
  ARCHIVE_OUTPUT_DIRECTORY ${FATODE_LIB_DIR}
  VERSION ${PROJECT_VERSION}
  SOVERSION ${PROJECT_VERSION_MAJOR}
)

if (USE_SUITESPARSE)
  target_compile_definitions(TLM PUBLIC $<$<BOOL:${USE_SUITESPARSE}>:SPARSE_UMF>)
endif()
if (USE_LAPACK)
  target_compile_definitions(TLM PUBLIC $<$<BOOL:${USE_LAPACK}>:FULL_ALGEBRA>)
endif()

add_subdirectory(ERK_TLM)
add_subdirectory(RK_TLM)
add_subdirectory(ROS_TLM)
add_subdirectory(SDIRK_TLM)