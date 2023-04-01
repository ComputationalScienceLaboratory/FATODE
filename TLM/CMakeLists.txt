set(TLM_MOD_DIR "${FATODE_MOD_DIR}/TLM")

add_library(TLM)
add_library(csl::TLM ALIAS TLM)

set_target_properties(TLM PROPERTIES 
  ARCHIVE_OUTPUT_DIRECTORY ${FATODE_LIB_DIR}
  VERSION ${PROJECT_VERSION}
  SOVERSION ${PROJECT_VERSION_MAJOR}
)

target_include_directories(TLM
  PUBLIC  
    $<BUILD_INTERFACE:${TLM_MOD_DIR}>
    $<INSTALL_INTERFACE:${FATODE_INSTALL_MOD_DIR}>
)

add_library(TLM_OBJECTS OBJECT)

set_target_properties(TLM_OBJECTS PROPERTIES
  Fortran_MODULE_DIRECTORY ${TLM_MOD_DIR}
)

if (USE_SUITESPARSE_LU)
  target_compile_definitions(TLM_OBJECTS PUBLIC $<$<BOOL:${USE_SUITESPARSE_LU}>:SPARSE_LU>)
endif()
if (USE_SUITESPARSE_UMF)
  target_compile_definitions(TLM_OBJECTS PUBLIC $<$<BOOL:${USE_SUITESPARSE_UMF}>:SPARSE_UMF>)
endif()
if (USE_LAPACK)
  target_compile_definitions(TLM_OBJECTS PUBLIC $<$<BOOL:${USE_LAPACK}>:FULL_ALGEBRA>)
endif()

add_subdirectory(ERK_TLM)
add_subdirectory(RK_TLM)
add_subdirectory(ROS_TLM)
add_subdirectory(SDIRK_TLM)

target_sources(TLM PRIVATE $<TARGET_OBJECTS:TLM_OBJECTS>)