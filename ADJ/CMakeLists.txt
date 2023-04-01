set(ADJ_MOD_DIR "${FATODE_MOD_DIR}/ADJ")

add_library(ADJ)
add_library(csl::ADJ ALIAS ADJ)

set_target_properties(ADJ PROPERTIES 
  ARCHIVE_OUTPUT_DIRECTORY ${FATODE_LIB_DIR}
  VERSION ${PROJECT_VERSION}
  SOVERSION ${PROJECT_VERSION_MAJOR}
)

target_include_directories(ADJ
  PUBLIC  
    $<BUILD_INTERFACE:${ADJ_MOD_DIR}>
    $<INSTALL_INTERFACE:${FATODE_INSTALL_MOD_DIR}>
)

add_library(ADJ_OBJECTS OBJECT)

set_target_properties(ADJ_OBJECTS PROPERTIES
  Fortran_MODULE_DIRECTORY ${ADJ_MOD_DIR}
)

if (USE_SUITESPARSE_LU)
  target_compile_definitions(ADJ_OBJECTS PUBLIC $<$<BOOL:${USE_SUITESPARSE_LU}>:SPARSE_LU>)
endif()
if (USE_SUITESPARSE_UMF)
  target_compile_definitions(ADJ_OBJECTS PUBLIC $<$<BOOL:${USE_SUITESPARSE_UMF}>:SPARSE_UMF>)
endif()
if (USE_LAPACK)
  target_compile_definitions(ADJ_OBJECTS PUBLIC $<$<BOOL:${USE_LAPACK}>:FULL_ALGEBRA>)
endif()

add_subdirectory(ERK_ADJ)
add_subdirectory(RK_ADJ)
add_subdirectory(ROS_ADJ)
add_subdirectory(SDIRK_ADJ)

target_sources(ADJ PRIVATE $<TARGET_OBJECTS:ADJ_OBJECTS>)