################################################################################
# Add subdirectories


if(USE_SUITESPARSE_UMF)
  add_subdirectory(UMFPACK)
endif()

if(USE_SUPER_LU)
  add_subdirectory(SUPERLU)
endif()