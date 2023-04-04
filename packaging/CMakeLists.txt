include(CMakePackageConfigHelpers)

install(
  TARGETS 
    FATODE ADJ TLM FWD
  EXPORT 
    fatode_Exports
  LIBRARY DESTINATION ${INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}
)

# install the mod files
install(
  DIRECTORY 
    ${FATODE_MOD_DIR}/
  DESTINATION 
    ${INSTALL_MOD_DIR}
  FILES_MATCHING PATTERN "*.mod"
)

# install the cmake config files
set(cmake_config_install_location ${INSTALL_PREFIX}/cmake/fatode-${PROJECT_VERSION})

install(
  EXPORT 
    fatode_Exports
  DESTINATION 
    ${cmake_config_install_location}
  NAMESPACE csl::
)

configure_package_config_file(
  "${PROJECT_SOURCE_DIR}/cmake/fatodeConfig.cmake.in"
  "${PROJECT_BINARY_DIR}/fatodeConfig.cmake"
  INSTALL_DESTINATION
    ${cmake_config_install_location}
)

write_basic_package_version_file(
  "${PROJECT_BINARY_DIR}/fatodeConfigVersion.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion
)

install(
  FILES
    ${PROJECT_BINARY_DIR}/fatodeConfig.cmake
    ${PROJECT_BINARY_DIR}/fatodeConfigVersion.cmake
  DESTINATION
    ${cmake_config_install_location}
)

######################################################################
# uninstall target

# https://gitlab.kitware.com/cmake/community/-/wikis/FAQ#can-i-do-make-uninstall-with-cmake
if(NOT TARGET uninstall)
  configure_file(
    "${PROJECT_SOURCE_DIR}/cmake/cmake_uninstall.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
    IMMEDIATE @ONLY)

  add_custom_target(uninstall
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
endif()