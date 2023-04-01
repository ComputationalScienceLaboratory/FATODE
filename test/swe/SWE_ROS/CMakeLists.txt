create_standard_test(
  NAME swe_ros
  SOURCES swe2D_ros_dr.F90 swe2D_upwind.F90
  DATA_FILES ${CMAKE_CURRENT_SOURCE_DIR}/swe_ros_lsode_sol.txt
)