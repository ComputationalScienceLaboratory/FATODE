create_standard_test(
  NAME swe_rk
  SOURCES swe2D_rk_dr.F90 swe2D_upwind.F90
  DATA_FILES ${CMAKE_CURRENT_SOURCE_DIR}/swe_rk_lsode_sol.txt
)