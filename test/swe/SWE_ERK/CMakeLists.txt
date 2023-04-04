create_standard_test(
  NAME swe_erk
  SOURCES swe2D_erk_dr.F90 swe2D_upwind.F90
  DATA_FILES ${CMAKE_CURRENT_SOURCE_DIR}/swe_lsode_sol.txt
)