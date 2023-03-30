create_standard_test(
  NAME cbm4_rk 
  SOURCES cbm4_rk_dr.F90 cbm4_parameters.F90
  LIBRARIES csl::FWD LAPACK::LAPACK
)