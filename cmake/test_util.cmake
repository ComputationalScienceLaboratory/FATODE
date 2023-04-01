################################################################################
# build and add a standard test (one linked to the fatode library)

function(create_standard_test)
  set(prefix TEST)
  set(singleValues NAME WORKING_DIRECTORY)
  set(multiValues SOURCES LIBRARIES DATA_FILES)
  include(CMakeParseArguments)
  cmake_parse_arguments(${prefix} " " "${singleValues}" "${multiValues}" ${ARGN})
  add_executable(test_${TEST_NAME} ${TEST_SOURCES})

  target_link_libraries(test_${TEST_NAME} PUBLIC csl::FATODE)

  # link additional libraries
  foreach(library ${TEST_LIBRARIES})
    target_link_libraries(test_${TEST_NAME} PUBLIC ${library})
  endforeach()


  if(NOT DEFINED TEST_WORKING_DIRECTORY)
    set(TEST_WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/test_binaries")
  endif()

  # Place the test executables in their own directory
  set_target_properties(test_${TEST_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${TEST_WORKING_DIRECTORY})

  add_fatode_test(${TEST_NAME} test_${TEST_NAME} "" ${TEST_WORKING_DIRECTORY} ${TEST_DATA_FILES})
endfunction(create_standard_test)

################################################################################
# Add a test

function(add_fatode_test test_name test_binary test_args working_dir)
  add_test(NAME ${test_name}
            COMMAND ${test_binary} ${test_args}
            WORKING_DIRECTORY ${working_dir})
  set_tests_properties(${test_name} PROPERTIES TIMEOUT 60)  # seconds

  # copy data files to test working directory
  set(data_files ${ARGN})
  foreach(data_file ${data_files})
    add_custom_command(TARGET ${test_binary} POST_BUILD
                       COMMAND ${CMAKE_COMMAND} -E copy
                       ${data_file}
                       ${working_dir}/${data_file})
  endforeach()
endfunction(add_fatode_test)

################################################################################
