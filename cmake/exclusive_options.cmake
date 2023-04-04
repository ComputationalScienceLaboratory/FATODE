function(check_mutually_exclusive_options)
  set(on_count 0)

  set(on_options)

  # Loop through all arguments passed to the function
  foreach(option ${ARGN})
    # Check if the option is a boolean and is set to ON
    if(${option})
      # If the option is already set to ON, increment the counter
      math(EXPR on_count "${on_count} + 1")
      list(APPEND on_options ${option})
    endif() 
  endforeach()

  if(on_count GREATER 1)
    message(FATAL_ERROR "The following options are set to on ${on_options}. Only one of ${ARGV} can be set to ON. Please choose one.")
  endif()
  if(on_count EQUAL 0)
    message(FATAL_ERROR "One of ${ARGV} must be set to ON. You have not chosen any options. Please choose one.")
  endif()
endfunction(check_mutually_exclusive_options)
