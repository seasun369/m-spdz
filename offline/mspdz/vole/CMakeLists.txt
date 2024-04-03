#Testing macro
macro (add_test_executable_with_lib _name)
	add_executable(test_${_name} "${_name}.cpp")
	target_link_libraries(test_${_name} ${EMP-OT_LIBRARIES})
endmacro()

macro (add_test_case _name)
	add_test_executable_with_lib(${_name})
  	add_test(NAME ${_name} COMMAND "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_${_name}" WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/")
endmacro()

macro (add_test_case_with_run _name)
	add_test_executable_with_lib(${_name})
	add_test(NAME ${_name} COMMAND "./run" "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_${_name}" WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/")

endmacro()

# Test cases
#add_test_case_with_run(vole_triple)
add_test_case_with_run(vole_triple_new)
# add_test_case_with_run(cope)
# add_test_case_with_run(base_svole)
add_test_case_with_run(base_vole)
#add_test_case_with_run(lpn)
# add_test_case_with_run(vole_f2k_triple)
# add_test_case_with_run(pre)

