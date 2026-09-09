execute_process(COMMAND "${TEST_EXECUTABLE}" "${CASE}"
    RESULT_VARIABLE result OUTPUT_VARIABLE output ERROR_VARIABLE error)
if(result EQUAL 0 OR NOT "${output}${error}" MATCHES
        "electrostatic_periodic: resonance not found")
    message(FATAL_ERROR "Nonresonant mode was not rejected correctly: ${output}${error}")
endif()
