execute_process(COMMAND "${TEST_EXECUTABLE}" reject_bparallel
    RESULT_VARIABLE result OUTPUT_VARIABLE output ERROR_VARIABLE error)
if(result EQUAL 0 OR NOT "${output}${error}" MATCHES
        "periodic Bparallel drive requires the linear compression response")
    message(FATAL_ERROR "Unsupported Bparallel was not rejected correctly: ${output}${error}")
endif()
