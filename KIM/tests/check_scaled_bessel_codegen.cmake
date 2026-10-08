# Returning to the fortnum_special umbrella import must fail this check on
# affected GNU Fortran builds. The existing Fourier-kernel test covers values.
set(helper "__flr2_fourier_kernel_m_MOD_scaled_bessel_pair")
execute_process(
    COMMAND "${OBJDUMP}" -dr "--disassemble=${helper}" "${KIM_LIBRARY}"
    RESULT_VARIABLE result
    OUTPUT_VARIABLE disassembly
    ERROR_VARIABLE error
    TIMEOUT 30)
if(NOT result EQUAL 0)
    message(FATAL_ERROR
        "Cannot inspect scaled_bessel_pair with ${OBJDUMP} (${result}): ${error}")
endif()
if(NOT disassembly MATCHES "<${helper}>:")
    message(FATAL_ERROR "scaled_bessel_pair was not found in ${KIM_LIBRARY}")
endif()
if(disassembly MATCHES "_gfortran_ieee_procedure_(entry|exit)")
    message(FATAL_ERROR
        "scaled_bessel_pair contains IEEE environment save/restore overhead; "
        "import bessel_in directly from fortnum_special_bessel")
endif()
