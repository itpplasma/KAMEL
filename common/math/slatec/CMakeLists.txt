file(GLOB libslatec_files "*.f")
add_library(KIM_slatec STATIC ${libslatec_files})
# Provide legacy name 'slatec' via output name
set_target_properties(KIM_slatec PROPERTIES OUTPUT_NAME slatec)