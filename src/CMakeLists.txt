set(KIM_mods plas_parameter_mod.f90)
add_library(KIM_lib "${KIM_mods}")

add_executable(KIM_exe kim_main.f90)
set_target_properties(KIM_exe PROPERTIES 
                        OUTPUT_NAME KIM_exe
                        RUNTIME_OUTPUT_DIRECTORY exe)
target_link_libraries(KIM_exe PRIVATE KIM_lib)