# Test directory CMakeLists.txt
# Contains Python-based integration and regression tests

# Find Python for running pytest (prefer the current environment's python3)
find_program(PYTHON3_FROM_PATH python3)
if(PYTHON3_FROM_PATH)
    set(Python3_EXECUTABLE "${PYTHON3_FROM_PATH}" CACHE FILEPATH "Python interpreter" FORCE)
endif()
find_package(Python3 COMPONENTS Interpreter REQUIRED)

execute_process(
    COMMAND ${Python3_EXECUTABLE} -c "import f90nml"
    RESULT_VARIABLE F90NML_STATUS
)
if(NOT F90NML_STATUS EQUAL 0)
    message(FATAL_ERROR "f90nml is required; install it for ${Python3_EXECUTABLE}.")
endif()

# Add subdirectories for different test categories
add_subdirectory(ql-balance)
