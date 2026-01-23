# Golden record regression tests for QL-Balance
#
# These tests compare ql-balance.x output against a golden record generated
# from the main branch. The golden record is automatically created if it
# doesn't exist.
#
# Prerequisites:
#   - Python 3 with pytest, numpy, scipy, h5py
#   - Built ql-balance.x executable
#
# Usage:
#   ctest -R golden_record        # Run golden record tests
#   ctest -R golden_record -V     # Run with verbose output

# Add the golden record test
# This test depends on ql-balance.x being built
add_test(
    NAME ql-balance_golden_record
    COMMAND ${Python3_EXECUTABLE} -m pytest
        ${CMAKE_CURRENT_SOURCE_DIR}/test_golden_record.py
        -v
        --tb=short
        --color=yes
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

# Set test properties
set_tests_properties(ql-balance_golden_record PROPERTIES
    LABELS "golden_record;ql-balance;integration"
    TIMEOUT 600  # 10 minutes timeout (includes potential KiLCA runs)
    ENVIRONMENT "PYTHONPATH=${CMAKE_SOURCE_DIR}/python"
)

# Make the test depend on ql-balance.x being built
# This ensures the executable exists before running the test
set_tests_properties(ql-balance_golden_record PROPERTIES
    FIXTURES_REQUIRED ql-balance_executable
)

# Create a fixture test that ensures ql-balance.x is built
add_test(
    NAME ql-balance_build_check
    COMMAND ${CMAKE_COMMAND} -E true
)
set_tests_properties(ql-balance_build_check PROPERTIES
    FIXTURES_SETUP ql-balance_executable
    DEPENDS ql-balance.x
)
