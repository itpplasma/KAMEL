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

# Setup fixture: clone + build the main-branch reference and generate
# golden.h5. This is the expensive step (a full from-scratch cmake build of
# main_ref), so it runs as its own fixture test outside the comparison timer.
# ql-balance_golden_record below then only runs the comparison.
add_test(
    NAME ql-balance_golden_setup
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/ensure_golden.py
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)
set_tests_properties(ql-balance_golden_setup PROPERTIES
    LABELS "golden_record;ql-balance;integration"
    TIMEOUT 1800  # building main_ref from scratch can take several minutes
    ENVIRONMENT "PYTHONPATH=${CMAKE_SOURCE_DIR}/python"
    FIXTURES_SETUP ql-balance_golden
)

# Add the golden record test. With golden.h5 prebuilt by the fixture above, the
# test only sets up a runfolder, runs the current-branch ql-balance.x, and
# compares against golden.h5; no main_ref rebuild happens inside this timer.
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

# The comparison needs both the current-branch ql-balance.x and the prebuilt
# golden record.
set_tests_properties(ql-balance_golden_record PROPERTIES
    FIXTURES_REQUIRED "ql-balance_executable;ql-balance_golden"
)

# Create a fixture test that ensures ql-balance.x is built
add_test(
    NAME ql-balance_build_check
    COMMAND ${CMAKE_COMMAND} -E true
)
set_tests_properties(ql-balance_build_check PROPERTIES
    FIXTURES_SETUP ql-balance_executable
)
