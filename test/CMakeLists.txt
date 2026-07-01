# Test directory CMakeLists.txt
#
# Golden-record regression tests for KIM, QL-Balance, and KiLCA used to live here
# (test/ql-balance/golden_record) and ran under CTest. They have moved to
# test/golden/ and now run ONLY in the dedicated "golden-record" GitHub Actions
# job (.github/workflows/golden.yml) -- never in local ctest / make test, so a
# normal build stays fast and offline.
#
# No CTest-registered tests remain in this directory.
