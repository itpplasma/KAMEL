#!/bin/bash
# Compare ref vs cur outputs at the code's golden tolerances (RTOL 1e-7, ATOL 1e-12).
# Honors GOLDEN_RECORD_SKIP_CASES. Usage: gr_compare.sh [cases...]
set -u
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; cd "$ROOT"
. "$ROOT/config.sh"
exec "$ROOT/bin/$GR_COMPARE" "$ROOT/out/ref" "$ROOT/out/cur" "$@"
