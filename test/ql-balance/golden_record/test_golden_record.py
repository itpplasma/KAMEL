#!/usr/bin/env python3
"""
Golden record regression tests for QL-Balance.

Runs ql-balance.x with the current branch and compares outputs against golden
record data generated from the main branch.

The golden.h5 file is automatically generated from the main branch if it does
not exist. This ensures tests compare against a known-good reference.
"""

import subprocess
from pathlib import Path

import pytest
from compare import QUANTITIES_TO_COMPARE, compare_hdf5_files, format_comparison_report
from ensure_golden import ensure_golden
from setup_runfolder import get_output_hdf5_path
from setup_runfolder import setup_runfolder as setup_runfolder_external

SCRIPT_DIR = Path(__file__).resolve().parent
GOLDEN_H5 = ensure_golden()


@pytest.fixture(scope="module")
def executable() -> Path:
    """Fixture providing path to ql-balance.x."""
    return SCRIPT_DIR.parents[2] / "build" / "install" / "bin" / "ql-balance.x"


@pytest.fixture(scope="module")
def test_output(executable: Path, tmp_path_factory) -> Path:
    """
    Fixture that runs ql-balance.x and returns path to output HDF5 file.

    This is module-scoped so ql-balance.x only runs once for all tests.
    """
    # Create a unique working directory
    work_dir = tmp_path_factory.mktemp("ql_balance_run")

    # Set up the runfolder
    setup_runfolder_external(work_dir)

    # Run ql-balance.x
    print(f"\nRunning ql-balance.x in {work_dir}...")
    result = subprocess.run(
        [str(executable)],
        cwd=work_dir,
        capture_output=True,
        text=True,
    )

    if result.returncode != 0:
        pytest.fail(
            f"ql-balance.x failed with return code {result.returncode}\n"
            f"stdout:\n{result.stdout}\n"
            f"stderr:\n{result.stderr}"
        )

    return get_output_hdf5_path(work_dir)


def test_golden_record(test_output: Path) -> None:
    """Test that ql-balance.x output matches golden record."""
    # Compare the output against golden record
    results = compare_hdf5_files(str(GOLDEN_H5), str(test_output))

    # Generate report
    report = format_comparison_report(results)
    print(f"\n{report}")

    # Check all comparisons passed
    failures = [r for r in results if not r.passed]
    if failures:
        failure_paths = [r.path for r in failures]
        pytest.fail(
            f"Golden record comparison failed for {len(failures)} quantities:\n"
            f"  {', '.join(failure_paths)}\n\n"
            f"{report}"
        )


@pytest.mark.parametrize(
    "quantity",
    [q.path for q in QUANTITIES_TO_COMPARE],
    ids=[q.path.replace("/", "_") for q in QUANTITIES_TO_COMPARE],
)
def test_golden_record_individual(test_output: Path, quantity: str) -> None:
    """Test individual quantities against golden record (for granular reporting)."""
    # Find the spec for this quantity
    spec = next((q for q in QUANTITIES_TO_COMPARE if q.path == quantity), None)
    if spec is None:
        pytest.skip(f"Quantity {quantity} not in comparison list")

    # Compare just this quantity
    results = compare_hdf5_files(str(GOLDEN_H5), str(test_output), [spec])
    result = results[0]

    if not result.passed:
        msg = f"Comparison failed for {quantity}"
        if result.error_message:
            msg += f": {result.error_message}"
        if result.max_abs_diff is not None:
            msg += f"\n  Max absolute diff: {result.max_abs_diff:.2e}"
        if result.max_rel_diff is not None:
            msg += f"\n  Max relative diff: {result.max_rel_diff:.2e}"
        pytest.fail(msg)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
