#!/usr/bin/env python3
"""
HDF5 comparison utilities for golden record tests.

This module provides functions to compare specific quantities in HDF5 files
with configurable absolute and relative tolerances.
"""

from dataclasses import dataclass
from typing import Any

import h5py
import numpy as np


@dataclass
class ComparisonResult:
    """Result of a single quantity comparison."""
    path: str
    passed: bool
    max_abs_diff: float | None = None
    max_rel_diff: float | None = None
    golden_shape: tuple | None = None
    actual_shape: tuple | None = None
    error_message: str | None = None


@dataclass
class QuantitySpec:
    """Specification for a quantity to compare."""
    path: str  # HDF5 path, e.g., "/group/dataset"
    rtol: float = 1e-8  # Relative tolerance
    atol: float = 1e-15  # Absolute tolerance


def get_nested_item(h5file: h5py.File, path: str) -> Any:
    """Get an item from HDF5 file by path string."""
    parts = path.strip("/").split("/")
    obj = h5file
    for part in parts:
        obj = obj[part]  # type: ignore
    return obj


def compare_arrays(
    golden: np.ndarray,
    actual: np.ndarray,
    rtol: float,
    atol: float,
) -> tuple[bool, float | None, float | None]:
    """
    Compare two numpy arrays with tolerances.

    Returns (passed, max_abs_diff, max_rel_diff).
    """
    if golden.shape != actual.shape:
        return False, None, None

    # Compute absolute difference
    abs_diff = np.abs(golden - actual)
    max_abs_diff = float(np.max(abs_diff))

    # Compute relative difference (avoiding division by zero)
    with np.errstate(divide='ignore', invalid='ignore'):
        rel_diff = np.where(
            np.abs(golden) > 0,
            abs_diff / np.abs(golden),
            np.where(abs_diff > 0, np.inf, 0.0)
        )
    max_rel_diff = float(np.max(rel_diff[np.isfinite(rel_diff)])) if np.any(np.isfinite(rel_diff)) else 0.0

    # Check if arrays are close
    passed = np.allclose(golden, actual, rtol=rtol, atol=atol)

    return passed, max_abs_diff, max_rel_diff


def compare_quantity(
    golden_file: h5py.File,
    actual_file: h5py.File,
    spec: QuantitySpec,
) -> ComparisonResult:
    """Compare a single quantity between golden and actual HDF5 files."""
    try:
        golden_data = get_nested_item(golden_file, spec.path)
        actual_data = get_nested_item(actual_file, spec.path)
    except KeyError as e:
        return ComparisonResult(
            path=spec.path,
            passed=False,
            error_message=f"Path not found: {e}",
        )

    # Handle datasets (arrays)
    if isinstance(golden_data, h5py.Dataset):
        golden_arr = golden_data[:]
        actual_arr = actual_data[:]

        if golden_arr.shape != actual_arr.shape:
            return ComparisonResult(
                path=spec.path,
                passed=False,
                golden_shape=golden_arr.shape,
                actual_shape=actual_arr.shape,
                error_message=f"Shape mismatch: golden {golden_arr.shape} vs actual {actual_arr.shape}",
            )

        passed, max_abs, max_rel = compare_arrays(
            golden_arr, actual_arr, spec.rtol, spec.atol
        )

        return ComparisonResult(
            path=spec.path,
            passed=passed,
            max_abs_diff=max_abs,
            max_rel_diff=max_rel,
            golden_shape=golden_arr.shape,
            actual_shape=actual_arr.shape,
            error_message=None if passed else f"Values differ beyond tolerance (rtol={spec.rtol}, atol={spec.atol})",
        )

    # Handle scalar attributes or other types
    return ComparisonResult(
        path=spec.path,
        passed=False,
        error_message=f"Unsupported data type: {type(golden_data)}",
    )


def compare_hdf5_files(
    golden_path: str,
    actual_path: str,
    quantities: list[QuantitySpec],
) -> list[ComparisonResult]:
    """
    Compare specific quantities between two HDF5 files.

    Args:
        golden_path: Path to the golden reference HDF5 file.
        actual_path: Path to the actual (test) HDF5 file.
        quantities: List of quantity specifications to compare.

    Returns:
        List of comparison results, one per quantity.
    """
    results = []

    with h5py.File(golden_path, "r") as golden_file, \
         h5py.File(actual_path, "r") as actual_file:
        for spec in quantities:
            result = compare_quantity(golden_file, actual_file, spec)
            results.append(result)

    return results


def format_comparison_report(results: list[ComparisonResult]) -> str:
    """Format comparison results as a human-readable report."""
    lines = []
    passed_count = sum(1 for r in results if r.passed)
    total_count = len(results)

    lines.append(f"Comparison Results: {passed_count}/{total_count} passed")
    lines.append("=" * 60)

    for result in results:
        status = "PASS" if result.passed else "FAIL"
        lines.append(f"\n[{status}] {result.path}")

        if result.golden_shape is not None:
            lines.append(f"  Shape: {result.golden_shape}")

        if result.max_abs_diff is not None:
            lines.append(f"  Max absolute diff: {result.max_abs_diff:.2e}")

        if result.max_rel_diff is not None:
            lines.append(f"  Max relative diff: {result.max_rel_diff:.2e}")

        if result.error_message and not result.passed:
            lines.append(f"  Error: {result.error_message}")

    return "\n".join(lines)


if __name__ == "__main__":
    # Example usage / quick test
    import sys

    if len(sys.argv) != 3:
        print("Usage: python compare_hdf5.py <golden.h5> <actual.h5>")
        sys.exit(1)

    # Import quantities from the config
    from quantities_to_compare import QUANTITIES_TO_COMPARE

    results = compare_hdf5_files(sys.argv[1], sys.argv[2], QUANTITIES_TO_COMPARE)
    print(format_comparison_report(results))

    # Exit with error if any comparison failed
    if not all(r.passed for r in results):
        sys.exit(1)
