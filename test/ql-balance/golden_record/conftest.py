"""Pytest configuration for golden record tests."""

import sys
from pathlib import Path

import pytest

# Add the python package directory to sys.path so imports work
_python_dir = Path(__file__).resolve().parents[3] / "python"
if str(_python_dir) not in sys.path:
    sys.path.insert(0, str(_python_dir))


def pytest_configure(config):
    """Register custom markers."""
    config.addinivalue_line(
        "markers", "golden_record: marks tests as golden record regression tests"
    )


def pytest_collection_modifyitems(config, items):
    """Apply golden_record marker to all tests in this directory."""
    for item in items:
        if "golden_record" in str(item.fspath):
            item.add_marker(pytest.mark.golden_record)
