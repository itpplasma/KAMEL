"""Pytest configuration for golden record tests."""

import pytest


def pytest_configure(config):
    """Register custom markers."""
    config.addinivalue_line(
        "markers", "golden_record: marks tests as golden record regression tests"
    )


def pytest_collection_modifyitems(config, items):
    """Apply golden_record marker to all tests in this directory."""
    for item in items:
        if "golden_record_balance" in str(item.fspath):
            item.add_marker(pytest.mark.golden_record)
