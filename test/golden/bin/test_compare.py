#!/usr/bin/env python3
"""Unit tests for the QL-Balance HDF5 comparator quantity list."""
import compare


def _paths():
    return {spec.path for spec in compare.QUANTITIES_TO_COMPARE}


def test_linear_profiles_include_all_time_steps():
    paths = _paths()
    for t in range(9):
        for q in compare._LINEAR_PROFILE_QUANTITIES:
            assert f"/f_6_2/LinearProfiles/{t}/{q}" in paths


def test_kin_profiles_include_initial_and_final_time_steps():
    paths = _paths()
    for t in (1000, 1008):
        for q in ("Te", "Ti", "n", "Er"):
            assert f"/f_6_2/KinProfiles/{t}/{q}" in paths
