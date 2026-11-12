#!/usr/bin/env python3
"""Unit tests for the QL-Balance HDF5 comparator quantity list."""

import compare
import h5py
import numpy as np


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


def test_time_series_quantities_are_listed():
    paths = _paths()
    for q in compare._TIME_SERIES_QUANTITIES:
        assert f"/f_6_2/{q}" in paths
    for q in (
        "br_abs_res",
        "br_formfactor",
        "timstep_evol.dat",
        "T_tot_phi_e",
        "br_abs_antenna_factor",
    ):
        assert f"/f_6_2/{q}" in paths


def _write(path, series, drop=None):
    with h5py.File(path, "w") as f:
        for q in compare._TIME_SERIES_QUANTITIES:
            if q == drop:
                continue
            data = series[q]
            f.create_dataset(f"/f_6_2/{q}", data=data)


def _series(n=5):
    rng = np.random.default_rng(0)
    out = {q: rng.random(n) + 1.0 for q in compare._TIME_SERIES_QUANTITIES}
    out["timstep_evol.dat"] = rng.random((3, n)) + 1.0
    out["br_formfactor"] = out["br_formfactor"] + 1j * rng.random(n)
    return out


def _run(tmp_path, golden, actual, drop=None):
    g, a = str(tmp_path / "g.h5"), str(tmp_path / "a.h5")
    _write(g, golden)
    _write(a, actual, drop=drop)
    specs = [compare.QuantitySpec(f"/f_6_2/{q}") for q in compare._TIME_SERIES_QUANTITIES]
    return {r.path: r for r in compare.compare_hdf5_files(g, a, specs)}


def test_equal_time_series_pass(tmp_path):
    s = _series()
    results = _run(tmp_path, s, {k: v.copy() for k, v in s.items()})
    assert len(results) == len(compare._TIME_SERIES_QUANTITIES)
    assert all(r.passed for r in results.values())


def test_perturbed_time_series_fails(tmp_path):
    s = _series()
    bad = {k: v.copy() for k, v in s.items()}
    bad["br_abs_res"][2] *= 1 + 1e-6
    bad["br_formfactor"][1] += 1e-3j
    results = _run(tmp_path, s, bad)
    assert not results["/f_6_2/br_abs_res"].passed
    assert not results["/f_6_2/br_formfactor"].passed
    assert results["/f_6_2/T_tot_phi_e"].passed


def test_shape_mismatch_fails(tmp_path):
    s = _series(5)
    shorter = _series(4)
    results = _run(tmp_path, s, shorter)
    assert not any(r.passed for r in results.values())
    assert all("Shape mismatch" in r.error_message for r in results.values())


def test_missing_dataset_fails(tmp_path):
    s = _series()
    results = _run(tmp_path, s, {k: v.copy() for k, v in s.items()}, drop="br_formfactor")
    assert not results["/f_6_2/br_formfactor"].passed
    assert "Path not found" in results["/f_6_2/br_formfactor"].error_message
    assert results["/f_6_2/br_abs_res"].passed


def test_compound_complex_dataset_compared_by_field(tmp_path):
    dt = np.dtype([("r", "f8"), ("i", "f8")])
    g, a = str(tmp_path / "g.h5"), str(tmp_path / "a.h5")
    base = np.array([(1.0, 2.0), (3.0, 4.0)], dtype=dt)
    off = base.copy()
    off["i"][1] += 1e-3
    for p, d in ((g, base), (a, off)):
        with h5py.File(p, "w") as f:
            f.create_dataset("/f_6_2/br_formfactor", data=d)
    spec = compare.QuantitySpec("/f_6_2/br_formfactor")
    res = compare.compare_hdf5_files(g, a, [spec])[0]
    assert not res.passed
    res = compare.compare_hdf5_files(g, g, [spec])[0]
    assert res.passed
