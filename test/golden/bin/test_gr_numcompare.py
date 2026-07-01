#!/usr/bin/env python3
"""Unit tests for the recursive numeric comparator."""
import os, subprocess, sys, struct
from pathlib import Path

HERE = Path(__file__).resolve().parent
CMP = HERE / "gr_numcompare.py"


def _run(a, b, *extra):
    r = subprocess.run([sys.executable, str(CMP), str(a), str(b), *extra],
                       capture_output=True, text=True)
    return r.returncode, r.stdout


def _write(p, text):
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(text)


def test_nested_numeric_match(tmp_path):
    for root in ("ref", "cur"):
        _write(tmp_path / root / "linear-data" / "m_6_n_2" / "EB.dat",
               "1.0 2.0\n3.0 4.0\n")
        _write(tmp_path / root / "background-data" / "b0z.dat", "5.0\n6.0\n")
    rc, out = _run(tmp_path / "ref", tmp_path / "cur")
    assert rc == 0, out
    assert "EB.dat" in out and "b0z.dat" in out          # nested files were found
    assert "FAIL" not in out


def test_nested_numeric_divergence_flags(tmp_path):
    _write(tmp_path / "ref" / "sub" / "EB.dat", "1.0\n")
    _write(tmp_path / "cur" / "sub" / "EB.dat", "1.1\n")   # 10% off >> rtol
    rc, out = _run(tmp_path / "ref", tmp_path / "cur")
    assert rc == 1, out
    assert "FAIL" in out


def test_uff_byte_identical(tmp_path):
    payload = struct.pack("<3d", 1.0, 2.0, 3.0)
    for root in ("ref", "cur"):
        p = tmp_path / root / "formfactors.flre.uff"
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_bytes(payload)
    rc, out = _run(tmp_path / "ref", tmp_path / "cur")
    assert rc == 0, out
    assert "formfactors.flre.uff" in out and "MATCH(bytes)" in out


def test_uff_byte_divergence_flags(tmp_path):
    for root, val in (("ref", 1.0), ("cur", 9.0)):
        p = tmp_path / root / "formfactors.flre.uff"
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_bytes(struct.pack("<d", val))
    rc, out = _run(tmp_path / "ref", tmp_path / "cur")
    assert rc == 1, out
    assert "DIFFER(bytes" in out


def test_near_zero_noise_not_flagged(tmp_path):
    # An oscillatory field (like KIM's rho): physical elements match, but the
    # zero-crossing elements are float noise (~1e-12) that differ hugely in
    # *relative* terms between two builds. Those must not fail the comparison.
    # Element 3 diverges by 1.6e-12 > atol (1e-12) on a ~4e-12 value: 33% rel.
    # The pre-existing atol floor does NOT catch it (d > atol); the near-zero
    # floor must. This is the exact KIM rho / itpplasma-KAMEL#164 scenario.
    _write(tmp_path / "ref" / "fields" / "rho.dat",
           "9.5780e1\n3.1738e-6\n3.9e-12\n")
    _write(tmp_path / "cur" / "fields" / "rho.dat",
           "9.5780e1\n3.1738e-6\n5.5e-12\n")
    rc, out = _run(tmp_path / "ref", tmp_path / "cur")
    assert rc == 0, out
    assert "FAIL" not in out
    assert "near-zero skipped" in out


def test_near_zero_floor_does_not_mask_physical_divergence(tmp_path):
    # A value above the floor that genuinely diverges must still fail, even when
    # the file also contains skippable near-zero noise.
    _write(tmp_path / "ref" / "rho.dat", "9.50e1\n-6.03e-14\n")
    _write(tmp_path / "cur" / "rho.dat", "9.60e1\n 3.04e-13\n")   # 1% on the big element
    rc, out = _run(tmp_path / "ref", tmp_path / "cur")
    assert rc == 1, out
    assert "FAIL" in out


def test_volatile_files_skipped(tmp_path):
    for root in ("ref", "cur"):
        _write(tmp_path / root / "run.log", f"started at {root}\n")
        _write(tmp_path / root / "EB.dat", "1.0\n")
    rc, out = _run(tmp_path / "ref", tmp_path / "cur")
    assert rc == 0, out
    assert "run.log" not in out          # volatile logs excluded even when nested
