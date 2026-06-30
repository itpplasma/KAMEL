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


def test_volatile_files_skipped(tmp_path):
    for root in ("ref", "cur"):
        _write(tmp_path / root / "run.log", f"started at {root}\n")
        _write(tmp_path / root / "EB.dat", "1.0\n")
    rc, out = _run(tmp_path / "ref", tmp_path / "cur")
    assert rc == 0, out
    assert "run.log" not in out          # volatile logs excluded even when nested
