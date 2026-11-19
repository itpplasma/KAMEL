#!/usr/bin/env python3
# Generic numeric comparator for codes without a code-specific one.
# Usage: gr_numcompare.py <dir_A> <dir_B> [rtol] [atol] [near_zero_floor]
# Recursively compares every regular file present in BOTH trees: whitespace-
# separated numeric columns where parseable (HDF5 files: all numeric datasets,
# except run-time metadata), else byte-equality. Reports max
# relative diff per file. Exit 1 if any file diverges beyond rtol.
#
# near_zero_floor (default 1e-9): an element where BOTH the reference and the
# current value are at or below this magnitude carries no physical signal and is
# treated as equal. Oscillatory fields (e.g. KIM's charge density rho) are
# O(1e-12) floating-point noise at their zero crossings, ~14 orders below their
# O(1e2) physical scale; a *relative* diff there is meaningless and flips on the
# reference build's last bits (a fortnum-vs-GSL last-bit wobble reads as tens of
# percent). This floors only genuine near-zero noise -- the smallest physically
# meaningful value in the current golden set is ~1e-6, far above it -- so it
# never masks a real divergence on a physical value. Tune per code via arg 5.
import os
import sys

A, B = sys.argv[1], sys.argv[2]
rtol = float(sys.argv[3]) if len(sys.argv) > 3 else 1e-7
atol = float(sys.argv[4]) if len(sys.argv) > 4 else 1e-12
floor = float(sys.argv[5]) if len(sys.argv) > 5 else 1e-9

# zone_*_poy_test_err.dat is the Poynting energy-balance residual (the solution's
# own numerical self-consistency error, |div S - (P_abs + jE)|/max(1,|div S|)). It
# is a diagnostic, never consumed downstream, and is a near-cancellation: its
# pointwise value amplifies last-bit differences ~1e8x, so two independent
# backends that agree bit-for-bit on every physical field still diverge tens of
# percent here (itpplasma-KAMEL#164/#172). It therefore cannot be compared
# relatively ref-vs-cur. Instead each build is checked against an absolute
# self-consistency bar (the healthy residual peaks at ~3e-2 at the plasma edge, so
# the default 1e-1 catches gross solve breakage without flagging last-bit noise).
POY_BAR = float(os.environ.get("GR_POY_TEST_ERR_BAR", "1e-1"))

SKIP = {"run.log", "exit_code.txt", "runtime_seconds.txt", "migrate.log", "prepare.log"}


# Run-dependent HDF5 datasets (wall-clock and thread count) carry no physics.
H5_IGNORED = {"runtime", "omp_num_threads"}
H5_SEP = "::"


def is_h5(path):
    return path.endswith((".h5", ".hdf5"))


def h5_dataset_names(p):
    """Names of the physics datasets in an HDF5 file; None if the file cannot be opened.

    A missing h5py raises ImportError (loud failure) instead of dropping the file.
    """
    import h5py

    names = []

    def visit(name, obj):
        if isinstance(obj, h5py.Dataset) and name not in H5_IGNORED:
            names.append(name)

    try:
        with h5py.File(p, "r") as f:
            f.visititems(visit)
    except OSError:
        return None
    return names


def h5_dataset_nums(p, name):
    """Flat float list of one dataset; compound (real, imag) datasets flatten to floats."""
    import h5py
    import numpy as np

    with h5py.File(p, "r") as f:
        obj = f[name]
        data = obj[()]
        if obj.dtype.names:
            out = []
            for field in obj.dtype.names:
                out.extend(np.ravel(data[field]).astype(float).tolist())
            return out
        if obj.dtype.kind in "iufb":
            return np.ravel(data).astype(float).tolist()
        return None  # strings etc.: compared as raw bytes of the value


def nums(p):
    out = []
    try:
        for line in open(p, errors="ignore"):
            for tok in line.split():
                try:
                    out.append(float(tok))
                except ValueError:
                    pass
    except Exception:
        return None
    return out


def last_col(p):
    """Residual column (last whitespace token per line) of a save_real_array file."""
    out = []
    try:
        for line in open(p, errors="ignore"):
            toks = line.split()
            if toks:
                try:
                    out.append(float(toks[-1]))
                except ValueError:
                    pass
    except Exception:
        return None
    return out


def rel_files(root):
    found = []
    for dirpath, _dirs, files in os.walk(root):
        for f in files:
            if f in SKIP:
                continue
            full = os.path.join(dirpath, f)
            if os.path.isfile(full):
                rel = os.path.relpath(full, root)
                if is_h5(f):
                    # one virtual entry per dataset: schema drift is reported per name
                    names = h5_dataset_names(full)
                    if names is None:  # unreadable: whole-file byte comparison
                        found.append(rel)
                    else:
                        found.extend(rel + H5_SEP + n for n in names)
                else:
                    found.append(rel)
    return set(found)


fail = 0
checked = 0
if os.path.isdir(A) and os.path.isdir(B):
    files_a = rel_files(A)
    files_b = rel_files(B)
    for rel in sorted(files_a - files_b):
        checked += 1
        fail += 1
        print(f"{rel}: MISSING(cur)")
    for rel in sorted(files_b - files_a):
        checked += 1
        fail += 1
        print(f"{rel}: MISSING(ref)")
    common = sorted(files_a & files_b)
else:
    common = []

for rel in common:
    pa, pb = os.path.join(A, rel), os.path.join(B, rel)
    if H5_SEP in rel:
        fa, _, dsname = pa.partition(H5_SEP)
        fb = pb.partition(H5_SEP)[0]
        na, nb = h5_dataset_nums(fa, dsname), h5_dataset_nums(fb, dsname)
        if na is None or nb is None:
            import h5py

            with h5py.File(fa, "r") as ha, h5py.File(fb, "r") as hb:
                same = ha[dsname][()].tobytes() == hb[dsname][()].tobytes()
            checked += 1
            print(f"{rel}: {'MATCH(value)' if same else 'DIFFER(value)'}")
            fail += 0 if same else 1
            continue
        if len(na) != len(nb):
            checked += 1
            fail += 1
            print(f"{rel}: DIFFER(shape) ref={len(na)} cur={len(nb)}")
            continue
        if not na:
            checked += 1
            print(f"{rel}: MATCH(empty)")
            continue
    if os.path.basename(rel).endswith("poy_test_err.dat"):
        ra, rb = last_col(pa) or [], last_col(pb) or []
        ma = max((abs(v) for v in ra), default=0.0)
        mb = max((abs(v) for v in rb), default=0.0)
        checked += 1
        ok = ma <= POY_BAR and mb <= POY_BAR
        print(
            f"{rel}: poy_test_err self-consistency max(ref)={ma:.3e} "
            f"max(cur)={mb:.3e} bar={POY_BAR:.1e} {'PASS' if ok else 'FAIL'}"
        )
        fail += 0 if ok else 1
        continue
    if H5_SEP not in rel:
        na, nb = nums(pa), nums(pb)
    if H5_SEP not in rel and (na is None or nb is None or len(na) != len(nb) or not na):
        same = open(pa, "rb").read() == open(pb, "rb").read()
        checked += 1
        print(f"{rel}: {'MATCH(bytes)' if same else 'DIFFER(bytes/shape)'}")
        fail += 0 if same else 1
        continue
    mx = 0.0
    skipped = 0
    for x, y in zip(na, nb):
        if abs(x) <= floor and abs(y) <= floor:
            skipped += 1  # near-zero noise: carries no signal
            continue
        d = abs(x - y)
        s = d / (abs(y) + atol)
        mx = max(mx, 0.0 if d <= atol else s)
    checked += 1
    ok = mx <= rtol
    note = f" ({skipped} near-zero skipped)" if skipped else ""
    print(f"{rel}: max_rel={mx:.3e} {'PASS' if ok else 'FAIL'}{note}")
    fail += 0 if ok else 1

print(f"-- {checked} files compared, {fail} over rtol={rtol:.1e}")
sys.exit(1 if fail else 0)
