#!/usr/bin/env python3
# Generic numeric comparator for codes without a code-specific one.
# Usage: gr_numcompare.py <dir_A> <dir_B> [rtol] [atol] [near_zero_floor]
# Recursively compares every regular file present in BOTH trees: whitespace-
# separated numeric columns where parseable, else byte-equality. Reports max
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
import sys, os

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
                found.append(os.path.relpath(full, root))
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
    if os.path.basename(rel).endswith("poy_test_err.dat"):
        ra, rb = last_col(pa) or [], last_col(pb) or []
        ma = max((abs(v) for v in ra), default=0.0)
        mb = max((abs(v) for v in rb), default=0.0)
        checked += 1
        ok = ma <= POY_BAR and mb <= POY_BAR
        print(f"{rel}: poy_test_err self-consistency max(ref)={ma:.3e} "
              f"max(cur)={mb:.3e} bar={POY_BAR:.1e} {'PASS' if ok else 'FAIL'}")
        fail += 0 if ok else 1
        continue
    na, nb = nums(pa), nums(pb)
    if na is None or nb is None or len(na) != len(nb) or not na:
        same = open(pa, "rb").read() == open(pb, "rb").read()
        checked += 1
        print(f"{rel}: {'MATCH(bytes)' if same else 'DIFFER(bytes/shape)'}")
        fail += 0 if same else 1
        continue
    mx = 0.0
    skipped = 0
    for x, y in zip(na, nb):
        if abs(x) <= floor and abs(y) <= floor:
            skipped += 1                       # near-zero noise: carries no signal
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
