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
    common = sorted(rel_files(A) & rel_files(B))
else:
    common = []

for rel in common:
    pa, pb = os.path.join(A, rel), os.path.join(B, rel)
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
