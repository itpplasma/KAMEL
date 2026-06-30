#!/usr/bin/env python3
# Generic numeric comparator for codes without a code-specific one.
# Usage: gr_numcompare.py <dir_A> <dir_B> [rtol] [atol]
# Recursively compares every regular file present in BOTH trees: whitespace-
# separated numeric columns where parseable, else byte-equality. Reports max
# relative diff per file. Exit 1 if any file diverges beyond rtol.
import sys, os

A, B = sys.argv[1], sys.argv[2]
rtol = float(sys.argv[3]) if len(sys.argv) > 3 else 1e-7
atol = float(sys.argv[4]) if len(sys.argv) > 4 else 1e-12

SKIP = {"run.log", "exit_code.txt", "runtime_seconds.txt", "migrate.log"}


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
    for x, y in zip(na, nb):
        d = abs(x - y)
        s = d / (abs(y) + atol)
        mx = max(mx, 0.0 if d <= atol else s)
    checked += 1
    ok = mx <= rtol
    print(f"{rel}: max_rel={mx:.3e} {'PASS' if ok else 'FAIL'}")
    fail += 0 if ok else 1

print(f"-- {checked} files compared, {fail} over rtol={rtol:.1e}")
sys.exit(1 if fail else 0)
