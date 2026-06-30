#!/usr/bin/env python3
# Generic numeric comparator for codes without a code-specific one.
# Usage: gr_numcompare.py <dir_A> <dir_B> [rtol] [atol]
# Compares every regular file present in BOTH dirs: whitespace-separated numeric
# columns where parseable, else byte-equality. Reports max relative diff per file.
import sys, os, math
A, B = sys.argv[1], sys.argv[2]
rtol = float(sys.argv[3]) if len(sys.argv) > 3 else 1e-7
atol = float(sys.argv[4]) if len(sys.argv) > 4 else 1e-12
def nums(p):
    out = []
    try:
        for line in open(p, errors="ignore"):
            for tok in line.split():
                try: out.append(float(tok))
                except ValueError: pass
    except Exception: return None
    return out
fail = 0; checked = 0
common = sorted(set(os.listdir(A)) & set(os.listdir(B))) if os.path.isdir(A) and os.path.isdir(B) else []
for f in common:
    pa, pb = os.path.join(A, f), os.path.join(B, f)
    if not (os.path.isfile(pa) and os.path.isfile(pb)): continue
    if f in ("run.log", "exit_code.txt", "runtime_seconds.txt"): continue
    na, nb = nums(pa), nums(pb)
    if na is None or nb is None or len(na) != len(nb) or not na:
        same = open(pa, 'rb').read() == open(pb, 'rb').read()
        checked += 1; print(f"{f}: {'MATCH(bytes)' if same else 'DIFFER(bytes/shape)'}");  fail += 0 if same else 1
        continue
    mx = 0.0
    for x, y in zip(na, nb):
        d = abs(x - y); s = d / (abs(y) + atol)
        mx = max(mx, 0.0 if d <= atol else s)
    checked += 1
    ok = mx <= rtol
    print(f"{f}: max_rel={mx:.3e} {'PASS' if ok else 'FAIL'}")
    fail += 0 if ok else 1
print(f"-- {checked} files compared, {fail} over rtol={rtol:.1e}")
sys.exit(1 if fail else 0)
