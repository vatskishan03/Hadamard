#!/usr/bin/env python3.11
"""extract_pairs.py w h WORKDIR [--swap-reduce] : print one representative of every equivalence class as +-1 strings on the cyclic group
Z_2w (CRT index y <-> (y mod 2, y mod w)); the PAF check here uses numpy FFT (independent of the exact integer check in pgplib)."""
import sys, os, subprocess, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pgplib import *
w, h, workdir = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
sr = "--swap-reduce" in sys.argv
types = parse_types(w, h, sr); P = types["P"]; m = P.m
pgp = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src", "pgp")
out = subprocess.run([pgp, "join", str(w), str(h), workdir] + (["--swap-reduce"] if sr else []), capture_output=True, text=True, check=True).stdout
sols = {}
for ln in out.splitlines():
    f = [int(x) for x in ln.split()[:12]]
    A, B, C, D = (mask_to_orb(f[2 * i], f[2 * i + 1], m) for i in range(4))
    s = assemble(P, A, B, C, D, f[8], f[9], f[10], f[11]); sols[sol_str(s)] = s
def tosol(s): return tuple([{'+': 1, '-': -1}[c] for c in x] for x in s)
unseen = set(sols); reps = []
while unseen:
    s = next(iter(unseen)); imgs = gen_group_images(tosol(s), m)
    keys = set(sol_str(v) for v in imgs.values()); unseen -= keys
    reps.append(min(keys))
def cyc(l0, l1):
    return [ (l1 if y % 2 else l0)[0 if y % w == 0 else 1 + P.cls[y % w]] for y in range(2 * w)]
res = []
for r in sorted(reps):
    a0, a1, b0, b1 = tosol(r)
    a, b = cyc(a0, a1), cyc(b0, b1)
    fa, fb = np.abs(np.fft.fft(a)) ** 2, np.abs(np.fft.fft(b)) ** 2
    ok = bool(np.allclose(fa + fb, 4 * w, atol=1e-6))           # PSD_a + PSD_b = 2 * (2w)  <=> periodic Golay
    ok2 = pgp_check(P, (a0, a1, b0, b1))
    res.append(dict(a="".join("+" if x > 0 else "-" for x in a), b="".join("+" if x > 0 else "-" for x in b), fft_ok=ok, exact_ok=ok2, sum_a=int(sum(a)), sum_b=int(sum(b))))
json.dump(dict(w=w, h=h, length=2 * w, group="Z_2 x Z_%d ~ Z_%d (CRT)" % (w, 2 * w), classes=res), sys.stdout, indent=1)
print()
