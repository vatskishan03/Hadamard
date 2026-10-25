#!/usr/bin/env python3
"""tablegen.py v ulen [REM_EVALS=3e9] -> measures, for the best row-sum pattern(s) of (v,|U|), the flat fraction p(Phi) (REM sampling)
and the SA generation rate g0(Phi) with flatbig; writes tg_<v>_<ulen>.json.  All runs single-threaded."""
import sys, os, json, math, subprocess, re
from math import comb, log10
import numpy as np
import gsflat
v, ulen = int(sys.argv[1]), int(sys.argv[2]); evals = float(sys.argv[3]) if len(sys.argv) > 3 else 3e9
here = os.path.dirname(os.path.abspath(__file__))
S = gsflat.build(v, ulen); m, C = S['m'], S['C']
tab = f"tab_{v}.txt"
if not os.path.exists(tab):
    with open(tab, "w") as f:
        f.write(f"{v} {m} {C} {ulen} {S['pm']}\n"); f.write(" ".join(repr(float(x)) for x in S['eta'].real) + "\n"); f.write(" ".join(repr(float(x)) for x in S['eta'].imag) + "\n")
opts = gsflat.rowsum_options(v, ulen)
pats = gsflat.patterns(v, ulen)
if not pats: print("no pattern"); json.dump(dict(v=v, ulen=ulen, patterns=[]), open(f"tg_{v}_{ulen}.json", "w")); sys.exit()
def role(r):
    f0, Ssum, npos = opts[r][0]
    return dict(r=r, f0=f0, npos=npos, logN=log10(comb(m, npos)))
scored = sorted(pats, key=lambda q: -min(role(r)['logN'] for r in q))
best = scored[0]; roles = [role(r) for r in best]
print(v, ulen, "m", m, "C", C, "patterns", len(pats), "best", best, [round(x['logN'], 2) for x in roles], flush=True)
binary = "./flatbig32" if C <= 32 else "./flatbig64"
rr = min(roles, key=lambda x: x['r'])      # measure on the smallest row sum
nbase = evals / (rr['npos'] * (m - rr['npos']))
out = subprocess.run([binary, "rem", tab, str(rr['f0']), str(rr['npos']), str(nbase), "7"], capture_output=True, text=True)
rem = []
for line in out.stdout.splitlines():
    mm = re.match(r"Phi<([0-9.]+)\s+(\d+)\s+p=([0-9.e+-]+)", line)
    if mm: rem.append((float(mm.group(1)), int(mm.group(2)), float(mm.group(3))))
print(out.stderr.strip(), flush=True)
gen = {}
for thr in (0.14, 0.18, 0.24):
    fn = f"w/_tg_{v}.bin"
    if os.path.exists(fn): os.remove(fn)
    o = subprocess.run([binary, "gen2", tab, str(rr['f0']), str(rr['npos']), str(thr), "0.4", "10000", "4000", "5", fn], capture_output=True, text=True)
    nprop, nd, dt = o.stdout.split(); gen[thr] = float(nd) / float(dt)
    print("  SA thr", thr, "distinct/s", round(gen[thr], 1), flush=True)
json.dump(dict(v=v, ulen=ulen, m=m, C=C, pattern=best, roles=roles, rem=rem, rem_evals=evals, gen=gen, patterns_all=[list(q) for q in pats]), open(f"tg_{v}_{ulen}.json", "w"))
