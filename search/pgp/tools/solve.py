#!/usr/bin/env python3
"""solve.py w h WORKDIR [--swap-reduce] [--oracle FILE] [--threads T]
 stage 1 (canonical keys) + stage 2 (join) + exact checks; with --oracle: the union of the group images of the found solutions must
 equal the oracle's complete list of pairs."""
import sys, subprocess, os, shutil, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pgplib import *
ap = argparse.ArgumentParser()
ap.add_argument("w", type=int); ap.add_argument("h", type=int); ap.add_argument("workdir")
ap.add_argument("--swap-reduce", action="store_true"); ap.add_argument("--oracle"); ap.add_argument("--threads", default="1")
ap.add_argument("--pgp", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src", "pgp"))
ap.add_argument("--keep", action="store_true"); ap.add_argument("--join-only", action="store_true"); ap.add_argument("--no-rev", action="store_true"); ap.add_argument("--closure", action="store_true")
a = ap.parse_args()
w, h = a.w, a.h
types = parse_types(w, h, a.swap_reduce, a.pgp)
P = types["P"]; m = P.m
if not a.keep and os.path.exists(a.workdir): shutil.rmtree(a.workdir)
extra = (["--swap-reduce"] if a.swap_reduce else []) + (["--no-rev"] if a.no_rev else [])
if not a.join_only:
    subprocess.run([a.pgp, "rows", str(w), str(h), a.workdir, "--threads", a.threads, "--quiet"] + extra, check=True)
out = subprocess.run([a.pgp, "join", str(w), str(h), a.workdir] + (["--swap-reduce"] if a.swap_reduce else []), capture_output=True, text=True, check=True)
sols = set(); nlines = 0
for ln in out.stdout.splitlines():
    f = ln.split()
    ap_, an, bp, bn, cp, cn, dp, dn = [int(x) for x in f[:8]]
    ainf, binf, cinf, dinf = [int(x) for x in f[8:12]]
    A, B, C, D = (mask_to_orb(x, y, m) for x, y in ((ap_, an), (bp, bn), (cp, cn), (dp, dn)))
    sol = assemble(P, A, B, C, D, ainf, binf, cinf, dinf)
    assert pgp_check(P, sol), "join produced a non-Golay pair"
    sols.add(sol_str(sol)); nlines += 1
print(f"(w,h)=({w},{h}) m={m}: {nlines} join matches, {len(sols)} distinct solutions as found, all verified exactly on Z_{2*w}")
# classes (connected components of the found representatives under the symmetry group) and, optionally, the full group closure
def tosol(s): return tuple([{'+': 1, '-': -1}[c] for c in x] for x in s)
unseen = set(sols); classes = []
while unseen:
    s = next(iter(unseen)); imgs = gen_group_images(tosol(s), m)
    members = [k for k in (sol_str(v) for v in imgs.values()) if k in unseen]
    unseen -= set(members); classes.append((min(sol_str(v) for v in imgs.values()), len(imgs), len(members)))
import hashlib
chash = hashlib.sha1(repr(sorted(c[0] for c in classes)).encode()).hexdigest()[:12]
print(f"  ordered pairs (a, b) in these classes (sum of orbit sizes under the symmetry group): {sum(c[1] for c in classes)}")
print(f"  equivalence classes among the found representatives: {len(classes)} (class-set hash {chash})"
      + (f"   (class sizes under the full group: {sorted(set(c[1] for c in classes))})" if classes else ""))
if classes:
    print("  class representatives (a0 a1 b0 b1; infinity point first, then orbits 0..m-1):")
    for c in classes[:6]: print("   ", " ".join(c[0]))
if a.closure or a.oracle:
    allimg = {}
    for s in sols:
        for k, v in gen_group_images(tosol(s), m).items(): allimg[sol_str(v)] = v
    print(f"  union of group images: {len(allimg)} ordered pairs")
    for s, v in allimg.items(): assert pgp_check(P, v)
    if a.oracle:
        ora = set(tuple(ln.split()) for ln in open(a.oracle))
        ok = set(allimg) == ora
        print(f"  oracle: {len(ora)} ordered pairs; equal={ok}; only_engine={len(set(allimg) - ora)} only_oracle={len(ora - set(allimg))}")
        sys.exit(0 if ok else 1)
