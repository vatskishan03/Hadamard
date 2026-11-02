#!/usr/bin/env python3
"""check_known.py w h WORKDIR CERT.json [--swap-reduce] : is the pair (a, b) of an existing certificate (cyclic group Z_2w) among the
group images of the solutions found by the engine in WORKDIR?  (a, b must be invariant under the multiplier group H of order h.)"""
import sys, os, json, subprocess
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pgplib import *
w, h, workdir, cert = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3], sys.argv[4]
sr = "--swap-reduce" in sys.argv
types = parse_types(w, h, sr); P = types["P"]; m = P.m
pgp = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src", "pgp")
out = subprocess.run([pgp, "join", str(w), str(h), workdir] + (["--swap-reduce"] if sr else []), capture_output=True, text=True, check=True).stdout
found = set()
for ln in out.splitlines():
    f = [int(x) for x in ln.split()[:12]]
    A, B, C, D = (mask_to_orb(f[2 * i], f[2 * i + 1], m) for i in range(4))
    found.add(sol_str(assemble(P, A, B, C, D, f[8], f[9], f[10], f[11])))
c = json.load(open(cert))
a, b = c["a"], c["b"]
assert len(a) == 2 * w
def layers(seq):
    """cyclic Z_2w sequence -> [inf, orbits] per layer (CRT index y <-> (y mod 2, y mod w)); also check H-invariance"""
    L = []
    for e in (0, 1):
        v = [None] * (m + 1)
        for y in range(2 * w):
            if y % 2 != e: continue
            r = y % w; slot = 0 if r == 0 else 1 + P.cls[r]
            if v[slot] is None: v[slot] = seq[y]
            elif v[slot] != seq[y]: raise SystemExit("sequence is not invariant under H of order %d" % h)
        L.append(v)
    return L
a0, a1 = layers(a); b0, b1 = layers(b)
tgt = (a0, a1, b0, b1)
assert pgp_check(P, tgt)
imgs = gen_group_images(tgt, m)
hit = [k for k in imgs if sol_str(imgs[k]) in found]
print(f"known pair of {cert}: invariant under H (order {h}), exact PAF check ok; engine found {len(found)} representatives; "
      f"{len(hit)} of its {len(imgs)} images are among them -> {'FOUND' if hit else 'NOT FOUND'}")
sys.exit(0 if hit else 1)
