"""Exact Z_2 compression test of the bordered array (search/bsz.py), faster than bsz.z2_ok: every block is compressed
to the sums (a, b) of its even and odd slice of Z_m (an orbit of H lies in one slice), a + b = r, and the four values
(a - b)^2 must add up to 4N. The two slices are summed separately and the four blocks are paired, so the cost is
that of the subset sums instead of their product."""
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import bsz


def z2_ok_fast(st, rows):
    """same result as bsz.z2_ok"""
    opts = []
    for kind, r in zip(bsz.KINDS, rows):
        sz, orb = st.sizes(kind), st.orb(kind)
        par = [int(np.flatnonzero(orb == k)[0]) % 2 for k in range(len(sz))]
        S = [{0}, {0}]
        for s, p in zip(sz, par):
            S[p] = {a + e * s for a in S[p] for e in (1, -1)}
        opts.append({(2 * a - r) ** 2 for a in S[0] if (r - a) in S[1]})
    s01 = {x + y for x in opts[0] for y in opts[1]}
    s23 = {x + y for x in opts[2] for y in opts[3]}
    return any(4 * st.N - x in s23 for x in s01)
