"""Control of zs.py on an independently found solution: the zero-diagonal Szoellosi core of order 850 (hadamard/6796).

usage: python3 control850.py [CERTIFICATE.json]      (default: ../../hadamard/6796/certificate.json of the repository)
G = Z5^2 x Z17 (425 elements), H = <(M, 8)> with M = [[0,1],[1,2]] (order 24; the symmetry group of the certificate), plain array
(T = {1}, target 0 off 0, 4m - 2 = 1698 at 0). The solution of the certificate was found by other code (and is verified by
verify/verify.py); here zs.ZS must agree with it: the orbit structure (the sequences are constant on the orbits of H, A and B on those of
<H, -1> with A(0) = B(0) = 0), the target class values tvec = the PAF class values of the four sequences, the row sums
(0, 40, 7, 7) among st.patterns(), the dynamic-programming counts of sign vectors with those row sums nonzero, and
the sampler (st.sample: sequences with the requested row sum, invariant) and the DP total.
"""
import json, os, sys
import numpy as np
import tbio
from zq import *

path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(tbio.REPO, "hadamard", "6796", "certificate.json")
core = json.load(open(path))["core"]
SIGN = {"+": 1, "-": -1, "0": 0}
A, B, C, D = (np.array([SIGN[ch] for ch in x], dtype=np.int64) for x in core["sequences"])
M = mat([[0, 1], [1, 2]])
Hm = block_diag(M, mat([[8]]))
assert Hm.tolist() == core["symmetry_generators"][0]
st = ZS([5, 5, 17], [Hm], plain=True, name="850core")
print(st.info())
ok = True
def rep(name, cond):
    global ok
    print("  %-70s %s" % (name, "OK" if cond else "FAIL")); ok &= bool(cond)
rep("A, B, C, D constant on the H-orbits", all(np.array_equal(x, x[st.perms[0]]) for x in (A, B, C, D)))
rep("A, B symmetric, A(0) = B(0) = 0", all(np.array_equal(x, x[st.neg]) and x[0] == 0 for x in (A, B)) and (np.abs(C) == 1).all() and (np.abs(D) == 1).all())
tot = sum(st.paf_classes(x[None, :].astype(np.int8))[0] for x in (A, B, C, D))
rep("sum of the PAF class values == st.tvec (0 off 0, 1698 at 0)", np.array_equal(tot, st.tvec) and st.tvec[st.zero_class] == 1698)
rows = tuple(int(x.sum()) for x in (A, B, C, D))
rep("row sums %s in st.patterns()" % (rows,), rows in [tuple(p) for p in st.patterns()])
for kind, x, r in (("L", A, rows[0]), ("L", B, rows[1]), ("A", C, rows[2]), ("A", D, rows[3])):
    rep("DP count of block %s with row sum %d is positive (%g)" % (kind, r, st.counts(kind).get(r, 0)), st.counts(kind).get(r, 0) > 0)
rng = np.random.default_rng(3)
for kind, r in (("L", 0), ("L", 40), ("A", 7)):
    S = st.sample(kind, r, 500, rng); X = st.seq(kind, S)
    rep("sampler %s row sum %d: 500 sequences with that sum, H-invariant" % (kind, r), (X.sum(1) == r).all() and (X[:, st.perms[0]] == X).all())
orb, sz, free = st.kind_orbits("L")
tc = st.counts("L")
rep("DP counts of kind L sum to 2^(free orbits) = 2^%d" % len(free), sum(tc.values()) == 2.0 ** len(free))
print("control850:", "ALL OK" if ok else "FAILURES")
sys.exit(0 if ok else 1)
