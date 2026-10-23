"""cert_oct.py v gens "A B C D"   (hex orbit codes from oskew): expand, check at the U level, build the +-1 sequences by psi, certify with the builder of verify/verify.py
(mkcert.py: out/o<16v>/certificate.json).  usage: python3 cert_oct.py 61 "[12]" "09c27 0c325 2ed30 259e2"     -> order 976, all checks True (a control: 976 is known)"""
import sys, json, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from octic import *
from dcfeas import orbit_data
from mkcert import mkcert
def certify(v, gens, sol, name=None, verbose=True):
    H, orb, reps, sz, neg = orbit_data(v, gens)
    Us = [expand_codes(s, orb) for s in sol]
    # U-level checks
    tot = None
    for U in Us:
        P = gmul(U, gadj(U)); tot = P if tot is None else tot + P
    ok_U = bool(np.all(tot[1:] == 0)) and np.array_equal(tot[0], [16 * v, 0, 0, 0])
    A = Us[0]; skewU = bool(np.array_equal(A + gadj(A), np.array([[2, 0, 0, 0]] + [[0, 0, 0, 0]] * (v - 1))))
    seqs = [psi(U, v) for U in Us]
    assert all(set(np.unique(s)) <= {-1, 1} for s in seqs), "psi image not +-1"
    n = 4 * v; a = seqs[0]
    pal = a[0] == 1 and all(a[j] == a[n - j] for j in range(1, n))
    # round trip
    rt = all(np.array_equal(phi(s, v), U) for s, U in zip(seqs, Us))
    negpaf = [sum(nega_paf(s, k) for s in seqs) for k in range(n)]
    ok_seq = negpaf[0] == 4 * n and all(x == 0 for x in negpaf[1:])
    if verbose: print(f'v={v}: U-level sum UU* = 16v delta: {ok_U}; A skew in ring: {skewU}; A palindromic: {pal}; phi(psi(U)) = U: {rt}; sequence NPAF sum = 4n delta: {ok_seq}')
    d, c = mkcert([list(map(int, s)) for s in seqs], v, name or f'o{16*v}')
    return d, c, (ok_U and skewU and pal and rt and ok_seq)
if __name__ == "__main__":
    v = int(sys.argv[1]); gens = json.loads(sys.argv[2]); sol = sys.argv[3].split()
    d, c, ok = certify(v, gens, sol); print(d, c['order'], c['rowsums'], ok)
