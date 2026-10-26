"""Complex Yamada-Kiyasu step of the n = 955 construction: K478 (quaternary conference matrix of order 478 = (955 + 1)/2) ->
quaternary Hadamard matrix of order 1910 -> Turyn doubling -> Hadamard matrix of order 3820 = 4 * 955.

usage: python3 cyk955.py [--in K478.npy] [--out DIR] [--csv-gz FILE]
Reads K478.npy (assemble_K.py) from --out (default: current directory). Uses the repository's verify/verify.py for the Yamada
circulants (q = 953, r = 3, xi = (1, 2)) and for the canonical CSV. Checks exactly: C C* = 477 I over the integers; the Yamada
circulants M^2 + N^2 = q I; K K* = 1910 I (float32 products of small integers, exact); H H^T = 3820 I; prints the SHA-256 of the
int8 matrix and of its canonical CSV (the bytes of numpy.savetxt(H, fmt="%d", delimiter=","), what verify/check_matrix reads); saves
OUT/H3820.npy; --csv-gz writes the CSV gzip-compressed (mtime 0, so the file is reproducible too).
"""
import sys, os, json, hashlib, gzip, argparse, numpy as np, importlib.util
sys.dont_write_bytecode = True
from scipy.linalg import circulant
import tbio
ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
ap.add_argument("--in", dest="kin", default=None, help="K478.npy (default: OUT/K478.npy)")
ap.add_argument("--out", default=None, help="output directory (default: current directory)")
ap.add_argument("--csv-gz", default=None, help="also write the canonical CSV, gzip-compressed, to this file")
args = ap.parse_args()
out = tbio.outdir(args.out)
spec = importlib.util.spec_from_file_location("hv", os.path.join(tbio.REPO, "verify", "verify.py"))
hv = importlib.util.module_from_spec(spec); sys.modules["hv"] = hv
spec.loader.exec_module(hv)
K = np.load(args.kin or os.path.join(out, 'K478.npy'))
Cr = np.rint(K.real).astype(np.int64); Ci = np.rint(K.imag).astype(np.int64)
assert np.allclose(Cr + 1j*Ci, K)
c = Cr.shape[0]
# exact integer check C C* = (c-1) I
re = Cr@Cr.T + Ci@Ci.T; im = Ci@Cr.T - Cr@Ci.T
print("core order", c, "C C* = (c-1) I exactly:", bool((re == (c-1)*np.eye(c, dtype=np.int64)).all() and (im == 0).all()))
q = 953
assert hv._is_prime(q) and q % 8 == 1
n = (q+1)//2
assert c == n+1
# find r non-residue and primitive xi
psi = lambda a: pow(a % q, (q-1)//2, q)
r = next(a for a in range(2, q) if psi(a) == q-1)
xi = None
for xa in range(1, 30):
    for xb in range(1, 30):
        try:
            alpha, beta = hv._yamada_alpha_beta(q, r, (xa, xb)); xi = (xa, xb); break
        except ValueError:
            pass
    if xi: break
print("q", q, "r", r, "xi", xi)
s = np.arange(n)
assert np.array_equal(alpha, alpha[-s % n]) and np.array_equal(beta, beta[-s % n])
assert beta[0] == 0 and (beta[1:] != 0).all() and (alpha != 0).all()
paf = np.array([alpha @ np.roll(alpha, -k) + beta @ np.roll(beta, -k) for k in range(n)])
assert paf[0] == q and (paf[1:] == 0).all()
print("Yamada circulants: M^2+N^2 = qI OK")
zero = (Cr == 0) & (Ci == 0)
assert zero.sum(0).tolist() == [1]*c and zero.sum(1).tolist() == [1]*c
E = np.where(Cr == 1, 0, np.where(Ci == 1, 1, np.where(Cr == -1, 2, np.where(Ci == -1, 3, -1))))
E = E[:, np.argmax(E < 0, axis=1)]
nz = E >= 0
E = np.where(nz, (E - np.r_[0, E[0, 1:]][None, :]) % 4, -1)
E = np.where(nz, (E - np.r_[0, E[1:, 0]][:, None]) % 4, -1)
assert E[0, 0] == -1 and (E[0, 1:] == 0).all() and (E[1:, 0] == 0).all() and (np.diag(E) == -1).all()
Z = E[1:, 1:]
Zr = np.where(Z == 0, 1, np.where(Z == 2, -1, 0)).astype(np.int64)
Zi = np.where(Z == 1, 1, np.where(Z == 3, -1, 0)).astype(np.int64)
print("Z row/col sums zero:", not (Zr.sum(0).any() or Zr.sum(1).any() or Zi.sum(0).any() or Zi.sum(1).any()))
I = np.eye(n, dtype=np.int64); J2 = np.ones((2, 2), dtype=np.int64); P = np.array([[1, -1], [-1, 1]], dtype=np.int64)
A = np.kron(P, circulant(beta)) + np.kron(J2, I)
B = np.kron(P, circulant(alpha))
e, z = np.ones((1, 2*n), dtype=np.int64), np.zeros((1, 2*n), dtype=np.int64)
one, nil = np.ones((1, 1), dtype=np.int64), np.zeros((1, 1), dtype=np.int64)
Xr = np.block([[one, e], [-e.T, (A + B)//2]]); Xi = np.block([[nil, z], [z.T, (A - B)//2]])
Yr = np.block([[one, z], [z.T, np.kron(J2, Zr) + np.kron(P, I)]]); Yi = np.block([[nil, e], [-e.T, np.kron(J2, Zi)]])
Kr = np.block([[Xr, Yr], [Yr.T, -Xr.T]]).astype(np.int8); Ki = np.block([[Xi, Yi], [-Yi.T, Xi.T]]).astype(np.int8)
assert (np.abs(Kr.astype(np.int16)) + np.abs(Ki.astype(np.int16)) == 1).all()
N2 = Kr.shape[0]
print("quaternary matrix order", N2, "(= 2(q+2) =", 2*(q+2), ")")
Fr, Fi = Kr.astype(np.float32), Ki.astype(np.float32)
ok1 = bool(np.array_equal(Fr@Fr.T + Fi@Fi.T, N2*np.eye(N2, dtype=np.float32))) and bool(np.array_equal(Fi@Fr.T - Fr@Fi.T, np.zeros((N2, N2), dtype=np.float32)))
print("K K* = 1910 I exactly:", ok1)
Xh, Yh = Kr + Ki, Kr - Ki
H = np.block([[Xh, Yh], [-Yh, Xh]]).astype(np.int8)
Hf = H.astype(np.float32)
G = Hf@Hf.T
ok2 = bool(np.array_equal(G, (H.shape[0])*np.eye(H.shape[0], dtype=np.float32)))
print("H order", H.shape[0], "entries +-1:", bool(np.isin(H, (-1, 1)).all()), " H H^T = 3820 I exactly:", ok2)
np.save(os.path.join(out, 'H3820.npy'), H)
print("sha256(H int8 rowmajor):", hashlib.sha256(H.tobytes()).hexdigest())
csv = hashlib.sha256()
for chunk in hv.csv_chunks(H): csv.update(chunk)
print("sha256(canonical CSV):", csv.hexdigest())
if args.csv_gz:
    with open(args.csv_gz, "wb") as raw, gzip.GzipFile(fileobj=raw, mode="wb", filename="", mtime=0, compresslevel=9) as gz:
        for chunk in hv.csv_chunks(H): gz.write(chunk)
    print("wrote", args.csv_gz)
