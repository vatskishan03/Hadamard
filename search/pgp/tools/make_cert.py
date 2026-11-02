#!/usr/bin/env python3.11
"""make_cert.py -- from a periodic Golay pair (a, b) of length n = (q - 1)/2 to the symmetric Hadamard matrix of order 4q and a
certificate accepted by verify/verify.py (construction "symmetric-involution", Theorem F of symmetric/PERIODIC_GOLAY.md).

    make_cert.py --q 389 --pair JSON [--class 0] --out DIR          pair from extract_pairs.py output ("classes": [{"a": "+-..", "b": ".."}])
    make_cert.py --q 389 --a STRING --b STRING --out DIR           the two +- strings of length 194 (cyclic group Z_194 = Z_2 x Z_97, CRT index)
    make_cert.py --q 101 --a S --b S --group 2x5^2 --t 25 --out DIR   other groups, as in the repository's certificates
    make_cert.py --q 149 --from-cert symmetric/596/certificate.json --out DIR --compare symmetric/596/certificate.json

It never writes into the repository: DIR gets symmetric/<4q>/certificate.json (the family directory name "symmetric" is required by
verify.py).  Steps, all checked: (1) (a, b) are +-1 and PAF_a + PAF_b = 2n delta, exactly and again by FFT; t has order 2 in the group;
(2) the certificate is completed by the repository's builder (H and its SHA-256), which itself re-checks every identity of PERIODIC_GOLAY.md;
(3) H = H^T, entries +-1, H H^T = 4q I are re-checked here independently (float64 matrix product, exact for these sizes);
(4) the official verifier is run on the new directory."""
import argparse, hashlib, json, os, subprocess, sys
import numpy as np

REPO = os.environ.get("HADAMARD_REPO") or os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", ".."))
sys.path.insert(0, os.path.join(REPO, "verify"))
import verify  # noqa: E402  (read-only use of the repository's verifier / builder)

ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
ap.add_argument("--q", type=int, required=True)
ap.add_argument("--pair"); ap.add_argument("--class", dest="cls", type=int, default=0)
ap.add_argument("--a"); ap.add_argument("--b"); ap.add_argument("--from-cert")
ap.add_argument("--group"); ap.add_argument("--t", type=int)
ap.add_argument("--out", required=True); ap.add_argument("--compare")
args = ap.parse_args()

q = args.q
n = (q - 1) // 2
if args.pair:
    c = json.load(open(args.pair))["classes"][args.cls]
    a = [1 if ch == "+" else -1 for ch in c["a"]]; b = [1 if ch == "+" else -1 for ch in c["b"]]
elif args.from_cert:
    c = json.load(open(args.from_cert)); a, b = c["a"], c["b"]
    args.group = args.group or c.get("group"); args.t = args.t or c.get("t")
else:
    a = [1 if ch == "+" else -1 for ch in args.a]; b = [1 if ch == "+" else -1 for ch in args.b]
group = args.group or str(n)
t = args.t if args.t is not None else n // 2
assert len(a) == len(b) == n, f"pair length {len(a)} != (q-1)/2 = {n}"

# (1) the pair: exact PAF on the group, and FFT (cyclic groups only)
tg, add, neg = verify.abelian_group(group)
assert tg == n and add[t, t] == 0 and t != 0, "t must be an element of order 2"
aa, bb = np.array(a), np.array(b)
paf = np.array([int(aa @ aa[add[:, s]] + bb @ bb[add[:, s]]) for s in range(n)])
assert paf[0] == 2 * n and not paf[1:].any(), "not a periodic Golay pair"
if "x" not in group and "^" not in group:
    assert np.allclose(np.abs(np.fft.fft(aa)) ** 2 + np.abs(np.fft.fft(bb)) ** 2, 2 * n, atol=1e-6)
print(f"[1] (a, b) is a periodic Golay pair of length {n} over group {group} (row sums {aa.sum()}, {bb.sum()}; {aa.sum()**2 + bb.sum()**2} = 2n)")

# (2) certificate through the repository's builder
cert = {"order": 4 * q, "construction": "symmetric-involution", "symmetric": True, "q": q, "group": group, "t": t,
        "a": a, "b": b, "sha256": ""}
H, _ = verify.BUILDERS[cert["construction"]](cert)
sha = hashlib.sha256()
for chunk in verify.csv_chunks(H):
    sha.update(chunk)
cert["sha256"] = sha.hexdigest()
verify.build_checked(cert)
print(f"[2] builder: H of order {H.shape[0]}, sha256 {cert['sha256']}")

# (3) independent checks of H
Hf = H.astype(np.float64)
assert H.shape == (4 * q, 4 * q) and set(np.unique(H)) <= {-1, 1}, "entries must be +-1"
assert np.array_equal(H, H.T), "H not symmetric"
G = Hf @ Hf.T
assert np.array_equal(G, 4 * q * np.eye(4 * q)), "H H^T != 4q I"
print(f"[3] independent check: H = H^T, entries +-1, H H^T = {4 * q} I")

# (4) write and run the official verifier
d = os.path.abspath(os.path.join(args.out, "symmetric", str(4 * q)))
os.makedirs(d, exist_ok=True)
with open(os.path.join(d, "certificate.json"), "w") as f:
    json.dump(cert, f, indent=1); f.write("\n")
r = subprocess.run([sys.executable, os.path.join(REPO, "verify", "verify.py"), d], capture_output=True, text=True)
print("[4] verify/verify.py:", (r.stdout + r.stderr).strip())
assert r.returncode == 0, "verify.py failed"
if args.compare:
    ref = json.load(open(args.compare))
    same = ref["sha256"] == cert["sha256"]
    print(f"[5] compare with {args.compare}: sha256 {'identical' if same else 'DIFFERENT'}; "
          f"a,b {'identical' if (ref['a'], ref['b']) == (cert['a'], cert['b']) else 'different'}; group/t {'identical' if (ref.get('group'), ref['t']) == (cert['group'], cert['t']) else 'different'}")
    assert same or (ref["a"], ref["b"]) != (cert["a"], cert["b"]), "same pair but different matrix"
print("certificate written to", os.path.join(d, "certificate.json"))
