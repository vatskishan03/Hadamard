"""Certificate writer + verification by the repository's own code (verify/verify.py of this repository, imported in place).
Border conventions (verify.build_bordered_gs):  H = [[K, R], [C, M]],  R[x,(i,g)] = Rhat[x,i,coset(g)],  C[(i,g),x] = Chat[i,coset(g),x]
with Rhat = E (m x 4k, column index i*k + coset) and Chat = F (4k x m, row index i*k + coset)."""
import sys, os, json, hashlib
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.normpath(os.path.join(HERE, '..', '..', '..'))          # repository root (search/bgschain/e2e -> root)
sys.path.insert(0, os.path.join(REPO, 'verify'))
import verify as V

def sstr(row):
    return ''.join('+' if x > 0 else '-' for x in row)

def make_certificate(rec, K, E, F, t_in, t_out, order=None):
    """rec: engine solution record (group, subgroup_generators, gens, sequences, ...); K (m x m), E (m x 4k), F (4k x m): +-1 arrays"""
    mods = rec['group']; seqs = rec['sequences']
    cert = dict(order=int(rec['order'] if order is None else order), construction='bordered-goethals-seidel', group=[int(x) for x in mods],
                subgroup_generators=rec['subgroup_generators'], symmetries=rec['gens'], autocorrelation=dict(subgroup=int(t_in), elsewhere=int(t_out)),
                rowsums=[int(sum(1 if c == '+' else -1 for c in s)) for s in seqs], sequences=seqs,
                border=dict(X=[sstr(r) for r in K], Rhat=[sstr(r) for r in E], Chat=[sstr(r) for r in F]))
    return cert

def verify_with_repo(cert, out_path):
    """Build the matrix with the repository's own builder (verify.BUILDERS), record the SHA-256 of its canonical CSV, write the certificate to the
    file `out_path` and re-check it with verify.build_checked (order, entries +-1, H H^T = nI, SHA-256).  Returns (order, digest, certificate)."""
    H, seqs = V.BUILDERS[cert['construction']](cert)
    sha = hashlib.sha256()
    for ch in V.csv_chunks(H): sha.update(ch)
    cert = dict(cert); cert['sha256'] = sha.hexdigest()
    d = os.path.dirname(os.path.abspath(out_path))
    os.makedirs(d, exist_ok=True)
    with open(out_path, 'w') as f: json.dump(cert, f)
    H, seqs, digest = V.build_checked(json.load(open(out_path)))         # raises ValueError on any failure
    return H.shape[0], digest, cert
