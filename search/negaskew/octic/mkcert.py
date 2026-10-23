"""mkcert.py: writes a 'goethals-seidel-negacyclic' certificate of order 4n = 16 v (four +-1 sequences of length n = 4v, negaperiodic autocorrelations
summing to 4n delta, the first palindromic) with the builder of ../../../verify/verify.py, which checks it, into out/<name>/certificate.json and
sequences.txt (certificate with SHA-256).  module: mkcert(seqs (four lists of +-1), v, name) -> (directory, certificate); bits_to_seq(a, n)."""
import json, sys, hashlib, os
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', '..', 'verify'))
import verify as V
def mkcert(seqs, v, name=None, family='out'):
    n = len(seqs[0]); order = 4 * n
    cert = dict(order=order, construction="goethals-seidel-negacyclic", skew=True,
                group=f"negacyclic Z[x]/(x^{n}+1) = Z[zeta_8][Z_{v}]", v=v, n=n,
                rowsums=[int(sum(s)) for s in seqs],
                sequences=["".join('+' if x > 0 else '-' for x in s) for s in seqs])
    H, _ = V.build_goethals_seidel_negacyclic(cert)      # raises ValueError unless the sequences are complementary, A palindromic
    sha = hashlib.sha256()
    for ch in V.csv_chunks(H): sha.update(ch)
    cert['sha256'] = sha.hexdigest()
    d = f'{family}/{name or order}'
    os.makedirs(d, exist_ok=True)
    json.dump(cert, open(d + '/certificate.json', 'w'), indent=1)
    open(d + '/sequences.txt', 'w').write("\n".join(cert['sequences']) + "\n")
    return d, cert
def bits_to_seq(a, n): return [-1 if (a >> i) & 1 else 1 for i in range(n)]
