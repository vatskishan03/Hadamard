"""Skew-Hadamard matrices of orders 16, 48 and 80 from the exhaustive octic solutions: takes the first SOL line of `./exsmall n 1 0` (n = 4v = 4, 12, 20),
turns its four bit masks into +-1 sequences of length n (negacyclic blocks, A palindromic) and certifies them with the builder of ../../../verify/verify.py
(mkcert.py: out/o<16v>/certificate.json with the SHA-256).  usage: python3 cert_exh.py n         e.g.  python3 cert_exh.py 12   -> order 48
(n = 20 runs the whole exhaustive join of exsmall first: about 7 minutes).  The orders are known; a control of the octic ring, not a new matrix."""
import subprocess, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mkcert import mkcert, bits_to_seq
n = int(sys.argv[1]); v = n // 4
out = subprocess.run(['./exsmall', str(n), '1', '0'], capture_output=True, text=True).stdout
sol = [l for l in out.split('\n') if l.startswith('SOL')][0].split()[1:]
seqs = [bits_to_seq(int(x), n) for x in sol]
d, c = mkcert(seqs, v, f'o{16 * v}')
print('order', c['order'], 'rowsums', c['rowsums'], 'sha256', c['sha256'], '->', d)
