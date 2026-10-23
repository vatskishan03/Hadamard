"""Certificate of a skew Paley plug-in solution found by qp (negacyclic family, order 8v).
usage: python3 plugcert.py m q "[H_K generators]" "SOL line of qp (four blocks x0|x1|x2, one digit per K-orbit)"
The solution over Z_m x F_q (gcd(m, q) = 1) is mapped to Z_v, v = m q, by the Chinese remainder theorem; its orbit
exponents (in the orbit order of qdrv.shape(v, gens)) go to certify.py, which checks the quaternary and negacyclic
conditions and writes out/cert_<8v>.json; the SHA-256 of the matrix is added with the builder of verify/verify.py."""
import sys, os, json, hashlib, subprocess
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', '..', 'verify'))
from qshape import kshape, cyc_gens
import qdrv
import verify as V

m, q, hK, sol = int(sys.argv[1]), int(sys.argv[2]), json.loads(sys.argv[3]), sys.argv[4]
v = m * q
sh = kshape([m], hK)
blocks = sol.split()
if blocks[0] == 'SOL':
    blocks = blocks[1:]
assert len(blocks) == 4
squares = {x * x % q for x in range(1, q)}
E = []
for b in blocks:
    x0, x1, x2 = b.split('|')
    E.append([int((x0 if g % q == 0 else x1 if g % q in squares else x2)[sh['orb'][g % m]]) for g in range(v)])
gens = cyc_gens(m, q, hK)
d = qdrv.shape(v, gens)
strs = []
for e in E:
    for o in range(d['n']):
        assert len({e[g] for g in range(v) if d['orb'][g] == o}) == 1, 'not invariant under the multiplier group'
    strs.append(''.join(str(e[r]) for r in d['reps']))
print('multiplier generators', gens, 'orbits', d['n'])
print('orbit exponents:', ' '.join(strs))
r = subprocess.run([sys.executable, os.path.join(HERE, 'certify.py'), str(v), json.dumps(gens), ' '.join(strs)],
                   cwd=HERE, capture_output=True, text=True)
print(r.stdout.strip())
if r.returncode:
    sys.exit(r.stderr)
path = os.path.join(HERE, 'out', f'cert_{8 * v}.json')
cert = json.load(open(path))
H, _ = V.BUILDERS[cert['construction']](cert)
h = hashlib.sha256()
for ch in V.csv_chunks(H):
    h.update(ch)
cert['sha256'] = h.hexdigest()
json.dump(cert, open(path, 'w'))
print('sha256', cert['sha256'], '->', path)
