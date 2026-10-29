"""odd_spence_cert.py -- independent exact check of a certificate of the construction "spence" over a group G = Z_m1 x ... x Z_mk of odd order v
with pairwise coprime moduli (README, "Order 5912"; shares no code with verify/verify.py or the engines).
(1) The four sequences (index np.ravel_multi_index of the coordinates) have row sums (+-1)^4, and the sum of their periodic autocorrelations is
4v at 0 and -4 elsewhere (integer arithmetic, every shift).
(2) G is mapped to the cyclic group Z_v by the Chinese remainder theorem and H = [[-S, S x 1^T], [S x 1, M]] is built with circulant blocks, M the
Goethals-Seidel array of the four sequences and S that of their row sums; H is +-1 and H H^T = (4v + 4) I (float32 products of small integers: exact).
usage: odd_spence_cert.py certificate.json"""
import sys, json
import numpy as np


def need(cond, msg):
    """a check that python -O cannot remove"""
    if not cond: sys.exit('FAIL: ' + msg)


def paf(x, mods):
    """exact periodic autocorrelation over G of x (flat array in C order)"""
    X = x.reshape(mods); axes = tuple(range(len(mods)))
    return np.array([int((X * np.roll(X, [-c for c in g], axes)).sum()) for g in np.ndindex(*mods)])


def crt(mods):
    """cyclic index in Z_v of every element of G (C order)"""
    v = int(np.prod(mods)); coords = np.array(np.unravel_index(np.arange(v), mods))
    return sum(c * (v // m) * pow(v // m, -1, m) for c, m in zip(coords, mods)) % v


def build(seqs):
    """H of order 4v + 4 from four +-1 sequences on Z_v with row sums +-1"""
    v = len(seqs[0]); g = np.arange(v)
    A, B, C, D = [np.asarray(s, np.float32)[(g[None, :] - g[:, None]) % v] for s in seqs]      # X[g, h] = x(h - g)
    R = np.zeros((v, v), np.float32); R[g, -g % v] = 1
    a, b, c, d = (int(s.sum()) for s in seqs)
    M = np.block([[A, B @ R, C @ R, D @ R], [-B @ R, A, D.T @ R, -C.T @ R], [-C @ R, -D.T @ R, A, B.T @ R], [-D @ R, C.T @ R, -B.T @ R, A]])
    S = np.array([[a, b, c, d], [-b, a, d, -c], [-c, -d, a, b], [-d, c, -b, a]], np.float32); one = np.ones((v, 1), np.float32)
    return np.vstack([np.hstack([-S, np.kron(S, one.T)]), np.hstack([np.kron(S, one), M])])


def main(path):
    cert = json.load(open(path)); mods = cert['group']; v = int(np.prod(mods))
    need(cert['construction'] == 'spence' and v % 2 == 1 and all(np.gcd(a, b) == 1 for i, a in enumerate(mods) for b in mods[i + 1:]), 'odd order, coprime moduli')
    need(not set(''.join(cert['sequences'])) - {'+', '-'}, 'sequences must consist of + and - only')
    seqs = [np.array([1 if ch == '+' else -1 for ch in s], np.int64) for s in cert['sequences']]
    need(len(seqs) == 4 and all(len(s) == v for s in seqs) and cert['order'] == 4 * v + 4, 'four sequences of length v, order 4v + 4')
    rs = [int(s.sum()) for s in seqs]; print('order', cert['order'], 'group', mods, 'row sums', rs)
    need(all(abs(r) == 1 for r in rs), 'row sums must be +-1')
    tot = sum(paf(s, mods) for s in seqs)
    need(tot[0] == 4 * v and (tot[1:] == -4).all(), 'sum PAF is not 4v at 0 and -4 elsewhere')
    print(f'sum PAF = {4 * v} at 0 and -4 at the other {v - 1} elements: OK')
    p = crt(mods); cyc = []
    for s in seqs: t = np.zeros(v, np.int64); t[p] = s; cyc.append(t)
    H = build(cyc); n = len(H)
    need(n == cert['order'] and np.isin(H, (-1, 1)).all() and np.array_equal(H @ H.T, n * np.eye(n, dtype=np.float32)), 'H H^T is not nI')
    print(f'circulant build: H is +-1 of order {n}, H H^T = {n} I: OK')


if __name__ == '__main__': main(sys.argv[1])
