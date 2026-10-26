"""k = 8 border (K 28x28, E 28x32, F 32x28) for the spike-at-0 pattern y_i = a_i + d_i delta_0 on Q = Z2^3 (any odd/even a, d with |d|^2 = 28),
from the Z7-Singer border search z7search (z7search.c in this directory, compiled into a temporary directory on demand; the search is deterministic
and takes about 4 s per d).  Identities checked exactly:
  E E^T = 32 I,  F F^T = 32 I - 4 (I4 x J8),  K K^T = 28 I,  Mc E^T + F K^T = 0,  Mc = GS(a) x J8 + GS(d) x I8  (the compressed GS array)."""
import os, subprocess, tempfile
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
def gs4(d):
    a1, a2, a3, a4 = d
    return np.array([[a1, a2, a3, a4], [-a2, a1, a4, -a3], [-a3, -a4, a1, a2], [-a4, a3, -a2, a1]], dtype=np.int64)
def z7search_rows(d):
    """run z7search for the amplitudes d (first solution of the exhaustive search); returns the 8 rows (4 generators e, 4 generators f) of 32 entries"""
    with tempfile.TemporaryDirectory(prefix='z7search_') as tmp:
        exe = os.path.join(tmp, 'z7search')
        subprocess.run(['gcc', '-O2', os.path.join(HERE, 'z7search.c'), '-o', exe], check=True)
        out = subprocess.run([exe, *[str(x) for x in d]], capture_output=True, text=True, check=True).stdout.split('\n')
    if 'SOL' not in out: raise RuntimeError('z7search found no border for d = %s' % (list(d),))
    i = out.index('SOL')
    return [[int(x) for x in l.split()] for l in out[i + 1:i + 9]]
def make_border(a, d, rows=None):
    if rows is None: rows = z7search_rows(d)
    egen = [np.array(r) for r in rows[:4]]; fgen = [np.array(r) for r in rows[4:]]
    k = 8; b = 28
    I8 = np.eye(k, dtype=np.int64); J8 = np.ones((k, k), dtype=np.int64)
    Mc = np.kron(gs4(a), J8) + np.kron(gs4(d), I8)
    def shiftv(v, u):
        out = v.copy()
        for ii in range(4):
            for bb in range(7): out[ii * 8 + 1 + bb] = v[ii * 8 + 1 + ((bb + u) % 7)]
        return out
    E = np.array([shiftv(egen[j], aa) for j in range(4) for aa in range(7)])
    F = np.array([shiftv(fgen[j], aa) for j in range(4) for aa in range(7)]).T
    Kc = -(F.T @ Mc @ E.T)
    assert np.all(np.abs(Kc) == 32), 'compatibility fails'
    K = (Kc // 32).T
    assert np.array_equal(K @ K.T, b * np.eye(b, dtype=np.int64))
    assert np.array_equal(E @ E.T, 4 * k * np.eye(b, dtype=np.int64))
    assert np.array_equal(F @ F.T, 4 * k * np.eye(4 * k, dtype=np.int64) - 4 * np.kron(np.eye(4, dtype=np.int64), J8))
    assert np.array_equal(Mc @ E.T + F @ K.T, np.zeros((4 * k, b), dtype=np.int64))
    return K, E, F
if __name__ == '__main__':
    for a, d in (([2, 4, 4, 4], [4, -2, -2, 2]), ([2, 4, 4, 4], [-4, -2, 2, 2]), ([0, 4, 4, 4], [2, -2, 2, 4])):
        K, E, F = make_border(a, d); print('border OK for a', a, 'd', d, 'shapes', K.shape, E.shape, F.shape)
