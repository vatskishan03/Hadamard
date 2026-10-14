"""Vectorised exhaustive search for bordered T-matrices over small abelian groups G (|G| = t - 1)."""
import os
import sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import spence_lib as L
from bordered_od import od_matrix

def search(mods, limit=3):
    G = L.AbGroup(mods); v = G.v
    tr = np.stack([G.translate(g) for g in range(v)])
    sols = []
    digits = 8 ** np.arange(v)
    for start in range(0, 8 ** v, 1 << 18):
        idx = np.arange(start, min(start + (1 << 18), 8 ** v), dtype=np.int64)
        cfg = (idx[:, None] // digits) % 8
        which = cfg % 4; sign = np.where(cfg < 4, 1, -1)
        T = np.zeros((len(idx), 4, v), dtype=np.int64)
        np.put_along_axis(T, which[:, None, :], sign[:, None, :], axis=1)
        rs = np.sort(np.abs(T.sum(2)), axis=1)
        ok = (rs[:, 0] == 0) & (rs[:, 1] == 0) & (rs[:, 2] == 0) & (rs[:, 3] == 1)
        if not ok.any(): continue
        T = T[ok]
        P = np.zeros((len(T), v), dtype=np.int64)
        for g in range(1, v):
            P[:, g] = (T * T[:, :, tr[g]]).sum((1, 2))
        good = np.all(P[:, 1:] == -1, axis=1)
        for i in np.flatnonzero(good):
            sols.append(T[i])
            if len(sols) >= limit: return G, sols
    return G, sols

if __name__ == '__main__':
    rng = np.random.default_rng(0)
    for spec in sys.argv[1:]:
        mods = [int(m) for m in spec.split(',')]
        G, sols = search(mods)
        print(f"G=Z{mods} (t={G.v+1}): {len(sols)} found", flush=True)
        for T in sols[:1]:
            print('   T rows:', [''.join('+' if e > 0 else '-' if e < 0 else '0' for e in t) for t in T])
            ok = True
            for x in rng.standard_normal((5, 4)):
                H = od_matrix(G, T.astype(float), x)
                ok &= np.allclose(H @ H.T, (G.v + 1) * (x @ x) * np.eye(H.shape[0]))
            print(f"   OD(4t; t^4) identity for random real x: {ok}", flush=True)
