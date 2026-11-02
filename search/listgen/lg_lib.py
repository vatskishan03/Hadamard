"""Helpers for the list generator lg_pt (export of a block as a quadratic form in +-1 class signs, moments, tau)."""
import sys, struct, math, os, time
# cp_lowe, cyclo_product, workdir and verify are found through PYTHONPATH=search (see search/listgen/README.md)
import numpy as np
import cyclo_product as CP
import cp_lowe as CL

def make_structure(v, e, cfg):
    S = CP.Structure(v, v, e, skew=True)
    quads = CL.stage1_prime(S); y = quads[cfg]
    N = CL.cyclotomic_numbers(S)
    return S, y, N

def block_form(S, N, kind, f):
    """(m, const[D], lin[m][D], SI[m][m][D]) of the block: key_c(y) = const_c + sum lin y + sum_{i<j} SI y_i y_j.
    kind 'A' (skew, y = first e/2 signs, x = (y, -y)) or 'P' (e signs, value f at 0)."""
    e, h, nb = S.e, S.e // 2, S.nb
    N = N.astype(np.int64)
    if kind == 'A':
        m = h
        Mm = N[:, :h, :h] + N[:, h:, h:] - N[:, :h, h:] - N[:, h:, :h]
        lin = np.zeros((m, nb), dtype=np.int64)
    else:
        m = e
        Mm = N
        lin = np.zeros((m, nb), dtype=np.int64)
        for c in range(nb):
            lin[c, c] += f; lin[(c + h) % e, c] += f
    const = np.array([np.trace(Mm[c]) for c in range(nb)], dtype=np.int64)
    SI = np.zeros((m, m, nb), dtype=np.int64)
    for c in range(nb):
        Sc = Mm[c] + Mm[c].T
        np.fill_diagonal(Sc, 0)
        SI[:, :, c] = Sc
    assert np.abs(SI).max() < 30000
    return m, const, lin, SI

def keys_of(m, const, lin, SI, Y):
    """Keys of sign rows Y [n x m]."""
    n = len(Y); nb = len(const)
    K = np.tile(const, (n, 1)).astype(np.int64) + Y @ lin
    for c in range(nb):
        K[:, c] += ((Y @ SI[:, :, c]) * Y).sum(1) // 2
    return K

def sample_members(S, kind, s, n, rng):
    """Uniform member sign vectors (class signs x [n x e]) of a block: kind 'A' skew, 'P' with class sum s."""
    e, h = S.e, S.e // 2
    if kind == 'A':
        half = rng.integers(0, 2, size=(n, h)) * 2 - 1
        return np.concatenate([half, -half], axis=1)
    npos = (e + s) // 2
    X = -np.ones((n, e), dtype=np.int64)
    idx = np.argsort(rng.random((n, e)), axis=1)[:, :npos]
    np.put_along_axis(X, idx, 1, axis=1)
    return X

def moments(S, N, y, nsamp=400000, seed=1, chunk=100000):
    """Mean and covariance of the keys of uniform members of each of the four blocks (A, B, C, D)."""
    rng = np.random.default_rng(seed)
    out = []
    for k in range(4):
        if k == 0: kind, f, s = 'A', 1, 0
        else: f, s = CP.fs_of(S, y[k][0]); kind = 'P'
        m, const, lin, SI = block_form(S, N, kind, f)
        K = []
        for i in range(0, nsamp, chunk):
            X = sample_members(S, kind, s, min(chunk, nsamp - i), rng)
            Yv = X[:, :m] if kind == 'A' else X
            K.append(keys_of(m, const, lin, SI, Yv))
        K = np.concatenate(K).astype(float)
        out.append(dict(kind=kind, f=f, s=s, mu=K.mean(0), Sig=np.cov(K.T)))
    return out

def energy_setup(mom, k):
    """ck, W (Z x D), Z for block k: conditional mean given that the four keys sum to 0."""
    mus = [b['mu'] for b in mom]; sigs = [b['Sig'] for b in mom]
    Sp = np.linalg.pinv(sum(sigs), rcond=1e-10); mutot = sum(mus)
    ck = mus[k] + sigs[k] @ Sp @ (0 - mutot)
    W, Z = CL.whitening(sigs[k])
    return ck, W, Z

def export_block(path, S, N, mom, k, e=None, y=None, ck=None, W=None):
    kind = mom[k]['kind']; f = mom[k]['f']; s = mom[k]['s']
    m, const, lin, SI = block_form(S, N, kind, f)
    if ck is None: ck, W, Z = energy_setup(mom, k)
    Z, D = W.shape
    e = S.e
    mode = 1 if kind == 'A' else 0
    npos = (e + s) // 2 if kind == 'P' else 0
    with open(path, 'wb') as fo:
        fo.write(struct.pack('12i', m, D, Z, mode, npos, e, 1 if kind == 'A' else 0, 0, 0, 0, 0, 0))
        fo.write(np.ascontiguousarray(W, dtype='<f8').tobytes()); fo.write(np.asarray(ck, dtype='<f8').tobytes())
        fo.write(const.astype('<i4').tobytes()); fo.write(lin.astype('<i4').tobytes())
        fo.write(np.ascontiguousarray(SI.transpose(0, 1, 2), dtype='<i2').tobytes())
    return dict(m=m, D=D, Z=Z, mode=mode, npos=npos, kind=kind)

def class_mask_to_x(cm, e):
    """Sign vectors [n x e] from class masks (bit c set: sign -1)."""
    cm = np.asarray(cm, dtype=np.uint64)
    return 1 - 2 * ((cm[:, None] >> np.arange(e, dtype=np.uint64)[None, :]) & np.uint64(1)).astype(np.int64)

def canon(cm, e):
    h = e // 2
    full = np.uint64((1 << e) - 1) if e < 64 else np.uint64(0xFFFFFFFFFFFFFFFF)
    cm = np.asarray(cm, dtype=np.uint64)
    r = ((cm >> np.uint64(h)) | (cm << np.uint64(h))) & full
    return np.minimum(cm, r)

def read_records(path, D):
    dt = np.dtype([('mask', '<u8'), ('key', '<i2', (D,))])
    return np.fromfile(path, dtype=dt)
