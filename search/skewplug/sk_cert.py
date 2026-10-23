"""certificates for skew plug-in solutions: K-level (s,t,e) blocks -> +-1 sequences over G = K x F_q (mixed radix, index k*q + a) -> 'goethals-seidel-group' certificate, "skew": true,
built and checked with the repository's verify/verify.py builders (H H^T = nI, H + H^T = 2I, SHA-256) and written to out/skew_<order>_<tag>.json."""
import os, sys, json, hashlib, importlib.util
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from sk import SShape
REPO = os.path.abspath(os.path.join(HERE, '..', '..'))


def seq_str(x): return ''.join('+' if v == 1 else '-' for v in x)


def primitive_root(q):
    for g in range(2, q):
        if len({pow(g, i, q) for i in range(1, q)}) == q - 1: return g
    return 1


def automorphisms(sh):
    q = sh.q; k = len(sh.Kmods); mats = []
    if q > 3:
        M = np.eye(k + 1, dtype=int); M[k, k] = pow(primitive_root(q), 2, q); mats.append(M.tolist())
    for u in sh.hunits:
        M = np.eye(k + 1, dtype=int)
        if len(u) and not np.isscalar(u[0]):
            for i in range(k):
                for j in range(k): M[i, j] = u[i][j]
        else:
            for i in range(k): M[i, i] = u[i]
        mats.append(M.tolist())
    return mats


def build_cert(sh, blocks):
    seqs = sh.glevel_blocks(blocks)
    assert sh.g_target_skew_ok(seqs), 'G-level skew check failed'
    seqs = [np.asarray(x, dtype=np.int64) for x in seqs]
    cert = dict(order=4 * sh.G, construction='goethals-seidel-group', skew=True, group=sh.gmods(), sequences=[seq_str(x) for x in seqs],
                rowsums=[int(x.sum()) for x in seqs])
    cert['automorphisms'] = automorphisms(sh)
    return cert


def certify(sh, blocks, tag='', repo=REPO):
    cert = build_cert(sh, blocks)
    sys.path.insert(0, os.path.join(repo, 'verify'))
    spec = importlib.util.spec_from_file_location('verify_' + str(abs(hash(repo))), os.path.join(repo, 'verify', 'verify.py'))
    V = importlib.util.module_from_spec(spec); spec.loader.exec_module(V)
    H, _ = V.BUILDERS[cert['construction']](cert)
    order = H.shape[0]
    assert order == cert['order']
    assert V.is_pm1(H) and V.gram_is_scalar(H, order), 'H H^T != nI'
    assert np.array_equal(H + H.T, 2 * np.eye(order, dtype=H.dtype)), 'not skew'
    h = hashlib.sha256()
    for ch in V.csv_chunks(H): h.update(ch)
    cert['sha256'] = h.hexdigest()
    os.makedirs(os.path.join(HERE, 'out'), exist_ok=True)
    path = os.path.join(HERE, 'out', f'skew_{order}_{tag}.json')
    json.dump(cert, open(path, 'w'))
    return dict(order=order, file=path, sha256=cert['sha256'], verify_ok=True)
