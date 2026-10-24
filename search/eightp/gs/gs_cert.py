"""Certificates (verify/verify.py) for K-level solutions of the real GS-type plug-ins.
 mode gs     -> construction "goethals-seidel-group"   (order 4|G|)
 mode spence -> construction "spence"                  (order 4|G| + 4; row sums normalised to (2,0,0,0) by negations/reordering; P, Q as in hadamard/2292)
 mode bord   -> construction "bordered-goethals-seidel" (order 4(|G| + t); border solved with search/bordered.py (Q = Z_3, T = all characters: border3.c); other (Q, T): unsupported -> raises)
ODD |G| (search/NOTES.md section 92): row sums (+-1, +-1, +-1, +-1), P = Q = S, K = -S built by verify.py after oddspence/verify_odd_spence.patch (the repository's verify.py accepts only (2, 0, 0, 0)).
The certificate is written into a SCRATCH copy of the repository (`repo_copy` in this directory, see README.md; it needs verify/ and search/), built and checked there with python3 verify/verify.py <order>; the repository itself is never modified.
"""
import os, sys, json, hashlib, subprocess, importlib.util, itertools
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from gsk import *
REPO = os.path.join(HERE, 'repo_copy')
os.environ.setdefault('HADAMARD_WORK', os.path.join(HERE, 'work'))
SPENCE_P = [[-1, 1, 1, 1], [1, -1, 1, 1], [1, 1, -1, 1], [1, 1, 1, -1]]
SPENCE_Q = [[1, 1, 1, 1], [1, -1, 1, -1], [1, 1, -1, -1], [1, -1, -1, 1]]


def _verify_module(repo):
    sys.path.insert(0, os.path.join(repo, 'verify'))
    spec = importlib.util.spec_from_file_location('verify_' + str(abs(hash(repo))), os.path.join(repo, 'verify', 'verify.py'))
    V = importlib.util.module_from_spec(spec); spec.loader.exec_module(V)
    return V


def seq_str(x): return ''.join('+' if v == 1 else '-' for v in x)


def automorphism_matrices(sh):
    """integer matrices on the coordinates of G = K x F_r (prime r only) for the generators of H: (u_1..u_k, c) with c = square (phi = +1) or non-square (phi = -1), plus {1} x S"""
    from idx import GF, factor_prime_power
    p, e = factor_prime_power(sh.r)
    if e != 1: return None
    g = next(a for a in range(2, sh.r) if all(pow(a, (sh.r - 1) // f, sh.r) != 1 for f in range(2, sh.r) if (sh.r - 1) % f == 0 and all(f % d for d in range(2, f)))) if sh.r > 3 else 2
    k = len(sh.Kmods); mats = []
    sq = (g * g) % sh.r; nsq = g % sh.r
    M = np.eye(k + 1, dtype=int); M[k, k] = sq; mats.append(M.tolist())
    for (u, ph) in sh.hgens_raw:
        uu = [int(x) for x in (u if not np.isscalar(u) else [u] * k)]
        M = np.eye(k + 1, dtype=int)
        for i in range(k): M[i, i] = uu[i]
        M[k, k] = sq if ph == 1 else nsq
        mats.append(M.tolist())
    return mats


def build_cert(sh, blocks):
    ok, seqs = verify_g(sh, blocks)
    assert ok, 'G-level verification failed'
    mods = sh.gmods(); v = sh.G
    seqs = [np.asarray(x, dtype=np.int64) for x in seqs]
    if sh.mode == 'gs':
        cert = dict(order=4 * v, construction='goethals-seidel-group', group=mods, sequences=[seq_str(x) for x in seqs], rowsums=[int(x.sum()) for x in seqs])
        am = automorphism_matrices(sh)
        if am: cert['automorphisms'] = am
        return cert
    if sh.mode == 'spence' and all(abs(int(x.sum())) == 1 for x in seqs):
        # ODD |G|: row sums (+-1)^4, no normalisation needed; the patched verify.py builds P = Q = S (search/NOTES.md section 92)
        return dict(order=4 * v + 4, construction='spence', group=mods, P=[[0] * 4] * 4, Q=[[0] * 4] * 4, sequences=[seq_str(x) for x in seqs])
    if sh.mode == 'spence':
        sums = [int(x.sum()) for x in seqs]
        idx = [i for i, s in enumerate(sums) if abs(s) == 2]
        assert len(idx) == 1 and sorted(abs(s) for s in sums) == [0, 0, 0, 2], f'row sums {sums} are not a permutation of (+-2,0,0,0)'
        i0 = idx[0]; first = seqs[i0] * (1 if sums[i0] > 0 else -1)
        rest = [seqs[i] for i in range(4) if i != i0]
        seqs = [first] + rest
        cert = dict(order=4 * v + 4, construction='spence', group=mods, P=SPENCE_P, Q=SPENCE_Q, sequences=[seq_str(x) for x in seqs])
        return cert
    # bordered
    t = sh.t_border; Kg = sh.Kg; m = np.array(sh.Kmods)
    # L = common kernel of the characters
    Lel = [k for k in range(sh.N) if all(abs(((Kg.el[k] * np.array(c)) / m).sum() - round(((Kg.el[k] * np.array(c)) / m).sum())) < 1e-9 for c in sh.chars)]
    qidx = sh.N // len(Lel)
    gens = []
    nk = len(sh.Kmods)
    p_, e_ = __import__('idx').factor_prime_power(sh.r)
    for k in Lel:
        if k: gens.append([int(x) for x in Kg.el[k]] + [0] * e_)
    for j in range(e_): gens.append([0] * nk + [1 if i == j else 0 for i in range(e_)])
    # target values off K' must be constant
    nonL = [k for k in range(sh.N) if k not in set(Lel)]
    vals = {int(sh.phi[k]) for k in nonL}
    assert len(vals) <= 1, 'Gamma is not constant off L: the verify builder supports only a two-valued target'
    sub = int(sh.phi[0]); els = vals.pop() if vals else 0
    cert = dict(order=4 * (v + t), construction='bordered-goethals-seidel', group=mods, subgroup_generators=gens)
    am = automorphism_matrices(sh)
    if am: cert['symmetries'] = am
    cert['autocorrelation'] = dict(subgroup=sub, elsewhere=els)
    cert['rowsums'] = [int(x.sum()) for x in seqs]
    cert['sequences'] = [seq_str(x) for x in seqs]
    # border (Q = Z_3, T = all characters, border of size 12)
    assert qidx == 3 and t == 3 and sub == -12 and els == 0, 'border solver implemented for Q = Z_3 with T = all characters only'
    coords = np.array(np.unravel_index(np.arange(v), mods)).T
    # coset numbering of the verify builder: cosets of K' in order of smallest element
    Kp = set()
    for k in Lel:
        for a in range(sh.r): Kp.add(k * sh.r + a)
    cos = -np.ones(v, dtype=np.int64); q = 0
    for g in range(v):
        if cos[g] < 0:
            # coset of g = g + K'
            for k in Lel:
                kk = int(Kg.add[g // sh.r, k])
                for a in range(sh.r): cos[kk * sh.r + a] = q
            q += 1
    yc = [[int(x[cos == a].sum()) for a in range(3)] for x in seqs]
    sys.path.insert(0, os.path.join(REPO, 'search'))
    import bordered as BD
    if all(abs(v) == 1 for y in yc for v in y):
        # |K'| odd: all compressions are +-1, Mc (the Goethals-Seidel array of the compressions over Z3 in the layout of verify.build_bordered_gs) is a Hadamard matrix of order 12
        # and X = Paley(12), Chat = X, Rhat = -Mc solves the border problem:  X Chat^T + Rhat Mc^T = 12 I - 12 I = 0,  Chat Chat^T = Rhat Rhat^T = 12 I
        sys.path.insert(0, os.path.join(REPO, 'search'))
        import bordered as BD
        yy = [np.array(y, dtype=np.int64) for y in yc]
        idx = np.arange(3); one = (idx[None, :] - idx[:, None]) % 3; two = (idx[:, None] + idx[None, :]) % 3; ng = (-idx) % 3
        A = yy[0][one]; B, C, D = (y[two] for y in yy[1:]); Br, Cr, Dr = (y[ng][two] for y in yy[1:])
        Mc = np.block([[A, B, C, D], [-B, A, Dr, -Cr], [-C, -Dr, A, Br], [-D, Cr, -Br, A]])
        assert np.array_equal(Mc @ Mc.T, 12 * np.eye(12, dtype=np.int64)), 'compression array is not Hadamard'
        X = BD.paley12(); Rh = -Mc; Ch = X
        cert['border'] = dict(X=[seq_str(r) for r in X], Rhat=[seq_str(r) for r in Rh], Chat=[seq_str(r) for r in Ch])
        return cert
    sys.path.insert(0, os.path.join(REPO, 'search'))
    import bordered as BD
    # bordered.gs_array uses A[g,h] = a(g-h) (index() transposes), the verify builder A[g,h] = a(h-g): reflect the first compression
    yc = [[yc[0][(-a) % 3] for a in range(3)]] + yc[1:]
    b = BD.border(yc)
    if b is None: raise RuntimeError(f'no border for compressions {yc}')
    X, Rh, Ch = b
    cert['border'] = dict(X=[seq_str(r) for r in X], Rhat=[seq_str(r) for r in Rh], Chat=[seq_str(r) for r in Ch])
    return cert


def certify_solution(sh, blocks, tag='', repo=REPO):
    cert = build_cert(sh, blocks)
    V = _verify_module(repo)
    H, _ = V.BUILDERS[cert['construction']](cert)
    order = H.shape[0]
    assert order == cert['order'], (order, cert['order'])
    assert V.is_pm1(H) and V.gram_is_scalar(H, order), 'H H^T != nI'
    h = hashlib.sha256()
    for ch in V.csv_chunks(H): h.update(ch)
    cert['sha256'] = h.hexdigest()
    name = str(order)
    d = os.path.join(repo, 'hadamard', name)
    if os.path.exists(d) and tag != '__force__': d = os.path.join(repo, 'hadamard', f'T{order}_{tag}')
    os.makedirs(d, exist_ok=True)
    json.dump(cert, open(os.path.join(d, 'certificate.json'), 'w'))
    os.makedirs(os.path.join(HERE, 'out'), exist_ok=True)
    json.dump(cert, open(os.path.join(HERE, 'out', f'{sh.mode}_{order}_{tag}.json'), 'w'))
    rel = os.path.relpath(d, repo)
    arg = os.path.basename(d) if os.path.basename(d).isdigit() else rel
    r = subprocess.run([sys.executable, 'verify/verify.py', arg], capture_output=True, text=True, cwd=repo)
    return dict(order=order, dir=d, sha256=cert['sha256'], verify_ok=(r.returncode == 0 and 'ok' in r.stdout.lower()), verify_out=(r.stdout.strip() + r.stderr.strip())[:200])
