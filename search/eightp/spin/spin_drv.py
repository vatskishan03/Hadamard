"""spin_drv.py -- d = 3 block-symmetric Spence structures over Z_v (v odd): blocks (A, B, mu B, mu^2 B), A invariant under K = <H, mu>, B invariant under H, [K:H] = 3.
Condition  PAF_A(g) + PAF_B(g) + PAF_B(mu g) + PAF_B(mu^2 g) = -4  for all g != 0, row sums rho_A = rho_B = +1 (after negation).  Exact: Gray-code enumeration of B (spin3.c), hash lookup of A.
usage: spin_drv.py v OBmax [OBmin]   (search/NOTES.md section 92; the C engine spin3 is built in this directory)
NOTE: the COUNT printed is the number of pairs (B, distinct PAF vector of an admissible A): admissible A with the same PAF vector are collapsed in the hash table, so it is not the number of ordered (A, B);
zero / nonzero is exact."""
import sys, os, math, subprocess, itertools, time, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(os.path.dirname(HERE), 'plain'))
import plain_spence3 as ps

def pair_forms(v, lab, O, g):
    idx = np.arange(v); sh = (idx + g) % v
    M = np.zeros((O, O), dtype=np.int64); np.add.at(M, (lab, lab[sh]), 1)
    return M

def setup(v, H, mu):
    K = ps.closure(list(H) + [mu], v)
    labH, repsH = ps.orbits(v, H); labK, repsK = ps.orbits(v, K)
    OH = len(repsH); OK = len(repsK)
    sizeH = np.bincount(labH, minlength=OH); sizeK = np.bincount(labK, minlength=OK)
    # classes of nonzero K-orbits under negation
    negl = labK[(-np.array(repsK)) % v]
    seen = set(); cls = []
    for a in range(1, OK):
        if a in seen: continue
        seen.add(a); seen.add(int(negl[a])); cls.append(a)
    gs = [repsK[a] for a in cls]
    return K, labH, OH, sizeH, labK, OK, sizeK, gs

def build_input(v, H, mu, path):
    K, labH, OH, sizeH, labK, OK, sizeK, gs = setup(v, H, mu)
    C = len(gs)
    SB = np.zeros((C, OH, OH), dtype=np.int64)
    for c, g in enumerate(gs):
        for j in range(3):
            SB[c] += pair_forms(v, labH, OH, (pow(mu, j, v) * g) % v)
    # A list: all x in {+-1}^OK with row sum +1
    XA = ps.enum_list(OK, sizeK, rowsum=1)
    MA = [pair_forms(v, labK, OK, g) for g in gs]
    XAi = XA.astype(np.int64)
    VA = np.stack([np.einsum('ni,ij,nj->n', XAi, M, XAi) for M in MA], 1)
    VAu, first = np.unique(VA, axis=0, return_index=True)       # distinct vectors (first representative)
    T = -4 - VAu
    with open(path, 'w') as f:
        f.write(f"{C} {OH} {len(VAu)}\n")
        f.write(' '.join(map(str, sizeH.tolist())) + '\n')
        f.write(' '.join(map(str, SB.ravel().tolist())) + '\n')
        for row in T: f.write(' '.join(map(str, row.tolist())) + '\n')
    return dict(C=C, OH=OH, OK=OK, NA=len(VAu), XA=XA, first=first, labH=labH, labK=labK, sizeH=sizeH, sizeK=sizeK, K=K)

def verify_solution(v, info, gray, aidx, mu):
    OH = info['OH']; bits = [(gray >> i) & 1 for i in range(OH)]
    xB = np.array([1 - 2 * b for b in bits], dtype=np.int64)
    xA = info['XA'][info['first'][aidx]].astype(np.int64)
    A = xA[info['labK']]; B = xB[info['labH']]
    h = np.arange(v)
    seqs = [A, B, B[(mu * h) % v], B[(mu * mu * h) % v]]
    tot = np.zeros(v, dtype=np.int64)
    for s in seqs:
        for g in range(v): tot[g] += int(s @ s[(h + g) % v])
    ok = tot[0] == 4 * v and np.all(tot[1:] == -4) and all(int(s.sum()) == 1 for s in seqs)
    return bool(ok), seqs

def rho_feasible(size):
    tot = int(sum(size)); off = tot
    dp = np.zeros(2 * tot + 1); dp[off] = 1
    for s in size:
        s = int(s); new = np.zeros_like(dp); new[s:] += dp[:-s]; new[:-s] += dp[s:]; dp = new
    return dp[off + 1] > 0

def run_pair(v, H, mu, tmp=None):
    tmp = tmp or f'tmp_in_{os.getpid()}.txt'
    path = os.path.join(HERE, os.path.basename(tmp))
    K, labH, OH, sizeH, labK, OK, sizeK, gs = setup(v, H, mu)
    if not rho_feasible(sizeK) or not rho_feasible(sizeH):
        return 0, [], dict(OH=OH, OK=OK, C=len(gs), NA=0)
    info = build_input(v, H, mu, path)
    if info['NA'] == 0: return 0, [], info
    r = subprocess.run([os.path.join(HERE, 'spin3'), path], capture_output=True, text=True)
    count = 0; sols = []
    for line in r.stdout.splitlines():
        if line.startswith('COUNT'): count = int(line.split()[1])
        elif line.startswith('SOL'):
            _, g, a = line.split(); sols.append((int(g), int(a)))
    return count, sols, info

def brute_pair(v, H, mu):
    K, labH, OH, sizeH, labK, OK, sizeK, gs = setup(v, H, mu)
    XB = ps.enum_list(OH, sizeH, rowsum=1); XA = ps.enum_list(OK, sizeK, rowsum=1)
    h = np.arange(v); cnt = 0
    for xa in XA:
        A = xa.astype(np.int64)[labK]
        pa = np.array([int(A @ A[(h + g) % v]) for g in range(v)])
        for xb in XB:
            B = xb.astype(np.int64)[labH]
            tot = pa.copy()
            for j in range(3):
                Bj = B[(pow(mu, j, v) * h) % v]
                tot += np.array([int(Bj @ Bj[(h + g) % v]) for g in range(v)])
            if tot[0] == 4 * v and np.all(tot[1:] == -4): cnt += 1
    return cnt

def small_gens(H, v):
    gens = []; cur = frozenset([1])
    for e in sorted(H):
        if e in cur: continue
        gens.append(e); cur = ps.closure(gens, v)
        if len(cur) == len(H): break
    return gens

def pairs_for(v, OBmax, OBmin=2):
    subs = ps.all_subgroups(v)
    U = ps.units(v)
    out = {}
    for H in subs:
        # quick orbit count without the python loops: number of orbits = sum over divisors d of v of phi(d) / |image of H in Z_d^*|
        labH, repsH = ps.orbits(v, H); OH = len(repsH)
        if not (OBmin <= OH <= OBmax): continue
        Hs = set(H); gH = small_gens(H, v)
        cosets = set()
        for mu in U:
            if mu in Hs: continue
            if pow(mu, 3, v) not in Hs: continue
            K = ps.closure(gH + [mu], v)
            key = (H, K)
            if key not in out: out[key] = mu
    return [(H, mu) for (H, K), mu in out.items()]

if __name__ == '__main__':
    v = int(sys.argv[1]); OBmax = int(sys.argv[2]); OBmin = int(sys.argv[3]) if len(sys.argv) > 3 else 2
    prs = pairs_for(v, OBmax, OBmin)
    print(f'v = {v}: {len(prs)} pairs (H, mu) with {OBmin} <= O_H <= {OBmax}', flush=True)
    res = []
    for H, mu in sorted(prs, key=lambda hm: (len(ps.orbits(v, hm[0])[1]), -len(hm[0]))):
        lab, reps = ps.orbits(v, H); OH = len(reps)
        t0 = time.time()
        try:
            count, sols, info = run_pair(v, H, mu)
        except Exception as ex:
            print('ERR', OH, len(H), mu, repr(ex)[:100], flush=True); continue
        rec = dict(v=v, OH=OH, h=len(H), mu=mu, OK=info['OK'], C=info['C'], NA=info['NA'], count=count, secs=round(time.time() - t0, 1))
        if count:
            okv, seqs = verify_solution(v, info, sols[0][0], sols[0][1], mu)
            rec['verified'] = okv; rec['sol'] = [s.tolist() for s in seqs] if okv else None
            print('FOUND', rec, flush=True)
        print({k: val for k, val in rec.items() if k != 'sol'}, flush=True)
        res.append(rec)
