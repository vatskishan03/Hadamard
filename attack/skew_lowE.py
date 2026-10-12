"""Skew Goethals-Seidel search for one structure and row-sum pattern with lists of unequal sizes.

skew_ktree.py gives every role the same list size L. With the skew block A enumerated completely (N_A
sequences) the k-tree work is N_A N_D + N_B N_C and the yield is proportional to N_B N_C N_D, so for a
fixed memory budget it is better to take N_B N_C close to N_A N_D: a long list for D and shorter ones for
B and C. B, C and D are low-energy candidates from skew_gen.c (local minima of the distance of the
autocorrelations from their mean), which narrows each class distribution about twofold.

  python3 attack/skew_lowE.py v g[,g2] b,c,d --nbc N --nd M [--steps S] [--maxtau T] [--top K] [--seed R] [--amax N]

Roles with equal row sums share one list. Before matching, the Gaussian estimate of the number of solutions
in the four lists is printed (computed from the lists' own class means and variances; v prime). Hits are
verified (H H^T = 4vI and H + H^T = 2I) and appended to $HADAMARD_WORK/skew_hits.jsonl.
"""
import argparse, json, math, os, subprocess, sys, tempfile, time
import numpy as np
from pathlib import Path as _P; sys.path[:0] = [str(_P(__file__).resolve().parent), str(_P(__file__).resolve().parent.parent / 'verify')]
import workdir as W
import spence_lib as L
import skew_ktree as K
import verify


def log10_yield(Vs, pm, mult, csize):
    """Gaussian estimate of log10 #{quadruples with sum_i v_i = -pm on every class}; mult[i] = list size.
    The row sums fix sum_c csize_c (sum_i P_i(c)), so one class is implied and left out (as in skew_rank.py)."""
    lg = sum(math.log10(m) for m in mult)
    terms = []
    for c in range(Vs[0].shape[1]):
        mu = sum(V[:, c].astype(np.float64).mean() for V in Vs) + pm
        var = sum(V[:, c].astype(np.float64).var() for V in Vs)
        terms.append((csize[c] * math.sqrt(var), (-mu * mu / (2 * var)) / math.log(10) - 0.5 * math.log10(2 * math.pi * var)))
    terms.sort()
    return lg + sum(t for _, t in terms[:-1])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('v', type=int); ap.add_argument('gens'); ap.add_argument('pattern')
    ap.add_argument('--nbc', type=float, default=2e6); ap.add_argument('--nd', type=float, default=2e7)
    ap.add_argument('--steps', type=int, default=8); ap.add_argument('--maxtau', type=int, default=1000)
    ap.add_argument('--top', type=int, default=8); ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--cap', type=float, default=2.5e7); ap.add_argument('--sols', type=int, default=4)
    ap.add_argument('--amax', type=float, default=2 ** 22, help='enumerate all skew A up to this many, else sample this many')
    a = ap.parse_args()
    v = a.v; gens = [int(g) for g in a.gens.split(',')]; b, c, d = map(int, a.pattern.split(','))
    G = L.AbGroup([v]); st = L.Structure(G, [G.perm_from_matrix(np.array([[g]])) for g in gens])
    rng = np.random.default_rng(a.seed); pm = v % 4
    log = lambda s: print(s, flush=True)
    log(f"v={v} H={gens} orbits={st.O} classes={st.C} rows=(1,{b},{c},{d}) seed={a.seed} nbc={a.nbc:.3g} nd={a.nd:.3g} steps={a.steps}")
    rep = np.array([np.flatnonzero(st.orb == o)[0] for o in range(st.O)])
    partner = st.orb[G.neg()[rep]]
    pairs = [(o, partner[o]) for o in range(1, st.O) if o < partner[o]]
    if 2 ** len(pairs) <= a.amax:                                  # every skew A
        s = 1 - 2 * ((np.arange(2 ** len(pairs))[:, None] >> np.arange(len(pairs))) & 1)
    else:                                                          # a uniform sample of distinct skew A
        s = np.unique(rng.choice(np.array([-1, 1], dtype=np.int8), size=(int(a.amax * 1.02), len(pairs))), axis=0)[:int(a.amax)]
    YA = np.zeros((len(s), st.O), dtype=np.int8); YA[:, 0] = 1
    for j, (o, q) in enumerate(pairs): YA[:, o] = s[:, j]; YA[:, q] = -s[:, j]
    VA, ok = K.values(st, YA, pm); YA, VA = YA[ok], VA[ok]
    t0 = time.time()
    # one list per distinct row sum, as large as the largest role that uses it; roles with equal row sums share it
    want = {}
    for r, n in ((b, a.nbc), (c, a.nbc), (d, a.nd)): want[r] = max(want.get(r, 0), int(n))
    lists = {}
    for r, n in want.items():
        Y = K.gen_pool(st, r, n, a.steps, rng, v); V, ok = K.values(st, Y, pm); lists[r] = (Y[ok], V[ok])
    (YB, VB), (YC, VC), (YD, VD) = lists[b], lists[c], lists[d]
    log(f"lists A={len(YA)} B={len(YB)} C={len(YC)} D={len(YD)} [{time.time()-t0:.0f}s]; full pools "
        f"{st.pool(b):.3g} {st.pool(c):.3g} {st.pool(d):.3g}")
    log(f"Gaussian estimate of solutions inside the lists: 10^{log10_yield([VA, VB, VC, VD], pm, [len(YA), len(YB), len(YC), len(YD)], st.csize):.2f}")
    roles = [(YA, VA), (YB, VB), (YC, VC), (YD, VD)]          # match4 pairs (A, D) and (B, C)
    S = K.clamp((VA, VD), (VB, VC), a.cap); S = max(1, min(S, st.C - 2))
    tmp = tempfile.mkdtemp(dir=W.tmpdir())
    paths = [os.path.join(tmp, f"{k}.bin") for k in range(4)]
    K.write(paths[0], VA)
    byrow = {}
    for k, r in ((1, b), (2, c), (3, d)):
        if r in byrow: paths[k] = byrow[r]
        else: K.write(paths[k], roles[k][1]); byrow[r] = paths[k]
    cmd = [W.bin('match4'), '-s', str(S), '-t', str(-pm), '-n', str(a.maxtau), '-k', str(a.sols), '-m', str(a.top)] + paths + [str(st.C)]
    log(f"clamp S={S}; {' '.join(cmd[:11])} ...")
    t0 = time.time()
    errlog = W.WORK / 'logs' / f'match4_{v}_{a.seed}.log'          # progress of the matcher, one line per tau
    errlog.parent.mkdir(parents=True, exist_ok=True)
    with open(errlog, 'w') as ef:
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=ef, text=True)
        out, _ = proc.communicate()
    err = errlog.read_text()
    for p in set(paths): os.remove(p)
    os.rmdir(tmp)
    for line in err.strip().splitlines()[-4:]: log("   " + line)
    found = []
    for line in out.splitlines():
        if not line.startswith('SOL'): continue
        ids = [int(x, 16) for x in line.split()[1:]]
        seqs = [roles[k][0][ids[k]].astype(np.int64)[st.orb] for k in range(4)]
        tot = verify.paf_sum(seqs)
        if not (tot[0] == 4 * v and not tot[1:].any()): continue
        H = verify.goethals_seidel(*map(verify.circulant, seqs)).astype(np.float64)
        had = np.array_equal(H @ H.T, 4 * v * np.eye(4 * v)); skew = np.array_equal(H + H.T, 2 * np.eye(4 * v))
        rs = [int(x.sum()) for x in seqs]
        log(f"VERIFIED rows={rs}: Hadamard={had} skew={skew} [{time.time()-t0:.0f}s]")
        found.append(dict(n=v, gens=gens, pattern=rs, hadamard=bool(had), skew=bool(skew), seed=a.seed,
                          seqs=[''.join('+' if x > 0 else '-' for x in sq) for sq in seqs]))
    if found:
        with open(W.hits(), 'a') as f:
            for h in found: f.write(json.dumps(h) + '\n')
    log(f"done: {len(found)} verified [{time.time()-t0:.0f}s matching]")


if __name__ == '__main__':
    main()
