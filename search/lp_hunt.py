"""Hunt for Legendre pairs of length l = 2n - 1 for the unresolved orders 4n of FRONTIER.md (a Legendre pair of length l
gives a Hadamard matrix of order 2l + 2 = 4n), over the H-invariant structures of lp_lowe.py.

For every factorisation l = m q (q prime, q not dividing m), every even e | q - 1 with F = (q - 1)/e >= 2, every cyclic
K_m = <k> <= Z_m^* and every twist t (|K_m| t = 0 mod e) whose blocks are small enough (at most 16 position orbits and
C(e, e/2)^orbits <= 10^14), stage 1 is listed exactly (lp_lowe.lp_stage1). For each stage-1 class (up to --maxcls):
  - blocks of at most --exact members each are enumerated in full and joined exactly (a complete search of the class);
  - larger classes get the Gaussian estimate of the number of pairs and of the share kept by low-energy lists of --L
    members; nothing is searched (lp_lowe.py does that).
One JSON line per class in $HADAMARD_WORK/lp_hunt.jsonl; verified pairs also go to lp_hits.jsonl.

  python3 search/lp_hunt.py [n1,n2,..] [--exact 2e7] [--L 3e7] [--maxcls 40] [--samples 50000] [--maxopts 2e6]
"""
import argparse, json, math, os, re, sys, time
import numpy as np
from math import comb, gcd
from pathlib import Path as _P
HERE = _P(__file__).resolve().parent
sys.path[:0] = [str(HERE)]
import workdir as W
import cyclo_product as CP
import cp_lowe as CL
import cp_sweep as SW
import lp_lowe as LP


def frontier():
    ns = []
    for l in open(HERE.parent / 'FRONTIER.md'):
        m = re.match(r'\|\s*(\d+)\s*\|\s*(\d+)\s*\|', l)
        if m: ns.append(int(m.group(2)))
    return ns


def factor(n):
    out, d = {}, 2
    while d * d <= n:
        while n % d == 0: out[d] = out.get(d, 0) + 1; n //= d
        d += 1
    if n > 1: out[n] = out.get(n, 0) + 1
    return out


def korbits(m, k):
    seen, orbs = set(), []
    for a in range(m):
        if a in seen: continue
        o, x = [], a
        while x not in o: o.append(x); x = x * k % m
        seen |= set(o); orbs.append(o)
    return orbs


def lp_structures(l, maxorb=16, maxlog=14):
    for q, mult in sorted(factor(l).items()):
        if mult > 1: continue
        m = l // q
        if m < 2 or m > 4000: continue
        subs = None
        for e in range(2, min(q, 41), 2):
            if (q - 1) % e or (q - 1) // e < 2: continue
            if subs is None: subs = SW.cyclic_subgroups(m)
            for k in subs:
                orbs = korbits(m, k)
                if len(orbs) > maxorb or len(orbs) * math.log10(comb(e, e // 2)) > maxlog: continue
                Kord = next(len(o) for o in orbs if gcd(o[0], m) == 1)
                for t in range(e):
                    if Kord * t % e or (t and Kord == 1): continue
                    yield q, m, e, k, t, len(orbs)


def search_exact(S, O, P, rng):
    """Enumerate both blocks in full and join exactly; returns the verified pairs (as +- strings)."""
    D = P['D']
    Rv = np.random.default_rng(11).integers(1, 2 ** 62, size=D)
    with np.errstate(over='ignore'):
        ht = int((P['t'].astype(np.int64) * Rv).sum())
    res = []
    for k in range(2):
        blk = P['blocks'][k][0]
        Wk = np.eye(D)[:1]; ck = np.zeros(D)                 # tau = inf: every member is kept
        masks, hs, bits, msg = LP.enum_hashes(blk, Wk, ck, 3.0e38, Rv, None)
        res.append((masks, hs, bits))
    (mA, hA, bA), (mB, hB, bB) = res
    o = np.argsort(hB, kind='stable'); hBs = hB[o]
    with np.errstate(over='ignore'):
        need = np.int64(ht) - hA
    lo = np.searchsorted(hBs, need, 'left'); hi = np.searchsorted(hBs, need, 'right')
    out = []
    for i in np.flatnonzero(hi > lo):
        for j in o[lo[i]:hi[i]]:
            A = S.sequence(CL.decode(P['blocks'][0][0], mA[i], bA)); B = S.sequence(CL.decode(P['blocks'][1][0], mB[j], bB))
            if LP.is_legendre_pair(A, B) and A.sum() == 1 and B.sum() == 1:
                out.append((''.join('+' if x > 0 else '-' for x in A), ''.join('+' if x > 0 else '-' for x in B)))
                if len(out) >= 3: return out
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('ns', nargs='?', default='')
    ap.add_argument('--exact', type=float, default=2e7); ap.add_argument('--L', type=float, default=3e7)
    ap.add_argument('--maxcls', type=int, default=40); ap.add_argument('--samples', type=int, default=50000)
    ap.add_argument('--maxopts', type=float, default=2e6)
    a = ap.parse_args()
    ns = [int(x) for x in a.ns.split(',')] if a.ns else frontier()
    out = W.WORK / 'lp_hunt.jsonl'
    done = set()
    if out.exists():
        for l in open(out):
            r = json.loads(l); done.add((r['n'], r['q'], r['e'], r['k'], r['t'], r.get('cls', -1)))
    for n in ns:
        l = 2 * n - 1
        for q, m, e, k, t, R in lp_structures(l):
            if (n, q, e, k, t, -1) in done: continue
            base = dict(n=n, l=l, q=q, m=m, e=e, k=k, t=t, R=R)
            t0 = time.time()
            try:
                S = CP.Structure(l, q, e, skew=False, kgen=k); O = CL.Orbits(S, t)
                classes = LP.lp_stage1(S, O, a.maxopts)
            except (RuntimeError, MemoryError, AssertionError, ValueError) as ex:
                rec = dict(base, cls=-1, error=repr(ex)[:100], sec=round(time.time() - t0))
                print(json.dumps(rec), flush=True)
                with open(out, 'a') as f: f.write(json.dumps(rec) + '\n')
                continue
            rec = dict(base, cls=-1, classes=len(classes), sec=round(time.time() - t0))
            print(json.dumps(rec), flush=True)
            with open(out, 'a') as f: f.write(json.dumps(rec) + '\n')
            for ci, ys in enumerate(classes[:a.maxcls]):
                if (n, q, e, k, t, ci) in done: continue
                t1 = time.time()
                try:
                    P = LP.prepare_lp(S, O, ys, a.samples, np.random.default_rng(1))
                    sizes = P['sizes']; logN = sum(math.log10(x) for x in sizes)
                    est, rank, off = LP.gauss2(P['mus'], P['sigs'], P['t'], logN)
                    rec = dict(base, cls=ci, y=ys, sizes=[int(x) for x in sizes], logE=round(est, 2), rank=rank, off=round(off, 2))
                    if max(sizes) <= a.exact:
                        found = search_exact(S, O, P, np.random.default_rng(1))
                        rec.update(searched='exact', found=len(found))
                        for A, B in found:
                            hit = dict(l=l, order=2 * l + 2, q=q, e=e, kgen=k, twist=t, config=ci, method='lp_hunt exact', A=A, B=B)
                            with open(W.WORK / 'lp_hits.jsonl', 'a') as fh: fh.write(json.dumps(hit) + '\n')
                    else:
                        cs, Ws, condcov = LP.filters(P)
                        taus = LP.taus_for(P, [a.L, a.L], cs, Ws)
                        share = LP.kept_share(P, cs, Ws, condcov, taus)
                        rec.update(share=round(share, 4), logE_lists=round(est + math.log10(max(share, 1e-300)), 2))
                except (RuntimeError, MemoryError, AssertionError, ValueError, np.linalg.LinAlgError) as ex:
                    rec = dict(base, cls=ci, error=repr(ex)[:100])
                rec['sec'] = round(time.time() - t1)
                print(json.dumps(rec), flush=True)
                with open(out, 'a') as f: f.write(json.dumps(rec) + '\n')


if __name__ == '__main__':
    main()
