"""Validation battery: K-level exhaustive engine (gs_exh: gl1.c + gl2.c) versus the independent G-level brute force (gbf.count_G) on small shapes.
usage: validate.py MODE out.jsonl [tmax_per_shape] ; MODE in gs | spence | bord.  Resume-safe (shapes already in out.jsonl are skipped)."""
import sys, os, json, time, itertools, traceback
import numpy as np
from gsk import *
from gs_exh import exhaust
from gbf import count_G
from pk import signed_subgroups, KGrp


def shapes(mode):
    out = []
    Ks = [[2], [3], [4], [2, 2], [5], [6], [7], [8], [9], [3, 3], [10], [11], [13]]
    Tsets = {
        '[2]': [[(1,)], [(0,), (1,)]], '[3]': [[(0,), (1,), (2,)], [(1,), (2,)]], '[4]': [[(0,), (2,)], [(1,), (3,)], [(0,), (1,), (2,), (3,)]],
        '[2, 2]': [[(1, 0), (0, 1), (1, 1)], [(0, 0), (1, 0)]], '[5]': [[(0,), (1,), (2,), (3,), (4,)], [(1,), (2,), (3,), (4,)]], '[6]': [[(0,), (2,), (4,)], [(1,), (5,)], [(2,), (4,)]],
        '[7]': [[(0,), (1,), (2,), (3,), (4,), (5,), (6,)]], '[8]': [[(0,), (4,)], [(2,), (6,)]], '[9]': [[(0,), (3,), (6,)]], '[3, 3]': [[(0, 0), (1, 0), (2, 0)]],
        '[10]': [[(0,), (5,)]], '[11]': [], '[13]': []}
    for r in (3, 5, 7):
        for K in Ks:
            if mode == 'spence' and int(np.prod(K)) % 2: continue
            try: sub = signed_subgroups(K, maxgens=2)
            except Exception: sub = [([(1,) * len(K)], [])]
            for (G_, gens) in sub:
                hg = [(tuple(u), int(s)) for (u, s) in gens]
                if mode in ('gs', 'spence'): out.append(dict(mode=mode, K=K, r=r, hg=hg, T=None))
                else:
                    for T in Tsets.get(str(K), []): out.append(dict(mode='bord', K=K, r=r, hg=hg, T=[list(c) for c in T]))
    return out


def key(j): return json.dumps([j['mode'], j['K'], j['r'], j['hg'], j['T']])


def main():
    mode = sys.argv[1]; outp = sys.argv[2]; tmax = float(sys.argv[3]) if len(sys.argv) > 3 else 90.0
    want_match = int(sys.argv[4]) if len(sys.argv) > 4 else 30; want_nz = int(sys.argv[5]) if len(sys.argv) > 5 else 12
    done = set()
    if os.path.exists(outp):
        for ln in open(outp):
            try: done.add(key(json.loads(ln)['job']))
            except Exception: pass
    fo = open(outp, 'a')
    nmatch = 0; nnz = 0
    if os.path.exists(outp):
        for ln in open(outp):
            try:
                rr = json.loads(ln)
                if rr.get('status') == 'MATCH': nmatch += 1; nnz += (rr['G_count'] > 0)
            except Exception: pass
    jobs = sorted(shapes(mode), key=lambda j: (int(np.prod(j['K'])) * j['r'], len(j['hg'])))
    for job in jobs:
        if nmatch >= want_match and nnz >= want_nz: break
        if key(job) in done: continue
        hg = [(tuple(u), p) for u, p in job['hg']]
        rec = dict(job=job)
        try:
            sh = GShape(job['K'], job['r'], hg, job['mode'], [tuple(c) for c in job['T']] if job['T'] else None)
        except AssertionError as ex:
            rec['status'] = 'skip-invalid'; rec['why'] = str(ex)[:80]; fo.write(json.dumps(rec) + '\n'); fo.flush(); continue
        rec['nu'] = sh.nu; rec['Ntot'] = sh.N * sh.r
        if sh.nu > 7: rec['status'] = 'skip-big'; fo.write(json.dumps(rec) + '\n'); fo.flush(); continue
        t0 = time.time()
        try:
            gc, info = count_G(job['K'], job['r'], hg, job['mode'], [tuple(c) for c in job['T']] if job['T'] else None)
            rec['G_count'] = gc; rec['G_info'] = info; rec['t_G'] = round(time.time() - t0, 1)
        except Exception as ex:
            rec['status'] = 'skip-G'; rec['why'] = str(ex)[:80]; fo.write(json.dumps(rec) + '\n'); fo.flush(); continue
        t0 = time.time()
        try:
            tot, sols = exhaust(sh, tmax=tmax, max_inst=1_500_000, tl2=30.0)
            rec['K_count'] = tot['E_G']; rec['exhaustive'] = tot['exhaustive']; rec['N1_G'] = tot['N1_G']; rec['t_K'] = round(time.time() - t0, 1)
            rec['skipped'] = [(list(a), w) for a, w in tot['skipped']][:3]; rec['flags'] = tot['flags'][:3]
            if tot['exhaustive']:
                rec['status'] = 'MATCH' if tot['E_G'] == gc else 'MISMATCH'
                if rec['status'] == 'MATCH': nmatch += 1; nnz += (gc > 0)
            else: rec['status'] = 'partial'
        except Exception as ex:
            rec['status'] = 'error'; rec['why'] = repr(ex)[:200]
        fo.write(json.dumps(rec) + '\n'); fo.flush()
        print(rec['status'], job['K'], job['r'], job['hg'], job['T'], rec.get('G_count'), rec.get('K_count'), flush=True)

if __name__ == '__main__': main()
