#!/usr/bin/env python3
"""Exact joins of the jobs with ONE list too long for memory: that list is enumerated in 2^k slices (enum_dfs prefix
slices, which partition the list) and each slice is joined with the other three full lists, so nothing is sampled.
Takes the jobs of one N and mode that sched.py has not decided, whose largest list exceeds --maxslice and whose
second-largest is at most --second, in order of pair work up to --maxcost; each goes through the quotient screen first.

    python3 search/qhsweep/sliced.py N [--mode plain] [--maxcost 3e11] [--second 2e7] [--maxslice 1.2e7] [--budget SEC]
"""
import argparse, math, sys, time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import qhlib as Q
import qkill


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('N', type=int); ap.add_argument('--mode', default='plain', choices=Q.MODES)
    ap.add_argument('--maxcost', type=float, default=3e11); ap.add_argument('--second', type=float, default=2e7)
    ap.add_argument('--maxslice', type=float, default=1.2e7, help='largest slice of the big list'); ap.add_argument('--budget', type=float, default=1e9)
    a = ap.parse_args()
    N, mode = a.N, a.mode
    seen = Q.logged(mode, N)
    jobs = [j for j in Q.plan(mode, N, 64, 1e30, 1e30, 1e30) if j['z2'] and j['cost'] <= a.maxcost
            and max(j['sizes']) > a.maxslice and sorted(j['sizes'])[-2] <= a.second
            and seen.get((tuple(j['gens']), tuple(j['pat'])), {}).get('status') in (None, 'open')]
    print(f"N = {N} {mode}: {len(jobs)} single-big-list jobs within pair work {a.maxcost:g}, total {sum(j['cost'] for j in jobs):.3g}", flush=True)
    t0, cur, st = time.time(), None, None
    for j in jobs:
        if time.time() - t0 > a.budget:
            print('budget exhausted', flush=True); break
        killed, det = qkill.quotient_kill(mode, N, j['gens'], j['pat'])
        if killed:
            Q.record(mode, N, dict(j, status='quotient', detail=str(det)))
            print(f"{time.strftime('%H:%M:%S')} {Q.describe(j)}: QUOTIENT KILL {det}", flush=True); continue
        if tuple(j['gens']) != cur:
            if st is not None:
                Q.clear(st)
            st, cur = Q.struct_of(mode, N, j['gens']), tuple(j['gens'])
        nbits = max(1, math.ceil(math.log2(max(j['sizes']) / a.maxslice)))
        try:
            r = Q.run_job(j, st, nbits)
        except Exception as ex:
            print(f"{Q.describe(j)}: ERROR {ex!r}"[:300], flush=True); continue
        print(f"{time.strftime('%H:%M:%S')} {Q.describe(j)}: total {r['total']} in {r['slices']} slices enum {r['t_enum']}s join {r['t_join']}s"
              + (' HIT' if r.get('HIT') else ''), flush=True)
        if r.get('HIT'):
            print('HIT FOUND', flush=True); break
    if st is not None:
        Q.clear(st)


if __name__ == '__main__':
    main()
