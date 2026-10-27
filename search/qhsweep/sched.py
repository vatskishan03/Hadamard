#!/usr/bin/env python3
"""Cost-ordered exact sweep of the Szöllősi arrays of order 2N (qhlib.py: plain and Williamson-type arrays of
search/ssxx.py, bordered array of search/bsz.py) over several N at once.

Every job (multiplier group H with at most --maxo orbits, row-sum pattern) is taken in order of the pair work of its
exact join. Bordered jobs that fail the exact Z_2 test are dropped. From --screen pair work upwards a job first goes
through the exact quotient screen (qkill.py), whose kill is a final answer; the rest are enumerated (enum_dfs.c) and
joined (join4y.c), and a hit is verified exactly and written to hits_<mode>_<N>.jsonl (finish_hit.py takes it from
there). Resumable: every decided job is a line of WORK/qhsweep/<mode>_<N>.jsonl. --screen-only screens the jobs whose
lists are too long to enumerate and records the survivors as open.

    python3 search/qhsweep/sched.py plan N [N ...] [--modes plain,wt,bord] [--maxo 64] [--maxlist 2e7] [--show 25]
    python3 search/qhsweep/sched.py run N [N ...] [--modes plain,bord] [--cap 1e10] [--mincost 0] [--budget SEC]
                                        [--screen 1e8] [--screen-only]
    python3 search/qhsweep/sched.py status N [N ...] [--modes plain,wt,bord] [--maxo 64]      what the logs have decided
"""
import argparse, sys, time
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import qhlib as Q
import qkill


def show_plan(Ns, modes, maxo, maxlist, show):
    for N in Ns:
        for mode in modes:
            jobs = Q.plan(mode, N, maxo, maxlist, 1e300)
            live = [j for j in jobs if j['z2']]
            print(f"N = {N} {mode}: {len(jobs)} structure/pattern pairs, {len(live)} pass the Z_2 test, "
                  f"{sum(j['fits'] for j in live)} with lists within {maxlist:.3g}")
            for j in live[:show]:
                print(('  ' if j['fits'] else 'x ') + Q.describe(j))


def status(Ns, modes, maxo):
    for N in Ns:
        for mode in modes:
            seen = Q.logged(mode, N)
            jobs = [j for j in Q.plan(mode, N, maxo, 1e30, 1e300) if j['z2']]
            c = Counter(seen.get((tuple(j['gens']), tuple(j['pat'])), {}).get('status', 'not run') for j in jobs)
            print(f"N = {N} {mode}: {len(jobs)} jobs{' past the Z_2 test' if mode == 'bord' else ''}: {c['join']} joined, "
                  f"{c['quotient']} quotient kills, {c['open']} open after the screen, {c['not run']} not run; "
                  f"hits {sum(bool(r.get('HIT')) for r in seen.values())}")


def run(Ns, modes, cap, mincost, budget, screen, screen_only, maxo, maxlist):
    jobs = []
    for N in Ns:
        for mode in modes:
            seen = Q.logged(mode, N)
            for j in Q.plan(mode, N, maxo, maxlist, 1e300 if screen_only else cap):
                key = (tuple(j['gens']), tuple(j['pat']))
                if screen_only:
                    take = not j['fits'] and key not in seen
                else:
                    take = j['ok'] and j['cost'] >= mincost and seen.get(key, {}).get('status') in (None, 'open')
                if take and j['z2']:
                    jobs.append(j)
    jobs.sort(key=lambda j: (j['cost'], j['N'], -j['H']))
    print(f"{len(jobs)} jobs, total pair work {sum(j['cost'] for j in jobs):.3g}", flush=True)
    t0, cur, st = time.time(), None, None
    for j in jobs:
        if time.time() - t0 > budget:
            print('budget exhausted', flush=True); break
        mode, N = j['mode'], j['N']
        if screen_only or j['cost'] >= screen:
            tq = time.time()
            killed, det = qkill.quotient_kill(mode, N, j['gens'], j['pat'])
            if killed or screen_only:
                Q.record(mode, N, dict(j, status='quotient' if killed else 'open', detail=str(det), t_join=round(time.time() - tq, 1)))
                print(f"{time.strftime('%H:%M:%S')} {Q.describe(j)}: {'QUOTIENT KILL' if killed else 'SURVIVES'} {det}", flush=True)
                continue
        if (mode, N, tuple(j['gens'])) != cur:
            if st is not None:
                Q.clear(st)
            st, cur = Q.struct_of(mode, N, j['gens']), (mode, N, tuple(j['gens']))
        try:
            r = Q.run_job(j, st)
        except Exception as ex:
            print(f"{time.strftime('%H:%M:%S')} {Q.describe(j)}: ERROR {ex!r}"[:300], flush=True); continue
        Q.clear(st, 40e6)
        print(f"{time.strftime('%H:%M:%S')} {Q.describe(j)}: total {r['total']} enum {r['t_enum']}s join {r['t_join']}s"
              + (' HIT' if r.get('HIT') else ''), flush=True)
        if r.get('HIT'):
            print('HIT FOUND', flush=True); break
    if st is not None:
        Q.clear(st)
    print('sched finished', flush=True)


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest='cmd', required=True)
    for name in ('plan', 'run', 'status'):
        p = sub.add_parser(name); p.add_argument('N', type=int, nargs='+')
        p.add_argument('--modes', default='plain,bord' if name == 'run' else 'plain,wt,bord')
        p.add_argument('--maxo', type=int, default=64, help='most orbits of H')
        if name != 'status':
            p.add_argument('--maxlist', type=float, default=2e7, help='longest list that is enumerated')
    sub.choices['plan'].add_argument('--show', type=int, default=25)
    p = sub.choices['run']
    p.add_argument('--cap', type=float, default=1e10, help='largest pair work of a job')
    p.add_argument('--mincost', type=float, default=0); p.add_argument('--budget', type=float, default=1e9)
    p.add_argument('--screen', type=float, default=1e8, help='pair work from which the quotient screen runs first')
    p.add_argument('--screen-only', action='store_true', help='quotient-screen the jobs with lists too long, no joins')
    a = ap.parse_args()
    modes = a.modes.split(',')
    if a.cmd == 'plan':
        show_plan(a.N, modes, a.maxo, a.maxlist, a.show)
    elif a.cmd == 'status':
        status(a.N, modes, a.maxo)
    else:
        run(a.N, modes, a.cap, a.mincost, a.budget, a.screen, a.screen_only, a.maxo, a.maxlist)


if __name__ == '__main__':
    main()
