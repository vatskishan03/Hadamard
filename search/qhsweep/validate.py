#!/usr/bin/env python3
"""Validation of the qhsweep tools. Every check is exact and ends in OK or FAILED (exit status 1); the numbers are
the ones of README.md.

  enum      enum_dfs.c against search/bnw_enum.c: identical sets of (sign mask, class values) for plain lists and for
            bordered lists with y_0 fixed
  join      join4y.c (raw, lattice-normalised, normalised with the forced out-of-memory split) against search/join4x.c
            on random lattice instances and the five known structures, qhlib.solve and the exact verification of its
            first solution
  z2        z2fast.py against bsz.z2_ok
  screen    a quotient kill (qkill.py) always agrees with the exact join: total 0; known solutions are never killed
  slice     the prefix slices partition a list and sliced joins give the same totals
  pipeline  sched.py -> finish_hit.py end to end: N = 779 (q = 1559) and N = 11 in every mode, sliced.py on N = 163

    python3 search/qhsweep/validate.py [enum] [join] [z2] [screen] [slice] [pipeline]     (all when none is named)
"""
import json, math, os, struct, subprocess, sys, time
from math import gcd
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import qhlib as Q
import qkill
import bnw, bsz, ssxx, workdir as W
from bnw import Struct, subgroups
from z2fast import z2_ok_fast

HERE = Path(__file__).resolve().parent
KNOWN = [(779, (49,), (5, 7, 39, 39), 96), (371, (13,), (3, 3, 25, 29), 96), (143, (3,), (3, 9, 11, 19), 16),
         (163, (38,), (17, 19, 1, 1), 8352), (119, (8,), (1, 5, 15, 15), 4704)]       # N, H, rows, number of solutions
failed = []


def report(name, ok, detail):
    print(f"{name}: {'OK' if ok else 'FAILED'} ({detail})", flush=True)
    if not ok:
        failed.append(name)


def job_of(mode, st, pat):
    cnt = {k: st.counts(k) for k in 'AL'}
    sizes = [cnt[k].get(r, 0) // (2 if mode == 'bord' and r == 0 else 1) for k, r in zip(ssxx.KINDS, pat)]
    return dict(N=st.N, mode=mode, gens=list(st.gens), pat=list(pat), sizes=sizes)


# ---------------------------------------------------------------- enum

def same_lists(st, kind, r, fix0):
    """(number of rows, equal): enum_dfs against bnw_enum on one list"""
    d = Path(W.tmpdir())
    Q.write_input(st, kind, r, d / 'val.in', fix0)
    out = {}
    for tag, exe in (('a', W.bin('bnw_enum')), ('b', Q.build('enum_dfs'))):
        subprocess.run([exe, str(d / 'val.in'), str(d / f'val_{tag}.v'), str(d / f'val_{tag}.m')], check=True, capture_output=True)
        V = np.fromfile(d / f'val_{tag}.v', np.int16).reshape(-1, st.C); M = np.fromfile(d / f'val_{tag}.m', np.uint64)
        o = np.argsort(M); out[tag] = (M[o], V[o])
    for p in d.glob('val*'):
        p.unlink()
    n = len(out['a'][0])
    return n, n == len(out['b'][0]) and np.array_equal(out['a'][0], out['b'][0]) and np.array_equal(out['a'][1], out['b'][1])


def check_enum():
    plain = [(143, (3,)), (247, (9, 35)), (119, (8,)), (163, (38,)), (185, (9,)), (93, (4,)), (267, (4,)), (377, (9,)), (183, (34,))]
    bord = [(N, gens) for N in (19, 23, 27, 31, 35) for H, gens in subgroups(N - 1).items()
            if sum(gcd(h - 1, N - 1) for h in H) // len(H) <= 18]
    n = {'plain': [0, 0], 'bord': [0, 0], 'fix0': 0}                  # lists, structures
    bad = 0
    for mode, cases in (('plain', plain), ('bord', bord)):
        for N, gens in cases:
            st = Struct(N, gens, plain=mode == 'plain')
            used = 0
            for kind in 'LA':
                if len(st.sizes(kind)) > 22:
                    continue
                cnt = st.counts(kind)
                for r in sorted(cnt):
                    if not 0 <= r <= (6 if mode == 'bord' else 2 * math.isqrt(N)) or cnt[r] > 3e6 or (r == 0 and mode == 'plain'):
                        continue
                    fix0 = int(mode == 'bord' and r == 0)
                    rows, ok = same_lists(st, kind, r, fix0)
                    if not ok or rows != cnt[r] // (1 + fix0):
                        bad += 1; print(f'  MISMATCH {mode} N={N} H={gens} {kind} r={r}: {rows} rows, {cnt[r]} expected', flush=True)
                    n[mode][0] += 1; n['fix0'] += fix0; used += 1
            n[mode][1] += bool(used)
    report('enum', bad == 0, f"{n['plain'][0]} plain lists of {n['plain'][1]} structures and {n['bord'][0]} bordered lists of {n['bord'][1]} structures "
           f"({n['fix0']} with y_0 fixed) identical to bnw_enum, masks and class values, row counts equal Struct.counts")


def check_partition():
    """the 2^k prefix slices of a list partition it; the bordered r = 0 list has y_0 fixed"""
    cases = [('plain', 779, (49,), 'A', 39), ('plain', 371, (13,), 'A', 25), ('plain', 163, (38,), 'A', 19), ('bord', 31, subgroups(30)[frozenset({1})], 'L', 0)]
    bad, tot = 0, 0
    for mode, N, gens, kind, r in cases:
        st = Q.struct_of(mode, N, gens)
        V, M = Q.enum_list(st, kind, r)
        for k in (1, 2, 3):
            parts = [Q.enum_list(st, kind, r, k, j) for j in range(2 ** k)]
            allm = np.concatenate([np.asarray(p[1]) for p in parts]); allv = np.concatenate([np.asarray(p[0]) for p in parts])
            o1, o2 = np.argsort(M), np.argsort(allm)
            ok = len(allm) == len(M) and np.array_equal(np.asarray(M)[o1], allm[o2]) and np.array_equal(np.asarray(V)[o1], allv[o2])
            bad += not ok; tot += 1
            print(f"  {mode} N={N} H={gens} {kind} r={r}: {len(M)} rows, k={k}, slices {[len(p[1]) for p in parts]} {'OK' if ok else 'MISMATCH'}", flush=True)
        Q.clear(st)
    report('slice partition', bad == 0, f'{tot} partitions of {len(cases)} lists into 2, 4, 8 slices')


# ---------------------------------------------------------------- join

def forced_split(L, C):
    """join with a MemoryError forced on the first level: the halves must add up"""
    orig, state = Q._join, {'n': 0}

    def fake(*a):
        state['n'] += 1
        if state['n'] == 1:
            raise MemoryError('forced')
        return orig(*a)
    Q._join = fake
    try:
        return int(Q.join(L, C, stop=False)[1])
    finally:
        Q._join = orig


def raw_join(exe, lists, C, T):
    """total of a join binary (join4x or join4y) on the raw lists and the target vector T"""
    d = Path(W.tmpdir()); path = d / 'val_join.bin'
    with open(path, 'wb') as f:
        f.write(struct.pack('5i', *[len(x) for x in lists], C)); f.write(np.asarray(T, np.int16).tobytes())
        for x in lists:
            f.write(np.ascontiguousarray(x, np.int16).tobytes())
    out = subprocess.run([exe, str(path), '2', '8000000', '0'], capture_output=True, text=True, check=True).stdout
    path.unlink()
    return int([ln.split()[1] for ln in out.splitlines() if ln.startswith('TOTAL')][0])


def check_random(trials=60, seed=1):
    """random instances on a lattice (values g * v + offset, class-dependent g), targets that are hit by construction or
    missed by one: join4x, join4y raw and normalised give the same total"""
    rng = np.random.default_rng(seed)
    bad = hits = 0
    for t in range(trials):
        C = int(rng.integers(2, 7)); n = [int(rng.integers(5, 300)) for _ in range(4)]
        g, off = rng.integers(1, 4, C), rng.integers(-3, 4, (4, C))
        lists = [(rng.integers(-4, 5, (n[k], C)) * g + off[k]).astype(np.int16) for k in range(4)]
        T = sum(L[rng.integers(len(L))].astype(np.int64) for L in lists)
        if t % 3 == 2:
            T[rng.integers(C)] += 1
        x = raw_join(W.bin('join4x'), lists, C, T)
        y = raw_join(Q.build('join4y'), lists, C, T)
        z = int(Q.join(lists, C, stop=False, target=T)[1])
        bad += not x == y == z; hits += x > 0
        if not x == y == z:
            print(f'  MISMATCH trial {t}: join4x {x}, join4y {y}, normalised {z}', flush=True)
    report('join random', bad == 0, f'{trials} random lattice instances ({hits} with solutions): join4x, join4y and normalised join4y give the same totals')


def check_join():
    check_random()
    bad = 0
    for N, gens, pat, want in KNOWN:
        st = Struct(N, gens, plain=True)
        job = job_of('plain', st, pat)
        lists = [Q.enum_list(st, k, r)[0] for k, r in zip(ssxx.KINDS, pat)]
        order = ssxx.best_pairing([len(x) for x in lists]); L = [lists[k] for k in order]
        t = time.time(); x = bnw.join4(L, st.C, K=2, target=0)[1]; tx = time.time() - t
        t = time.time(); y = int(Q.join(L, st.C, stop=False, norm=False)[1]); ty = time.time() - t
        t = time.time(); z = int(Q.join(L, st.C, stop=False)[1]); tz = time.time() - t
        f = forced_split(L, st.C)
        total, hit, _, _ = Q.solve(job, st, stop=False)
        Q.verify('plain', st, hit)                                    # the first solution is a complex Hadamard matrix
        ok = x == y == z == f == total == want
        bad += not ok
        print(f"  N={N} H={gens} rows {pat}: lists {job['sizes']}, join4x {x} ({tx:.1f}s), join4y raw {y} ({ty:.1f}s), "
              f"normalised {z} ({tz:.1f}s), forced split {f}, solve {total}, expected {want} {'OK' if ok else 'MISMATCH'}", flush=True)
        Q.clear(st)
    report('join', bad == 0, f"totals {', '.join(str(w) for *_, w in KNOWN)} on N = {', '.join(str(c[0]) for c in KNOWN)} for join4x, join4y, normalised join4y, forced split and solve; first solutions verified")


# ---------------------------------------------------------------- z2

def check_z2():
    n = bad = 0
    for N in (23, 31, 47, 59, 71, 87, 91, 95, 111, 119, 143, 155, 183, 227):
        for H, gens in subgroups(N - 1).items():
            if sum(gcd(h - 1, N - 1) for h in H) // len(H) > 22:
                continue
            st = Struct(N, gens, plain=False)
            for _, pat in bsz.patterns(st):
                n += 1
                if bsz.z2_ok(st, pat) != z2_ok_fast(st, pat):
                    bad += 1; print(f'  MISMATCH N={N} H={gens} rows {pat}', flush=True)
    report('z2', bad == 0, f'z2fast agrees with bsz.z2_ok on {n} (structure, pattern) cases')


# ---------------------------------------------------------------- screen

def check_screen(per_n=24):
    bad = kills = surv = tot0 = 0
    for N, gens, pat, want in KNOWN:                                   # the known solutions are never killed
        killed, det = qkill.quotient_kill('plain', N, gens, pat)
        print(f"  N={N} H={gens} rows {pat} ({want} solutions): {'KILLED' if killed else 'survives'} {det}", flush=True)
        bad += killed
    for N in (119, 143, 163, 371):
        jobs = [j for j in Q.plan('plain', N, 64, 2e7, 1e7) if j['ok']]
        sample = jobs[::max(1, len(jobs) // per_n)][:per_n]
        k = s = 0
        for j in sample:
            killed, det = qkill.quotient_kill('plain', N, j['gens'], j['pat'])
            st = Struct(N, j['gens'], plain=True)
            total = Q.solve(j, st, stop=False)[0]
            Q.clear(st)
            if killed and total:
                bad += 1; print(f"  FAILED: killed job with {total} solutions: {Q.describe(j)}", flush=True)
            k += killed; s += not killed; tot0 += (not killed) and total == 0
            if not killed:
                print(f"  survivor: {Q.describe(j)}: {total} solutions", flush=True)
        kills += k; surv += s
        print(f'  N={N}: {len(sample)} of {len(jobs)} jobs within pair work 1e7: {k} killed (total 0), {s} survive', flush=True)
    report('screen', bad == 0, f'no kill of a job with solutions: {kills} kills agree with total 0, {surv} survivors ({tot0} with total 0)')


# ---------------------------------------------------------------- slice

def check_slice():
    bad = 0
    for N, gens, pat, want in KNOWN[3:] + [(143, (3,), (3, 9, 11, 19), 16)]:
        st = Struct(N, gens, plain=True)
        job = job_of('plain', st, pat)
        for nbits in (1, 2, 3):
            total = Q.solve(job, st, nbits, stop=False)[0]
            bad += total != want
            print(f"  N={N} H={gens} rows {pat}: sizes {job['sizes']}, 2^{nbits} slices of the largest list: total {total}, expected {want}", flush=True)
        total, hit, _, _ = Q.solve(job, st, 2, stop=True)
        Q.verify('plain', st, hit)
        Q.clear(st)
    check_partition()
    report('slice', bad == 0, 'sliced totals equal the unsliced ones for 2, 4, 8 slices; the stop mode finds a verified solution')


# ---------------------------------------------------------------- pipeline

def run(cmd, work):
    env = dict(os.environ, HADAMARD_WORK=str(work))
    t = time.time()
    r = subprocess.run([sys.executable, *cmd], capture_output=True, text=True, env=env)
    lines = (r.stdout + r.stderr).strip().splitlines()
    print(f"  $ {' '.join(Path(c).name if c.endswith('.py') else c for c in cmd)}   [{time.time() - t:.0f}s]", flush=True)
    for ln in lines[-4:]:
        print('    ' + ln, flush=True)
    return r.returncode, '\n'.join(lines)


def pipeline(mode, N, extra, check_repo=False):
    work = W.WORK / f'pipeline_{mode}_{N}'
    rc1, o1 = run([str(HERE / 'sched.py'), 'run', str(N), '--modes', mode, *extra], work)
    hits = work / 'qhsweep' / f'hits_{mode}_{N}.jsonl'
    ok = rc1 == 0 and 'HIT FOUND' in o1 and hits.exists()
    if ok:
        rc2, o2 = run([str(HERE / 'finish_hit.py'), str(hits)], work)
        ok = rc2 == 0 and f'K K* = {2 * N} I exact' in o2 and 'holds' in o2
        if check_repo:
            ref = Q.ROOT / 'complex-hadamard' / f'C{2 * N}.json'
            mine = json.loads((work / 'qhsweep' / 'cert' / f'C{2 * N}.json').read_text()); old = json.loads(ref.read_text())
            print(f"    repository certificate {ref.name}: same hit {mine['sha256'] == old['sha256']} (rows {mine['rowsums']} and {old['rowsums']}, H = <{mine['multipliers'][0]}> and <{old['multipliers'][0]}>)", flush=True)
    for p in (work / 'qhsweep' / 'lists').glob('*'):
        p.unlink()
    return ok


def check_pipeline():
    ok = pipeline('plain', 779, ['--cap', '1e7'], check_repo=True)
    report('pipeline plain N=779', ok, 'sched.py finds a hit, finish_hit.py: certificate, K K* = 1558 I exact, complex_miyamoto(q = 1559) of order 6236 with H H* = 6236 I')
    for mode in Q.MODES:
        ok = pipeline(mode, 11, [])
        report(f'pipeline {mode} N=11', ok, 'sched.py finds a hit, finish_hit.py: K K* = 22 I exact, complex_miyamoto(q = 23) of order 92')
    work = W.WORK / 'pipeline_sliced_163'
    rc, out = run([str(HERE / 'sliced.py'), '163', '--maxslice', '20000', '--maxcost', '1e8'], work)
    report('pipeline sliced N=163', rc == 0 and 'HIT FOUND' in out and '4 slices' in out, 'sliced.py finds a verified hit in 4 slices')


if __name__ == '__main__':
    todo = sys.argv[1:] or ['enum', 'join', 'z2', 'screen', 'slice', 'pipeline']
    for name in todo:
        globals()['check_' + name]()
    print('ALL OK' if not failed else f'FAILED: {failed}')
    sys.exit(1 if failed else 0)
