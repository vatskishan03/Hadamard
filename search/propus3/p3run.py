"""Driver for p3.c: problem files, running the engine, exact re-checks of hits, solution orbits.

  library: write_problem / run_p3 / exact_check / orbit_group_images; the command lines are p3problem.py (problem file),
  p3sweep.py (one structure, all row patterns), p3cert.py (hits -> certificate); the test suites are p3test.py,
  p3validate.py, p3plant.py, p3scale_test.py, p3cert_test.py, p3resume_test.py, p3lowe.py.

A hit is two things: the three orbit masks printed by p3, and (here) the three +-1 sequences on Z_v, whose periodic
autocorrelations are summed directly (FFT, rounded: exact for these sizes) and compared with 4 v delta.  Nothing in this
file trusts the C engine's hash or class-value arithmetic.
"""
import itertools, json, math, os, subprocess, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from p3struct import Structure, paf

CC = ['cc', '-O3', '-march=native', '-ffast-math', '-pthread']


def build_binary(out=None, extra=()):
    out = Path(out or HERE / 'p3')
    src = HERE / 'p3.c'
    if not out.exists() or out.stat().st_mtime < src.stat().st_mtime:
        subprocess.run([*CC, *extra, str(src), '-o', str(out), '-lm'], check=True)
    return str(out)


# ------------------------------------------------------------------------------------------------------------------
def write_problem(st, rows, path, *, canonB=True, canonD=True, psd=(True, True, True), seed=1, tgt=None, energy=None):
    """Write the text problem file read by p3.c (format documented there and in README.md).
    energy: None, or dict(B=(tau, m, W), D=(tau, m, W)) with m (C,), W (C, C): keep a member iff |W (F - m)|^2 <= tau."""
    v, C = st.v, st.C
    ra, rb, rd = rows
    tgt = np.zeros(C, dtype=np.int64) if tgt is None else np.asarray(tgt, dtype=np.int64)
    rng = np.random.RandomState(seed)
    R = [int(x) | 1 for x in rng.randint(0, 2 ** 62, size=C, dtype=np.int64) * 4 + rng.randint(0, 4, size=C)]
    R = [(r * 0x9E3779B97F4A7C15 + 0x1234567) % (1 << 64) | 1 for r in R]
    L = []
    P = L.append
    P('P3PROB 1')
    P(f'{v} {C} {st.A.O} {st.B.O}')
    P(f'{ra} {rb} {rd}')
    P(' '.join(map(str, tgt.tolist())))
    P(' '.join(map(str, st.A.sizes.tolist())))
    P(' '.join(map(str, st.B.sizes.tolist())))
    for blk in (st.A, st.B):
        P(' '.join(map(str, blk.const0.tolist())))
        P(' '.join(map(str, blk.T.reshape(-1).tolist())))
    P(' '.join(repr(float(x)) for x in st.theta.reshape(-1)))
    bA = 4 * v if psd[0] else -1
    bB = 2 * v if psd[1] else -1
    bD = 4 * v if psd[2] else -1
    P(f'{bA} {bB} {bD}')
    perms = [st.perm_of(st.B, t).tolist() for t in st.mult_reps[1:]]
    P(f'{int(canonB)} {int(canonD)}')
    P(str(len(perms) if canonB else 0))
    if canonB:
        for p in perms:
            P(' '.join(map(str, p)))
    P(' '.join(map(str, st.perm_of(st.B, st.neg).tolist())))
    P(' '.join(map(str, R)))
    for key in ('B', 'D'):
        if energy and energy.get(key) is not None:
            tau, m, Wm = energy[key]
            P(f'1 {float(tau)!r}')
            P(' '.join(repr(float(x)) for x in m))
            P(' '.join(repr(float(x)) for x in np.asarray(Wm).reshape(-1)))
        else:
            P('0 0.0')
            P(' '.join('0.0' for _ in range(C)))
            P(' '.join('0.0' for _ in range(C * C)))
    Path(path).write_text('\n'.join(L) + '\n')


def run_p3(problem, binary=None, threads=1, ns=2, shard=None, ckpt=None, resume=False, hits=None, count=False,
           listonly=False, alist=None, blist=None, dlist=None, extra=(), timeout=None):
    binary = binary or build_binary()
    cmd = [binary, str(problem), '-t', str(threads), '-ns', str(ns)]
    if shard: cmd += ['-shard', f'{shard[0]}/{shard[1]}']
    if ckpt: cmd += ['-ckpt', str(ckpt)]
    if resume: cmd += ['-resume']
    if hits: cmd += ['-hits', str(hits)]
    if count: cmd += ['-count']
    if listonly: cmd += ['-listonly']
    for flag, f in (('-alist', alist), ('-blist', blist), ('-dlist', dlist)):
        if f: cmd += [flag, str(f)]
    cmd += list(extra)
    t = time.time()
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    if r.returncode != 0:
        raise RuntimeError(f'p3 failed ({r.returncode}): {r.stderr[-2000:]}')
    hitl, done, lists = [], None, None
    for line in r.stdout.splitlines():
        if line.startswith('HIT'):
            hitl.append(tuple(int(x, 16) for x in line.split()[1:4]))
        elif line.startswith('DONE'):
            done = dict(kv.split('=') for kv in line.split()[1:])
        elif line.startswith('LISTS'):
            lists = tuple(int(x) for x in line.split()[1:4])
    return dict(hits=hitl, done=done, lists=lists, stderr=r.stderr, wall=time.time() - t)


# ------------------------------------------------------------------------------------------------------------------
def masks_to_seqs(st, ma, mb, md):
    ya = np.array([1 if (ma >> i) & 1 else -1 for i in range(st.A.O)], dtype=np.int64)
    yb = np.array([1 if (mb >> i) & 1 else -1 for i in range(st.B.O)], dtype=np.int64)
    yd = np.array([1 if (md >> i) & 1 else -1 for i in range(st.B.O)], dtype=np.int64)
    return st.A.seq(ya), st.B.seq(yb), st.B.seq(yd)


def seqs_to_masks(st, A, B, D):
    """Inverse of masks_to_seqs (sequences constant on the orbits are assumed)."""
    ma = sum(1 << i for i, r in enumerate(st.A.reps) if A[r] > 0)
    mb = sum(1 << i for i, r in enumerate(st.B.reps) if B[r] > 0)
    md = sum(1 << i for i, r in enumerate(st.B.reps) if D[r] > 0)
    return ma, mb, md


def exact_check(st, A, B, D):
    """PAF_A + 2 PAF_B + PAF_D == 4 v delta on Z_v, A symmetric; direct, independent of the engine."""
    v = st.v
    A, B, D = (np.asarray(x, dtype=np.int64) for x in (A, B, D))
    if not np.array_equal(A, A[(-np.arange(v)) % v]):
        return False
    tot = paf(A) + 2 * paf(B) + paf(D)
    return bool(tot[0] == 4 * v and np.all(tot[1:] == 0))


# ------------------------------------------------------------------------------------------------------------------
def orbit_group_images(st, ma, mb, md):
    """Gamma-orbit of a solution in orbit-mask space: multipliers (all three blocks), reversal of B alone and of D alone,
    translations of B alone and of D alone by Fix(H).  Returns the set of mask triples (closure under the generators)."""
    v = st.v
    fix = [f for f in range(v) if all(h * f % v == f for h in st.H)]

    def pm(blk, perm, m):                       # bit i of image = bit perm[i] of m
        return sum(1 << i for i in range(blk.O) if (m >> int(perm[i])) & 1)

    gens = []
    for t in st.mult_reps[1:]:
        pa, pb = st.perm_of(st.A, t), st.perm_of(st.B, t)
        gens.append(lambda x, pa=pa, pb=pb: (pm(st.A, pa, x[0]), pm(st.B, pb, x[1]), pm(st.B, pb, x[2])))
    pneg = st.perm_of(st.B, st.neg)
    gens.append(lambda x: (x[0], pm(st.B, pneg, x[1]), x[2]))
    gens.append(lambda x: (x[0], x[1], pm(st.B, pneg, x[2])))
    for f in fix:
        if f == 0:
            continue
        ptr = np.array([st.B.orb[(r + f) % v] for r in st.B.reps])
        gens.append(lambda x, ptr=ptr: (x[0], pm(st.B, ptr, x[1]), x[2]))
        gens.append(lambda x, ptr=ptr: (x[0], x[1], pm(st.B, ptr, x[2])))
    seen = {(ma, mb, md)}
    frontier = [(ma, mb, md)]
    while frontier:
        nxt = []
        for x in frontier:
            for g in gens:
                y = g(x)
                if y not in seen:
                    seen.add(y)
                    nxt.append(y)
        frontier = nxt
    return seen


def gamma_order(st):
    """Generic orbit size |Mult| * 2 * |Fix(H)|^2 of the group of NOTES.md section 30 (reversal of B and D separately)."""
    fix = sum(1 for f in range(st.v) if all(h * f % st.v == f for h in st.H))
    return len(st.mult_reps) * 2 * fix * fix


def is_canonical(st, ma, mb, md):
    """B-mask minimal under Mult; D-mask minimal under reversal of D alone."""
    def pm(perm, m):
        return sum(1 << i for i in range(st.B.O) if (m >> int(perm[i])) & 1)
    for t in st.mult_reps[1:]:
        if pm(st.perm_of(st.B, t), mb) < mb:
            return False
    return pm(st.perm_of(st.B, st.neg), md) >= md
