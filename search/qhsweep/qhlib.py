"""Shared parts of the qhsweep tools: the C helpers, the exact lists and join of a job, the job plan, the logs.

A job is a multiplier group H and a row-sum pattern of one of three arrays of complex Hadamard matrices of order 2N:
  plain  Szöllősi's array (search/ssxx.py), sequences on Z_N, -1 not in H, target 0;
  wt     the same with -1 in H (all four blocks symmetric, Williamson type), target 0;
  bord   the bordered array (search/bsz.py), sequences on Z_(N-1), target -4.
Its four lists (every sign vector constant on the orbits of H with the given row sum) are enumerated completely
(enum_dfs.c) and joined exactly (join4y.c, after an exact lattice normalisation of every class); a hit is verified
exactly (ssxx.complex_hadamard, bsz.build). Files go below WORK/qhsweep (WORK is $HADAMARD_WORK, search/workdir.py):
bin/, lists/ (cached lists, deleted as a sweep moves on), <mode>_<N>.jsonl (a line per decided job, status join,
quotient or open) and hits_<mode>_<N>.jsonl.
"""
import json, math, os, struct, subprocess, sys, tempfile, time
from math import gcd
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = Path(__file__).resolve().parents[2]                    # repository root
sys.path[:0] = [str(HERE), str(ROOT / 'search')]
import workdir as W
import bsz, ssxx
from bnw import Struct, subgroups, mask_to_y
from z2fast import z2_ok_fast

QW = W.WORK / 'qhsweep'
QW.mkdir(exist_ok=True)
MODES = ('plain', 'wt', 'bord')
FLAGS = {'enum_dfs': ['-O3', '-march=native'], 'join4y': ['-O3', '-march=native', '-fopenmp']}
CH = 1 << 21                                                  # rows per chunk when a list is scanned


def build(name):
    """compile qhsweep/<name>.c into WORK/qhsweep/bin on first use (or when the source is newer); the binary's path"""
    src, out = HERE / f'{name}.c', QW / 'bin' / name
    out.parent.mkdir(exist_ok=True)
    if not out.exists() or out.stat().st_mtime < src.stat().st_mtime:
        tmp = out.with_name(f'.{name}.{os.getpid()}')         # compile aside, then rename: a running copy keeps its inode
        subprocess.run(['cc', *FLAGS[name], str(src), '-o', str(tmp)], check=True)
        os.replace(tmp, out)
    return str(out)


def struct_of(mode, N, gens):
    return Struct(N, gens, plain=mode != 'bord')


def describe(j):
    return f"{j['mode']} N={j['N']} |H|={j['H']} H={tuple(j['gens'])} OA={j['OA']} rows {tuple(j['pat'])} sizes {j['sizes']} cost {j['cost']:.3g}"


# ---------------------------------------------------------------- plan

def plan(mode, N, maxo=64, maxlist=2e7, maxpairs=1e12, maxmem=1.8e9):
    """Every (H, row sums) job of the mode with at most maxo orbits of H, cheapest pair work first (the best pairing of
    ssxx.pair_cost). fits: the lists and their memory are within the limits; ok: fits and the pair work is within
    maxpairs; z2: the exact Z_2 compression test of the bordered array (z2fast.py), True for the other modes."""
    m = N - 1 if mode == 'bord' else N
    out = []
    for H, gens in subgroups(m).items():
        if mode != 'bord' and ((N - 1) in H) != (mode == 'wt'):
            continue
        if sum(gcd(h - 1, m) for h in H) // len(H) > maxo:        # orbit count by Burnside
            continue
        st = struct_of(mode, N, gens)
        cnt = {k: st.counts(k) for k in 'AL'}
        pats = [p for _, p in bsz.patterns(st)] if mode == 'bord' else ssxx.patterns(st)
        for pat in pats:
            sizes = [cnt[k].get(r, 0) // (2 if mode == 'bord' and r == 0 else 1) for k, r in zip(ssxx.KINDS, pat)]
            if min(sizes) == 0:
                continue
            cost = ssxx.pair_cost(sizes, ssxx.best_pairing(sizes))
            mem = sum(sizes) * (2 * st.C + 36)
            fits = max(sizes) <= maxlist and mem <= maxmem
            out.append(dict(N=N, mode=mode, gens=list(gens), H=len(H), OA=len(st.sizesA), OL=len(st.sizesL), C=st.C,
                            pat=[int(x) for x in pat], sizes=[int(x) for x in sizes], cost=float(cost), mem=float(mem),
                            fits=bool(fits), ok=bool(fits and cost <= maxpairs), z2=mode != 'bord' or bool(z2_ok_fast(st, pat))))
    out.sort(key=lambda r: (r['cost'], -r['H']))
    return out


# ---------------------------------------------------------------- lists

def write_input(st, kind, r, path, fix0=0):
    sz = st.sizes(kind)
    with open(path, 'wb') as f:
        f.write(struct.pack('4i', len(sz), st.C, r, fix0))
        f.write(np.asarray(sz, np.int32).tobytes()); f.write(np.asarray(st.S[kind], np.int32).tobytes())


def mm(path, dtype, C):
    """read-only memmap of a binary list file (an empty array for an empty file)"""
    n = os.path.getsize(path) // np.dtype(dtype).itemsize
    if C:
        n //= C
    if n == 0:
        return np.zeros((0, C) if C else (0,), dtype)
    return np.memmap(path, dtype=dtype, mode='r', shape=(n, C) if C else (n,))


def stem(st, kind, r, nbits=0, j=0):
    return f"{'p' if st.plain else 'b'}{st.N}_{'-'.join(map(str, st.gens))}_{kind}_{r}" + (f'_k{nbits}j{j}' if nbits else '')


def enum_list(st, kind, r, nbits=0, j=0):
    """(class values int16 [n, C], sign masks uint64 [n]) of every sign vector of a block with row sum r, cached and
    memory-mapped. Bordered with r = 0: y_0 = +1 (the global sign). nbits > 0: slice j of 2^nbits, the first nbits
    decisions of the DFS (largest orbits first) are the bits of j."""
    d = QW / 'lists'
    d.mkdir(exist_ok=True)
    base = d / stem(st, kind, r, nbits, j)
    vals, masks = base.with_suffix('.v'), base.with_suffix('.m')
    if not vals.exists():
        inp = base.with_suffix('.in')
        write_input(st, kind, r, inp, int(not st.plain and r == 0))
        subprocess.run([build('enum_dfs'), str(inp), f'{vals}.tmp', f'{masks}.tmp', '-1', str(nbits), str(j)],
                       check=True, capture_output=True)
        os.replace(f'{vals}.tmp', vals); os.replace(f'{masks}.tmp', masks)
        inp.unlink()
    return mm(vals, np.int16, st.C), mm(masks, np.uint64, 0)


def clear(st, above=0):
    """delete the cached lists of a structure (those of at least `above` bytes: the disk is shared and small)"""
    for p in (QW / 'lists').glob(f"{'p' if st.plain else 'b'}{st.N}_{'-'.join(map(str, st.gens))}_*"):
        if p.stat().st_size >= above:
            p.unlink()


# ---------------------------------------------------------------- join

def write_join_input(path, lists, C, target, norm=True):
    """Write the join4y input. With norm every class is reduced to its coset lattice (exact): the values of a list on a
    class lie on min + gZ, g the gcd of the differences, so subtracting the minima and dividing by the gcd of all
    four lists keeps the solutions. The lists are scanned in chunks. False (nothing written) when the lattice test
    shows there is no solution: a class constant on all lists that misses the target, or a target off the lattice."""
    T = np.array(np.broadcast_to(target, C), np.int64)            # the target: one number or one per class
    mins, gs = [], np.zeros(C, np.int64)
    if norm:
        for x in lists:
            lo = np.full(C, 1 << 30, np.int64); hi = np.full(C, -(1 << 30), np.int64)
            for a in range(0, len(x), CH):
                blk = np.asarray(x[a:a + CH])
                lo = np.minimum(lo, blk.min(0)); hi = np.maximum(hi, blk.max(0))
            mins.append(lo)
            pres = [np.zeros(int(hi[c] - lo[c]) + 1, bool) for c in range(C)]     # the distinct values per class
            for a in range(0, len(x), CH):
                blk = np.asarray(x[a:a + CH]).astype(np.int32) - lo.astype(np.int32)
                for c in range(C):
                    pres[c][blk[:, c]] = True
            for c in range(C):
                u = np.flatnonzero(pres[c])
                if len(u) > 1:
                    gs[c] = math.gcd(int(gs[c]), int(np.gcd.reduce(u)))
        base = np.sum(mins, axis=0)
        for c in range(C):
            g, rem = int(gs[c]), int(T[c] - base[c])
            if rem % g if g else rem:
                return False
            T[c] = rem // g if g else 0
    with open(path, 'wb') as f:
        f.write(struct.pack('5i', *[len(x) for x in lists], C))
        f.write(np.asarray(T, np.int16).tobytes())
        for k, x in enumerate(lists):
            for a in range(0, len(x), CH):
                ch = np.asarray(x[a:a + CH]).astype(np.int64)
                if norm:
                    ch = ch - mins[k]
                    nz = gs > 0
                    ch[:, nz] //= gs[nz]
                    ch[:, ~nz] = 0
                f.write(ch.astype(np.int16).tobytes())
    return True


def _join(lists, C, K, budget, stop, target, norm):
    fd, path = tempfile.mkstemp(suffix='.bin', dir=W.tmpdir()); os.close(fd)
    try:
        if not write_join_input(path, lists, C, target, norm):
            return [], '0'
        r = subprocess.run([build('join4y'), path, str(K), str(budget), '1' if stop else '0'], capture_output=True, text=True)
    finally:
        os.unlink(path)
    if r.returncode != 0:
        if 'out of memory' in r.stderr or r.returncode in (-9, -11, 137, 139):
            raise MemoryError(r.stderr[-200:])
        raise RuntimeError(f'join4y failed rc={r.returncode}: {r.stderr[-300:]}')
    sols = [tuple(map(int, ln.split()[1:5])) for ln in r.stdout.splitlines() if ln.startswith('SOL')]
    tot = [ln.split(None, 1)[1] for ln in r.stdout.splitlines() if ln.startswith('TOTAL')]
    return sols, (tot[0] if tot else '0')


def join(lists, C, K=2, budget=8000000, stop=True, target=0, norm=True, depth=0):
    """Exact join: (index tuples of the first K solutions, TOTAL as join4y prints it: the number of solutions, or
    '>= K (stopped)' with stop). When join4y runs out of memory the larger list of its build pair is halved and the
    halves are joined separately (up to four levels)."""
    if any(len(x) == 0 for x in lists):
        return [], '0'
    try:
        return _join(lists, C, K, budget, stop, target, norm)
    except MemoryError:
        sizes = [len(x) for x in lists]
        bp = (0, 1) if sizes[0] * sizes[1] <= sizes[2] * sizes[3] else (2, 3)
        k = bp[0] if sizes[bp[0]] >= sizes[bp[1]] else bp[1]
        if depth >= 4 or sizes[k] < 2:
            raise
        h = sizes[k] // 2
        sols_all, tots = [], 0
        for off, part in ((0, slice(0, h)), (h, slice(h, None))):
            sub = list(lists); sub[k] = lists[k][part]
            sols, tot = join(sub, C, K, budget, stop, target, norm, depth + 1)
            sols_all += [tuple(v + off if i == k else v for i, v in enumerate(s)) for s in sols]
            tots += int(tot) if tot.isdigit() else 0
            if sols and stop:
                break
        return sols_all, str(tots)


# ---------------------------------------------------------------- jobs

def solve(job, st, nbits=0, stop=True):
    """Exact join of a job: (total, hit, seconds enumerating, seconds joining). hit = the four sequences of the first
    solution found, None if there is none; with stop the join ends at the first solution and total is None. nbits > 0:
    the largest list (the one too big for memory) is enumerated in 2^nbits slices, each joined with the other three
    lists and deleted."""
    kinds, pat = ssxx.KINDS, job['pat']
    kb = int(np.argmax(job['sizes'])) if nbits else -1
    fixed = {k: enum_list(st, kinds[k], pat[k]) for k in range(4) if k != kb}
    total, hit, te, tj = 0, None, 0.0, 0.0
    for j in range(2 ** nbits):
        t = time.time()
        L = [enum_list(st, kinds[k], pat[k], nbits, j) if k == kb else fixed[k] for k in range(4)]
        te += time.time() - t
        sizes, tot = [len(V) for V, _ in L], '0'
        if min(sizes) > 0:
            order = ssxx.best_pairing(sizes)
            t = time.time()
            sols, tot = join([L[k][0] for k in order], st.C, stop=stop, target=st.target)
            tj += time.time() - t
            if sols and hit is None:
                idx = dict(zip(order, sols[0]))
                hit = [st.expand(kinds[k], mask_to_y(L[k][1][idx[k]], len(st.sizes(kinds[k])))) for k in range(4)]
        if nbits:
            for p in (QW / 'lists').glob(stem(st, kinds[kb], pat[kb], nbits, j) + '.*'):
                p.unlink()
        if hit is not None and stop:
            return None, hit, te, tj
        total += int(tot)
    return total, hit, te, tj


def logged(mode, N):
    """{(gens, rows): record} of the jobs in the log of (mode, N)"""
    p = QW / f'{mode}_{N}.jsonl'
    return {(tuple(r['gens']), tuple(r['pat'])): r for r in map(json.loads, p.read_text().splitlines())} if p.exists() else {}


def record(mode, N, rec, hit=False):
    with open(QW / f"{'hits_' if hit else ''}{mode}_{N}.jsonl", 'a') as f:
        f.write(json.dumps(rec) + '\n')


def verify(mode, st, seqs):
    """hit record of four sequences after the exact check K K* = 2N I (raises ValueError otherwise)"""
    rec = {'N': st.N, 'mode': mode, 'gens': list(st.gens), 'subgroup_order': len(st.H), 'rowsums': [int(x.sum()) for x in seqs],
           'sequences': [''.join('+' if v > 0 else '-' for v in x) for x in seqs]}
    if mode == 'bord':
        rec['border'] = [[int(z.real), int(z.imag)] for z in bsz.build(seqs)[2]]
    else:
        ssxx.complex_hadamard(seqs)
    return rec


def run_job(job, st, nbits=0):
    """solve the job, verify a hit (hits_<mode>_<N>.jsonl) and log the result"""
    total, hit, te, tj = solve(job, st, nbits)
    res = dict(job, status='join', total=total, nsol=int(hit is not None), t_enum=round(te, 1), t_join=round(tj, 1))
    if nbits:
        res['slices'] = 2 ** nbits
    if hit is not None:
        record(job['mode'], job['N'], verify(job['mode'], st, hit), hit=True)
        res['HIT'] = True
    record(job['mode'], job['N'], res)
    return res
