"""T-matrices of prime order t by a support split and the four-list k-tree matcher (match4.c).

T-matrices of order t are four circulant (0, +-1) matrices T1..T4 whose supports partition Z_t and whose
periodic autocorrelations sum to t at shift 0 and to 0 elsewhere. With Williamson matrices of order w they
give a Hadamard matrix of order 4tw (Cooper-Wallis).

The split. Let E = supp T1 u supp T2 and F = Z_t \\ E = supp T3 u supp T4, and put

    X1 = T1 + T2,  X3 = T1 - T2   (+-1 on E, 0 on F),
    X2 = T3 + T4,  X4 = T3 - T4   (+-1 on F, 0 on E).

Then PAF(X1) + PAF(X3) = 2 (PAF(T1) + PAF(T2)) and likewise on F, so the T_i are T-matrices if and only if

    PAF(X1) + PAF(X2) + PAF(X3) + PAF(X4) = 2t [k = 0].

Conversely, any four such sequences give T-matrices: T1 = (X1 + X3)/2 and T2 = (X1 - X3)/2 have disjoint
supports covering E, and the same holds on F. So once E is fixed, the four sequences are independent and
the problem is a four-list problem of Goethals-Seidel type on two fixed supports. (Base sequences are the
case where E is an interval.) All four sequences are taken constant on the orbits of a multiplier group H
of odd order, so E is a union of H-orbits and the autocorrelations are constant on the classes
H-orbits of {k, -k}.

Each list holds every +-1 assignment of the orbits on its support with a given row sum; for x_i the row
sum of X_i, x1^2 + x2^2 + x3^2 + x4^2 = 2t. Row-sum patterns are ranked by the Gaussian estimate of the
number of solutions in the four lists (t is prime here, where that estimate has been reliable).

  python3 search/tmatrix.py T GEN [--trials K] [--seed S] [--maxtau M] [--minlog L] [--cap C]

Hits go to $HADAMARD_WORK/tmatrix_hits.jsonl and T<t>.txt files there (0, +, - first rows, one per line).
"""
import argparse, json, math, os, subprocess, sys, tempfile, time
from pathlib import Path as _P
sys.path[:0] = [str(_P(__file__).resolve().parent)]
import numpy as np
import workdir as W


def group(t, gens):
    H = {1}
    while True:
        new = {h * g % t for h in H for g in gens} | H
        if new == H: return sorted(H)
        H = new


class Structure:
    def __init__(self, t, gens):
        self.t, self.H = t, group(t, gens)
        assert len(self.H) % 2 == 1, "H must have odd order (-1 not in H)"
        orb = -np.ones(t, dtype=np.int64); reps = []
        for x in range(t):
            if orb[x] < 0:
                for h in self.H: orb[x * h % t] = len(reps)
                reps.append(x)
        self.orb, self.O = orb, len(reps)
        self.size = np.bincount(orb, minlength=self.O)
        cls = -np.ones(t, dtype=np.int64); self.kc = []
        for k in range(1, t):
            if cls[k] < 0:
                for h in self.H: cls[k * h % t] = cls[(t - k * h) % t] = len(self.kc)
                self.kc.append(k)
        self.C = len(self.kc)
        self.csize = np.bincount(cls[1:], minlength=self.C)
        # N[c, o, p] = #{x in orbit o : x + k_c in orbit p}
        self.N = np.zeros((self.C, self.O, self.O), dtype=np.int64)
        x = np.arange(t)
        for c, k in enumerate(self.kc):
            np.add.at(self.N[c], (orb[x], orb[(x + k) % t]), 1)

    def side(self, P, rowsum, chunk=1 << 20):
        """All +-1 assignments y of the orbits P with sum_o y_o |o| = rowsum, and the autocorrelation classes
        of the sequence X = y on P, 0 elsewhere, encoded v = (PAF - (m mod 2)) / 2."""
        P = np.asarray(P); n = len(P); sz = self.size[P]
        Q = self.N[:, P][:, :, P].astype(np.float32)               # C x n x n
        m = self.N[:, P][:, :, P].sum((1, 2))                        # pairs inside the support
        Ys, Vs = [], []
        for lo in range(0, 1 << n, chunk):
            idx = np.arange(lo, min(lo + chunk, 1 << n), dtype=np.int64)
            Y = 1 - 2 * ((idx[:, None] >> np.arange(n)) & 1)
            keep = (Y * sz).sum(1) == rowsum
            Y = Y[keep]
            if not len(Y): continue
            Yf = Y.astype(np.float32)
            paf = np.einsum('ri,cij,rj->rc', Yf, Q, Yf, optimize=True)
            paf = np.rint(paf).astype(np.int64)
            Vs.append((paf - (m % 2)) // 2); Ys.append(Y.astype(np.int8))
        if not Ys: return np.zeros((0, n), np.int8), np.zeros((0, self.C), np.int64), m
        return np.concatenate(Ys), np.concatenate(Vs), m

    def seq(self, P, y):
        X = np.zeros(self.t, dtype=np.int64)
        for o, s in zip(P, y): X[self.orb == o] = s
        return X


def paf(X):
    F = np.fft.fft(X)
    return np.rint(np.fft.ifft(F * np.conj(F)).real).astype(np.int64)


def is_tmatrix(T):
    T = np.asarray(T); t = T.shape[1]
    if not np.isin(T, (-1, 0, 1)).all() or not ((T != 0).sum(0) == 1).all(): return False
    s = sum(paf(x) for x in T)
    return s[0] == t and not s[1:].any()


def write(path, V):
    rec = np.zeros((len(V), 8 + V.shape[1]), np.uint8)
    rec[:, :8] = np.arange(len(V), dtype=np.uint64).view(np.uint8).reshape(-1, 8)
    rec[:, 8:] = V.astype(np.int8).view(np.uint8)
    rec.tofile(path)


def logE(Vs):
    """Gaussian estimate of log10 of the number of quadruples with sum 0 on every class."""
    lg = sum(math.log10(len(V)) for V in Vs)
    for c in range(Vs[0].shape[1]):
        mu = sum(V[:, c].mean() for V in Vs); var = sum(V[:, c].var() for V in Vs) + 1e-9
        lg += (-mu * mu / (2 * var)) / math.log(10) - 0.5 * math.log10(2 * math.pi * var)
    return lg


def clampS(Vl, Vr, cap):
    S = 0
    for (X, Y) in (Vl, Vr):
        lp, s = math.log(len(X)) + math.log(len(Y)), 0
        while lp > math.log(cap) and s < X.shape[1]:
            hx = np.bincount(X[:, s] + 128, minlength=256) / len(X)
            hy = np.bincount(Y[:, s] + 128, minlength=256) / len(Y)
            lp += math.log(np.convolve(hx, hy).max()); s += 1
        S = max(S, s)
    return S


def patterns(st, PE, PF):
    """Row sums (x1, x3) on E and (x2, x4) on F with x1 >= x3 >= 0, x2 >= x4 >= 0 and sum of squares 2t.
    Negating one X_i or swapping X1, X3 (or X2, X4) maps solutions to solutions."""
    def sums(P):
        s = {0}
        for o in P: s = {a + b for a in s for b in (st.size[o], -st.size[o])}
        return s
    sE, sF = sums(PE), sums(PF); out = []
    for x1 in sorted(x for x in sE if x >= 0):
        for x3 in sorted(x for x in sE if 0 <= x <= x1):
            for x2 in sorted(x for x in sF if x >= 0):
                r = 2 * st.t - x1 * x1 - x3 * x3 - x2 * x2
                if r < 0: continue
                x4 = math.isqrt(r)
                if x4 * x4 == r and x4 <= x2 and x4 in sF: out.append((x1, x2, x3, x4))
    return out


def search_E(st, PE, PF, args, rng, log):
    t = st.t; found = []
    pats = patterns(st, PE, PF)
    cache = {}
    def side(P, r, key):
        if (key, r) not in cache: cache[(key, r)] = st.side(P, r)
        return cache[(key, r)]
    ranked = []
    for (x1, x2, x3, x4) in pats:
        (Y1, V1, mE), (Y2, V2, mF) = side(PE, x1, 'E'), side(PF, x2, 'F')
        (Y3, V3, _), (Y4, V4, _) = side(PE, x3, 'E'), side(PF, x4, 'F')
        if min(len(Y1), len(Y2), len(Y3), len(Y4)) == 0: continue
        fold = (mE % 2) + (mF % 2)
        Vs = [V1 + fold, V2, V3, V4]
        ranked.append((logE(Vs), (x1, x2, x3, x4), [Y1, Y2, Y3, Y4], Vs))
    ranked.sort(key=lambda r: -r[0])
    log(f"  E: {len(PE)} orbits ({int(st.size[PE].sum())} points), F: {len(PF)} orbits; patterns {len(ranked)}; "
        f"best log10 E {[(round(r[0], 2), r[1]) for r in ranked[:4]]}")
    for lg, pat, Ys, Vs in ranked[:args.patterns]:
        if lg < args.minlog: break
        if max(np.abs(V).max() for V in Vs) > 127:
            keep = [np.abs(V).max(1) <= 127 for V in Vs]
            Ys = [Y[k] for Y, k in zip(Ys, keep)]; Vs = [V[k] for V, k in zip(Vs, keep)]
        t0 = time.time()
        n = [len(V) for V in Vs]
        if max(n[0] * n[1], n[2] * n[3], n[0] * n[3], n[1] * n[2]) <= args.exact:
            # exact four-list join (find4b pairs slots (0, 1) and (2, 3)); every solution is found
            order = [0, 1, 2, 3] if n[0] * n[1] <= n[0] * n[3] else [0, 3, 1, 2]
            fd, path = tempfile.mkstemp(suffix='.bin', dir=W.tmpdir()); os.close(fd)
            with open(path, 'wb') as f:
                f.write(np.array([n[k] for k in order] + [st.C], dtype=np.int32).tobytes())
                f.write(np.zeros(st.C, dtype=np.int16).tobytes())
                for k in order: f.write(np.ascontiguousarray(Vs[k].astype(np.int16)).tobytes())
            res = subprocess.run([W.bin('find4b'), path, '4'], capture_output=True, text=True)
            os.remove(path)
            sol_lines = [list(map(int, l.split()[1:])) for l in res.stdout.splitlines() if l.startswith('SOL')]
            tail = res.stdout.strip().splitlines()[-1] if res.stdout.strip() else res.stderr.strip()[-200:]
            log(f"    pattern {pat} sizes {n} log10E {lg:.2f} exact: {tail} [{time.time()-t0:.0f}s]")
        else:
            # match4 pairs (A, D) and (B, C): A = X1, D = X2 (left); B = X3, C = X4 (right)
            order = [0, 2, 3, 1]
            S = clampS((Vs[0], Vs[1]), (Vs[2], Vs[3]), args.cap); S = max(1, min(S, st.C - 2))
            tmp = tempfile.mkdtemp(dir=W.tmpdir())
            paths = [os.path.join(tmp, f"{k}.bin") for k in range(4)]
            for p, k in zip(paths, order): write(p, Vs[k])
            cmd = [W.bin('match4'), '-s', str(S), '-t', '0', '-n', str(args.maxtau), '-k', '4', '-m', str(args.top)] + paths + [str(st.C)]
            res = subprocess.run(cmd, capture_output=True, text=True)
            for p in paths: os.remove(p)
            os.rmdir(tmp)
            err = res.stderr.strip().splitlines()
            log(f"    pattern {pat} sizes {n} log10E {lg:.2f} S={S}: {err[-1] if err else ''} [{time.time()-t0:.0f}s]")
            sol_lines = [[int(x, 16) for x in l.split()[1:]] for l in res.stdout.splitlines() if l.startswith('SOL')]
        for ids in sol_lines:
            k = {order[j]: ids[j] for j in range(4)}
            X = [st.seq(PE if i in (0, 2) else PF, Ys[i][k[i]]) for i in range(4)]
            T = np.array([(X[0] + X[2]) // 2, (X[0] - X[2]) // 2, (X[1] + X[3]) // 2, (X[1] - X[3]) // 2])
            if is_tmatrix(T):
                found.append(T); log(f"    T-matrices of order {t} found, pattern {pat}")
                return found
    return found


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('t', type=int); ap.add_argument('gens', type=int, nargs='+')
    ap.add_argument('--trials', type=int, default=20); ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--maxtau', type=int, default=400); ap.add_argument('--top', type=int, default=6)
    ap.add_argument('--patterns', type=int, default=3); ap.add_argument('--minlog', type=float, default=-3.0)
    ap.add_argument('--cap', type=float, default=2e7)
    ap.add_argument('--exact', type=float, default=3e8, help='use the exact join find4b when every pair product is below this')
    ap.add_argument('--half', type=int, default=None, help='number of nonzero orbits in E (default O//2)')
    args = ap.parse_args()
    st = Structure(args.t, args.gens)
    rng = np.random.default_rng(args.seed)
    log = lambda s: print(s, flush=True)
    log(f"t={args.t} H={args.gens} |H|={len(st.H)} orbits={st.O} classes={st.C} seed={args.seed}")
    nz = np.arange(1, st.O); half = args.half or (st.O - 1) // 2
    for trial in range(args.trials):
        sel = rng.permutation(nz)
        PE, PF = list(sel[:half]), list(sel[half:])
        (PE if rng.random() < 0.5 else PF).insert(0, 0)
        PE, PF = sorted(PE), sorted(PF)
        log(f"trial {trial}: E = orbits {PE}")
        T = search_E(st, PE, PF, args, rng, log)
        if T:
            T = T[0]; t = args.t
            rows = [''.join({0: '0', 1: '+', -1: '-'}[int(v)] for v in r) for r in T]
            out = W.WORK / f"T{t}.txt"; out.write_text('\n'.join(rows) + '\n')
            with open(W.hits('tmatrix_hits.jsonl'), 'a') as f:
                f.write(json.dumps(dict(t=t, gens=args.gens, seed=args.seed, trial=trial, E=[int(o) for o in PE],
                                        rows=rows)) + '\n')
            log(f"wrote {out}")
            return
    log("no T-matrices found")


if __name__ == '__main__':
    main()
