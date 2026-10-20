#!/usr/bin/env python3
"""Spence-Paley products: Hadamard matrices of order 4n, n = m q + 1, from Spence's bordered construction over
G = Z_m x GF(q) with q = 3 (mod 4) prime, gcd(m, q) = 1.

Every block is X = u (x) delta + w (x) J + z (x) chi over Z_m x GF(q) (chi the quadratic character, J the all-ones
function): at each position a either w(a) = 0 and u(a), z(a) = +-1 (a z-position: value u at (a, 0), z chi on the
fibre), or z(a) = 0, w(a) = +-1 and u(a) in {0, -2w(a)} (a w-position: the fibre is constant w, the value at (a, 0) is
u + w). With Spence's target sum_k X_k X_k^* = (4mq + 4) delta - 4 J (sum of PAFs -4 off 0), the products J^2 = qJ,
J chi = 0, chi chi^* = q delta - J and chi^* = -chi split the condition over Z[Z_m] into

  stage 1:  sum_k y_k y_k^* = (4mq + 4) delta - 4q J,    y = u + q w   (the compressions to Z_m; row sums 2, 0, 0, 0)
  stage 2:  sum_k z_k z_k^* = R := sum_k (u_k w_k^* + w_k u_k^*) + q sum_k w_k w_k^* + 4 J,
            sum_k (z_k u_k^* - u_k z_k^*) = 0.

Comparing identity coefficients in the J-equation gives (q + 1) N_w - 4 N_- = 4(m - 1) for the numbers of
w-positions and of those with u = -2w. Stage 1 is found by tabu search on the compressions (few equations, many
solutions); stage 2 is, for each stage-1 solution, an exact four-list join over the signs z on the z-positions.
All sequences are constant on the orbits of a multiplier group K <= Z_m^* (times the squares of GF(q)^*).

    python3 search/spp.py n K [--stage1 N] [--seed S] [--log FILE]       (K: generators of K <= Z_m^*, or 1)
"""
import argparse, json, os, struct, subprocess, sys, time
from math import gcd
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import workdir as W


def closure(gens, m):
    S, todo = {1}, [1]
    while todo:
        x = todo.pop()
        for g in gens:
            y = x * g % m
            if y not in S:
                S.add(y); todo.append(y)
    return sorted(S)


def corr(x, y):
    """(x y^*)(b) = sum_a x(a + b) y(a), all b (cyclic, via FFT, exact for small integers)."""
    X, Y = np.fft.fft(x, axis=-1), np.fft.fft(y, axis=-1)
    return np.rint(np.fft.ifft(X * np.conj(Y), axis=-1).real).astype(np.int64)


class SPP:
    def __init__(self, n, kgens, q=None):
        m_q = [(m, q_) for q_ in range(3, n, 4) if all(q_ % d for d in range(2, int(q_ ** .5) + 1))
               for m in [(n - 1) // q_] if (n - 1) % q_ == 0 and gcd(m, q_) == 1]
        if q:
            m_q = [x for x in m_q if x[1] == q]
        self.m, self.q = m_q[0]
        m, q = self.m, self.q
        self.n, self.F = n, (q - 1) // 2
        self.K = closure(kgens, m)
        lab = -np.ones(m, np.int64); self.porb = []
        for a in range(m):
            if lab[a] < 0:
                o = sorted({a * k % m for k in self.K}); lab[o] = len(self.porb); self.porb.append(o)
        self.lab = lab
        self.sizes = np.array([len(o) for o in self.porb])
        L = sorted(set(self.K) | {(-k) % m for k in self.K})
        seen, self.reps, self.anti = set(), [], []
        for b in range(1, m):
            if b not in seen:
                o = {b * k % m for k in L}; seen |= o; self.reps.append(b)
                if (-b) % m not in {b * k % m for k in self.K}:
                    self.anti.append(b)
        self.C = len(self.reps)
        # stage-1 states per position orbit: (f, s) with y = f + 2F s/2 ... s in {-1, 0, 1} (w = s)
        self.states = [(f, w) for w in (0, 1, -1) for f in (1, -1)]   # (f, w): w = 0 z-position
        self.yval = np.array([f + (q - 1) * w for f, w in self.states])

    def describe(self):
        Nw = [(N, (self.q + 1) * N - 4 * (self.m - 1)) for N in range(0, 4 * self.m + 1)]
        ok = [N for N, r in Nw if r >= 0 and r % 4 == 0 and r // 4 <= N]
        return (f"n = {self.n} = {self.m} * {self.q} + 1, F = {self.F}, |K| = {len(self.K)}, position orbits "
                f"{len(self.porb)} (sizes {sorted(set(self.sizes.tolist()))}), shift classes {self.C}, "
                f"antisymmetric classes {len(self.anti)}; admissible N_w {ok[0] if ok else None}..{ok[-1] if ok else None}")

    # ---------------- stage 1: tabu search on the compressions ----------------
    def y_of(self, S):
        """S: 4 x O state indices -> y: 4 x m."""
        return self.yval[S][:, self.lab]

    def energy(self, Y):
        P = corr(Y, Y).sum(0)                       # sum_k PAF_{y_k}
        r = Y.sum(1)
        e = ((P[self.reps] + 4 * self.q) ** 2).sum()
        e += 50 * (int((r ** 2).sum()) - 4) ** 2
        return int(e)

    def stage1(self, rng, maxsteps=200000, tabu=10, restart=4000):
        """Tabu search over the orbit states of all four blocks; every move (block, orbit, state) is scored at once
        by batched FFTs. Energy: squared deviations of sum_k PAF_{y_k} from -4q on the shift classes, plus the row-sum
        condition sum r_k^2 = 4."""
        O, m, q = len(self.porb), self.m, self.q
        onehot = np.zeros((O, m)); onehot[self.lab, np.arange(m)] = 1           # orbit o -> its positions
        reps = np.array(self.reps)
        ytot = 4 * m * q + 4 - 4 * q                                           # sum_k PAF_{y_k}(0), the budget
        S = rng.integers(0, 2, size=(4, O))
        best_e, step0 = None, 0
        tabu_until = np.zeros((4, O, 6), np.int64)
        for step in range(maxsteps):
            Y = self.yval[S][:, self.lab].astype(np.float64)                  # 4 x m
            P = corr(Y, Y).sum(0)
            r = Y.sum(1)
            dv = self.yval[None, None, :] - self.yval[S][:, :, None]          # 4 x O x 6
            D = dv[..., None] * onehot[None, :, None, :]                      # 4 x O x 6 x m
            Yk = np.broadcast_to(Y[:, None, None, :], D.shape)
            dP = corr(D, Yk) + corr(Yk, D) + corr(D, D)
            newP = P[None, None, None, :] + dP
            e = ((newP[..., reps] + 4 * q) ** 2).sum(-1).astype(np.float64)
            dr = dv * self.sizes[None, :, None]
            rsq = (r ** 2).sum() - r[:, None, None] ** 2 + (r[:, None, None] + dr) ** 2
            e = e + 50 * (rsq - 4) ** 2 + (newP[..., 0] - ytot) ** 2
            cur = ((P[reps] + 4 * q) ** 2).sum() + 50 * ((r ** 2).sum() - 4) ** 2 + (P[0] - ytot) ** 2
            if cur == 0:
                return S
            if best_e is None or cur < best_e:
                best_e, step0 = cur, step
            if step - step0 > restart:
                S = rng.integers(0, 2, size=(4, O)); best_e = None; tabu_until[:] = 0; continue
            e[np.arange(4)[:, None], np.arange(O)[None, :], S] = np.inf
            e[tabu_until > step] = np.inf
            flat = np.flatnonzero(e == e.min())
            k, o, st = np.unravel_index(rng.choice(flat), e.shape)
            tabu_until[k, o, S[k, o]] = step + tabu
            S[k, o] = st
        return None

    # ---------------- stage 2: exact join over the z-signs ----------------
    def uwz(self, S):
        f = np.array([[self.states[s][0] for s in row] for row in S])[:, self.lab]
        w = np.array([[self.states[s][1] for s in row] for row in S])[:, self.lab]
        u = np.where(w == 0, f, f - w)
        return u, w

    def stage2_lists(self, S):
        """Per block: all z-sign vectors on its z-orbits and their keys (PAF_z at reps, antisymmetric cross terms)."""
        u, w = self.uwz(S)
        R = (corr(u, w) + corr(w, u)).sum(0) + self.q * corr(w, w).sum(0) + 4
        lists = []
        for k in range(4):
            zorb = [o for o in range(len(self.porb)) if self.states[S[k, o]][1] == 0]
            nz = len(zorb)
            signs = 1 - 2 * ((np.arange(1 << nz)[:, None] >> np.arange(nz)[None, :]) & 1)
            Zo = np.zeros((1 << nz, len(self.porb)), np.int64); Zo[:, zorb] = signs
            Z = Zo[:, self.lab]
            pz = corr(Z, Z)[:, self.reps]
            cross = (corr(Z, np.broadcast_to(u[k], Z.shape)) - corr(np.broadcast_to(u[k], Z.shape), Z))[:, self.anti]
            lists.append((zorb, signs, np.concatenate([pz, cross], 1)))
        target = np.concatenate([R[self.reps], np.zeros(len(self.anti), np.int64)])
        return lists, target

    def sequences(self, S, zs):
        """Full sequences on Z_m x GF(q), indexed (a, y) -> a * q + y (mixed radix of [m, q])."""
        q = self.q
        chi = -np.ones(q, np.int64); chi[np.arange(1, q) ** 2 % q] = 1; chi[0] = 0
        u, w = self.uwz(S)
        out = []
        for k in range(4):
            z = np.zeros(len(self.porb), np.int64); z[zs[k][0]] = zs[k][1]
            zf = z[self.lab]
            X = u[k][:, None] * (np.arange(q) == 0)[None, :] + w[k][:, None] + zf[:, None] * chi[None, :]
            out.append(X.reshape(-1))
        return out


def join_lists(V, target, d):
    """Exact four-list join (search/join4x.c); returns (i, j, k, l) or None."""
    path = W.WORK / f"spp_{os.getpid()}.bin"
    C = V[0].shape[1]
    with open(path, "wb") as f:
        f.write(struct.pack("5i", *[len(x) for x in V], C))
        f.write(np.asarray(target, np.int16).tobytes())
        for x in V:
            f.write(np.ascontiguousarray(x, np.int16).tobytes())
    out = subprocess.run([W.bin("join4x"), str(path), "1"], capture_output=True, text=True, check=True).stdout
    path.unlink()
    for ln in out.splitlines():
        if ln.startswith("SOL"):
            return tuple(map(int, ln.split()[1:5]))
    return None


def spence_check(seqs, m, q):
    """sum of PAFs over Z_m x GF(q) = -4 off 0 (via the isomorphism with Z_mq)."""
    v = m * q
    idx = np.array([(t % m) * q + (t % q) for t in range(v)])     # CRT: Z_mq -> Z_m x Z_q
    P = sum(corr(x[idx], x[idx]) for x in seqs)
    return P[0] == 4 * v and np.all(P[1:] == -4), [int(x.sum()) for x in seqs]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("n", type=int); ap.add_argument("K", type=int, nargs="+")
    ap.add_argument("--q", type=int); ap.add_argument("--stage1", type=int, default=50)
    ap.add_argument("--seed", type=int, default=1); ap.add_argument("--log")
    ap.add_argument("--maxlist", type=float, default=3e5)
    a = ap.parse_args()
    sp = SPP(a.n, a.K, a.q)
    print(sp.describe(), flush=True)
    rng = np.random.default_rng(a.seed)
    found = 0
    for it in range(a.stage1):
        t0 = time.time()
        S = sp.stage1(rng)
        if S is None:
            print(f"  stage 1 #{it}: not found", flush=True); continue
        Nw = int(sum(sp.sizes[o] for k in range(4) for o in range(len(sp.porb)) if sp.states[S[k, o]][1] != 0))
        lists, target = sp.stage2_lists(S)
        sizes = [len(x[1]) for x in lists]
        msg = f"  stage 1 #{it}: N_w = {Nw}, z-lists {sizes} [{time.time() - t0:.1f}s]"
        if max(sizes) > a.maxlist:
            print(msg + " (lists too large, skipped)", flush=True); continue
        hit = join_lists([x[2] for x in lists], target, None)
        print(msg + f", stage 2: {'SOLUTION' if hit else 'none'} [{time.time() - t0:.1f}s]", flush=True)
        if hit:
            zs = [(lists[k][0], lists[k][1][hit[k]]) for k in range(4)]
            seqs = sp.sequences(S, zs)
            ok, rows = spence_check(seqs, sp.m, sp.q)
            print(f"  Spence condition over Z_{sp.m} x GF({sp.q}): {ok}, row sums {rows}", flush=True)
            if ok and a.log:
                with open(a.log, "a") as f:
                    f.write(json.dumps({"n": a.n, "m": sp.m, "q": sp.q, "K": sp.K, "rows": rows,
                                        "sequences": ["".join("+" if v > 0 else "-" for v in x) for x in seqs]}) + "\n")
            found += ok
            if found:
                break


if __name__ == "__main__":
    main()
