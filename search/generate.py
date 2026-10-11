#!/usr/bin/env python3
"""Generate candidate sequences for the Goethals-Seidel search (JAX; TPU, GPU or CPU).

A candidate is a +-1 sequence on Z_p that is constant on the cosets of H = <g> and
takes the value +1 at 0, with a prescribed row sum. Starting from random sequences,
batched local search minimises

    F(x) = sum_c (P_x(c) - m_c)^2,

where P_x(c) is the periodic autocorrelation of x on class c (the classes are the
orbits of <H, -1> on Z_p \\ {0}) and m_c is the exact mean of P(c) over all
sequences with the same row sum. Moves swap the signs of one positive and one
negative coset, which keeps the row sum fixed. Three phases are run on every batch:
greedy pair moves, best-of-m pair moves, and iterated local search (random swaps
followed by re-descent, accepted when F decreases).

Output: binary records, one per candidate,
    uint64  mask   bit (o - 1) set when coset o (numbered from 1) is negative
    int8[C] v      v_c = (P_x(c) - (p mod 4)) / 4
which is the input format of dedupe and match.

    python3 search/generate.py 1123 60 --rowsum 67 --count 160000000 --out r67.bin
"""
import argparse
import sys
import time

import numpy as np
import jax
import jax.numpy as jnp


def orbit_structure(p, g):
    """Coset numbering (orbit 0 = {0}), autocorrelation classes, and the pair-count
    tensor S[c, o, q] = #{t : t in O_o, t + s_c in O_q} + (o <-> q)."""
    orb = -np.ones(p, dtype=np.int64)
    n_orb = 0
    for j in range(p):
        if orb[j] >= 0:
            continue
        x = j
        while orb[x] < 0:
            orb[x] = n_orb
            x = x * g % p
        n_orb += 1
    cls = -np.ones(p, dtype=np.int64)
    reps = []
    for s in range(1, p):
        if cls[s] >= 0:
            continue
        x = s
        while cls[x] < 0:
            cls[x] = cls[p - x] = len(reps)
            x = x * g % p
        reps.append(s)
    S = np.zeros((len(reps), n_orb, n_orb))
    t = np.arange(p)
    for c, s in enumerate(reps):
        np.add.at(S[c], (orb, orb[(t + s) % p]), 1)
    return orb, reps, S + S.transpose(0, 2, 1)


def class_means(S, negatives):
    """Exact mean of P(c) over sequences with `negatives` negative cosets (x(0) = +1)."""
    n1 = S.shape[1] - 1
    a = n1 - 2.0 * negatives
    E = np.full((n1 + 1, n1 + 1), (a * a - n1) / (n1 * (n1 - 1)))
    E[0, :] = E[:, 0] = a / n1
    np.fill_diagonal(E, 1.0)
    return 0.5 * np.einsum("coq,oq->c", S, E)


def make_batch_fn(S, means, negatives, batch, greedy, polish, width, rounds, swaps, steps):
    C, N = S.shape[0], S.shape[1]
    f32 = jnp.float32
    S_ = jnp.array(S, f32)
    S_rows = jnp.array(S.transpose(1, 0, 2).reshape(N, C * N), f32)
    diag = jnp.array(np.stack([np.diag(S[c]) for c in range(C)]), f32)
    Q = jnp.array((S * S).sum(0), f32)
    m = jnp.array(means, f32)
    free = jnp.array(np.arange(N) > 0)

    def full(y):
        G = jnp.einsum("bq,coq->bco", y, S_, precision="highest")
        return G, 0.5 * jnp.einsum("bo,bco->bc", y, G, precision="highest")

    def flip(y, G, T, h, move):
        yo = jnp.sum(y * h, axis=1)
        Go = jnp.einsum("bco,bo->bc", G, h, precision="highest")
        So = (h @ S_rows).reshape(batch, C, N)
        T = jnp.where(move[:, None], T - 2.0 * yo[:, None] * Go + 2.0 * (h @ diag.T), T)
        G = jnp.where(move[:, None, None], G - 2.0 * yo[:, None, None] * So, G)
        y = jnp.where(move[:, None] & (h > 0), -y, y)
        return y, G, T

    def step(state, m_best):
        """Apply the best swap among the m_best most promising positive cosets."""
        y, G, T = state
        d = -2.0 * y[:, None, :] * G + 2.0 * diag[None]       # change of P(c) when coset o flips
        e = T - m[None]
        f = jnp.sum(d * (d + 2.0 * e[:, :, None]), axis=1)    # change of F for a single flip
        _, cand = jax.lax.top_k(-jnp.where((y > 0) & free[None], f, jnp.inf), m_best)
        best = jnp.full((batch,), jnp.inf)
        o1_best = jnp.zeros((batch,), jnp.int32)
        o2_best = jnp.zeros((batch,), jnp.int32)
        for j in range(m_best):
            o1 = cand[:, j]
            h1 = jax.nn.one_hot(o1, N, dtype=f32)
            d1 = jnp.einsum("bco,bo->bc", d, h1, precision="highest")
            S1 = (h1 @ S_rows).reshape(batch, C, N)
            dF = (jnp.sum(f * h1, axis=1, keepdims=True) + f + 2.0 * jnp.sum(d1[:, :, None] * d, axis=1)
                  + 16.0 * (h1 @ Q) - 8.0 * jnp.sum(S1 * (d1[:, :, None] + d + e[:, :, None]), axis=1))
            dF = jnp.where(y < 0, dF, jnp.inf)
            o2 = jnp.argmin(dF, axis=1)
            v = jnp.min(dF, axis=1)
            better = v < best
            best = jnp.where(better, v, best)
            o1_best = jnp.where(better, o1, o1_best)
            o2_best = jnp.where(better, o2, o2_best)
        move = best < -1e-3
        for o in (o1_best, o2_best):
            y, G, T = flip(y, G, T, jax.nn.one_hot(o, N, dtype=f32), move)
        return y, G, T

    def objective(T):
        return jnp.sum((T - m[None]) ** 2, axis=1)

    @jax.jit
    def run(key):
        k_init, k_ils = jax.random.split(key)
        ranks = jnp.argsort(jnp.argsort(jax.random.uniform(k_init, (batch, N - 1)), axis=1), axis=1)
        y = jnp.concatenate([jnp.ones((batch, 1)), jnp.where(ranks < negatives, -1.0, 1.0)], axis=1)
        state = (y,) + full(y)
        state = jax.lax.fori_loop(0, greedy, lambda _, s: step(s, 1), state)
        state = jax.lax.fori_loop(0, polish, lambda _, s: step(s, width), state)

        def ils(_, carry):
            (y, G, T), k = carry
            k, kp = jax.random.split(k)
            y2 = y
            for _s in range(swaps):
                kp, ka, kb = jax.random.split(kp, 3)
                up = jnp.where((y2 > 0) & free[None], jax.random.uniform(ka, (batch, N)), -1.0)
                dn = jnp.where(y2 < 0, jax.random.uniform(kb, (batch, N)), -1.0)
                fl = jax.nn.one_hot(jnp.argmax(up, 1), N) + jax.nn.one_hot(jnp.argmax(dn, 1), N)
                y2 = jnp.where(fl > 0, -y2, y2)
            trial = jax.lax.fori_loop(0, steps, lambda _, s: step(s, width), (y2,) + full(y2))
            accept = objective(trial[2]) < objective(T)
            y = jnp.where(accept[:, None], trial[0], y)
            G = jnp.where(accept[:, None, None], trial[1], G)
            T = jnp.where(accept[:, None], trial[2], T)
            return (y, G, T), k

        state, _ = jax.lax.fori_loop(0, rounds, ils, (state, k_ils))
        return state[0], state[2]

    return run


def main():
    ap = argparse.ArgumentParser(description="Generate Goethals-Seidel candidate sequences.")
    ap.add_argument("p", type=int, help="prime modulus")
    ap.add_argument("g", type=int, help="generator of the multiplier subgroup H")
    ap.add_argument("--rowsum", type=int, required=True, help="row sum of the generated sequences")
    ap.add_argument("--count", type=int, required=True, help="number of candidates to write")
    ap.add_argument("--out", required=True, help="output file")
    ap.add_argument("--batch", type=int, default=32768)
    ap.add_argument("--greedy", type=int, default=6, help="greedy pair moves")
    ap.add_argument("--polish", type=int, default=4, help="best-of-m pair moves")
    ap.add_argument("--width", type=int, default=6, help="m in the best-of-m moves")
    ap.add_argument("--rounds", type=int, default=6, help="iterated local search rounds")
    ap.add_argument("--swaps", type=int, default=3, help="random swaps per perturbation")
    ap.add_argument("--steps", type=int, default=6, help="re-descent moves per round")
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()

    orb, reps, S = orbit_structure(a.p, a.g)
    N, h = S.shape[1], (a.p - 1) // (S.shape[1] - 1)
    # x(0) = +1 forces r = 1 (mod 2h); -x has the same autocorrelation, so either sign of r will do
    for r in (a.rowsum, -a.rowsum):
        num = (N - 1) - (r - 1) // h
        if (r - 1) % h == 0 and num % 2 == 0 and 0 <= num // 2 <= N - 1:
            break
    else:
        sys.exit(f"row sum {a.rowsum} is not attainable (need |r| = +-1 mod {2 * h})")
    negatives = num // 2
    means = class_means(S, negatives)
    run = make_batch_fn(S, means, negatives, a.batch, a.greedy, a.polish, a.width,
                        a.rounds, a.swaps, a.steps)
    print(f"p={a.p} |H|={h} cosets={N - 1} classes={len(reps)} rowsum={r} "
          f"negative cosets={negatives} device={jax.devices()[0]}", file=sys.stderr, flush=True)

    weights = np.array([1 << i for i in range(N - 1)], dtype=np.uint64)
    key = jax.random.PRNGKey(a.seed)
    done, t0 = 0, time.time()
    with open(a.out, "wb") as out:
        while done < a.count:
            key, sub = jax.random.split(key)
            y, T = run(sub)
            y = np.asarray(y)
            P = np.rint(np.asarray(T)).astype(np.int64)
            rec = np.zeros((a.batch, 8 + len(reps)), np.uint8)
            rec[:, :8] = ((y[:, 1:] < 0).astype(np.uint64) * weights).sum(1).astype(np.uint64) \
                .view(np.uint8).reshape(a.batch, 8)
            rec[:, 8:] = ((P - a.p % 4) // 4).astype(np.int8).view(np.uint8)
            out.write(rec.tobytes())
            done += a.batch
            if done % (10 * a.batch) == 0 or done >= a.count:
                rate = done / (time.time() - t0)
                var = np.mean((P - means[None]) ** 2)
                print(f"{done} candidates, {rate:.3g}/s, mean class variance {var:.1f}",
                      file=sys.stderr, flush=True)


if __name__ == "__main__":
    main()
