"""Low-energy cuts for the B and D lists of the three-list propus search (heuristic; the exact search needs none).

Model.  A solution has F_A + 2 F_B + F_D = T (class values, T = 0).  Treat the class-value vectors of the three
(PSD-filtered) pools as Gaussian with means mu_X and covariances S_X.  Given B, the vector -F_A - F_D must equal 2 F_B:
its density is N(T - mu_A - mu_D, S_A + S_D) at 2 F_B, so
        E_B(F) = | W_B (F - m_B) |^2,   m_B = (T - mu_A - mu_D)/2,   W_B = 2 (S_A + S_D)^(+1/2);
given D, F_A + 2 F_B must equal T - F_D:
        E_D(F) = | W_D (F - m_D) |^2,   m_D = T - mu_A - 2 mu_B,   W_D = (S_A + 4 S_B)^(+1/2)
(pseudo-inverse: the row sums make the covariances singular).  A list keeps the members with E <= tau.  The moments are
those of the PSD-filtered pools and are averaged over the multiplier group (class permutations), so that E is invariant under
Mult and under the reversals: a canonical representative passes iff its whole orbit does.  The cut is a filter, not a
theorem: it can discard solutions.  `study` measures how many solution orbits survive a cut that keeps a fraction f of
each pool (exact pools on small structures), `checkc` verifies that the C engine's cut equals the numpy cut exactly.

  python3.11 p3lowe.py study      retention vs list fraction on small structures with many solutions
  python3.11 p3lowe.py checkc     C engine energy filter vs numpy, exact set equality
  python3.11 p3lowe.py plan V GEN --rows a,b,d [--samples N] [--fB f] [--fD f] --out PROBLEM   (sampled moments)
"""
import argparse, sys, tempfile, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from p3struct import Structure
import p3run as P
import p3test as T


def class_perm(st, t):
    return np.array([st.cls[(t * s) % st.v] for s in st.creps])


def symmetrize(st, mean, cov):
    """Average (mean, covariance) over the class permutations of Mult."""
    second = cov + np.outer(mean, mean)
    perms = [class_perm(st, t) for t in st.mult_reps]
    m = np.mean([mean[p] for p in perms], axis=0)
    M = np.mean([second[np.ix_(p, p)] for p in perms], axis=0)
    return m, M - np.outer(m, m)


def pinv_sqrt(M, rel=1e-9):
    w, U = np.linalg.eigh((M + M.T) / 2)
    keep = w > rel * w.max()
    return (U[:, keep] / np.sqrt(w[keep])) @ U[:, keep].T


def energies(F, m, W):
    d = (np.asarray(F, dtype=np.float64) - m) @ W.T
    return (d * d).sum(1)


def psd_pool_exact(st, blk, r, bound):
    got = T.all_vectors(blk, r, 1 << 22)
    masks, Y = got
    F = blk.values(Y)
    ok = (st.psd(F) <= bound + 0.5).all(1)          # same slack as the engine (psd_ok: bound + 0.5): the filter is only necessary
    return masks[ok], Y[ok], F[ok]


def moments_of(F):
    return F.mean(0), np.cov(F.T, bias=True)


def energy_model(st, momA, momB, momD, tgt=None):
    tgt = np.zeros(st.C) if tgt is None else np.asarray(tgt, dtype=np.float64)
    (mA, SA), (mB, SB), (mD, SD) = [symmetrize(st, *m) for m in (momA, momB, momD)]
    WB = 2 * pinv_sqrt(SA + SD)
    WD = pinv_sqrt(SA + 4 * SB)
    return dict(B=((tgt - mA - mD) / 2, WB), D=(tgt - mA - 2 * mB, WD))


def midpoint_tau(E, f):
    """tau keeping a fraction ~f of the members, placed between two distinct energies (no float ties)."""
    s = np.sort(E)
    k = max(1, min(len(s) - 1, int(round(f * len(s)))))
    while k < len(s) - 1 and s[k] - s[k - 1] < 1e-9 * max(1.0, s[k]):
        k += 1
    return 0.5 * (s[k - 1] + s[k])


STUDY = [(49, 18, (11, 5, 5)), (49, 18, (13, 1, 5)), (43, 6, (7, 1, 11)), (43, 6, (5, 7, 7)), (43, 6, (7, 7, 5)),
         (31, 5, (5, 5, 7)), (37, 10, (1, 7, 7)), (37, 10, (11, 1, 5)), (39, 16, (3, 7, 7)), (71, 5, (11, 9, 1))]


def study(cases=STUDY, fracs=(0.5, 0.3, 0.2, 0.1, 0.05, 0.02)):
    tot = {f: [0, 0] for f in fracs}
    for (v, u, rows) in cases:
        st = Structure(v, [u])
        if max(st.B.O, st.A.O) > 17 and v != 71:
            continue
        ref = T.reference_solutions(st, rows)
        if not ref:
            continue
        canon = sorted(s for s in ref if P.is_canonical(st, *s))
        mA, YA, FA = psd_pool_exact(st, st.A, rows[0], 4 * v)
        mB, YB, FB = psd_pool_exact(st, st.B, rows[1], 2 * v)
        mD, YD, FD = psd_pool_exact(st, st.B, rows[2], 4 * v)
        em = energy_model(st, moments_of(FA), moments_of(FB), moments_of(FD))
        EB, ED = energies(FB, *em['B']), energies(FD, *em['D'])
        # solution members
        def vals(blk, ms):
            Y = np.array([[1 if (m >> i) & 1 else -1 for i in range(blk.O)] for m in ms], dtype=np.int8)
            return blk.values(Y)
        sFB = vals(st.B, [s[1] for s in canon]); sFD = vals(st.B, [s[2] for s in canon])
        sEB, sED = energies(sFB, *em['B']), energies(sFD, *em['D'])
        print(f'v={v} |H|={len(st.H)} rows={rows}: PSD pools A {len(mA)} B {len(mB)} D {len(mD)}; {len(canon)} solution orbits '
              f'(canonical reps); |Gamma|={P.gamma_order(st)}', flush=True)
        line = '   f_B=f_D :  ' + '  '.join(f'{f:5.2f}' for f in fracs)
        ret = []
        for f in fracs:
            tB, tD = midpoint_tau(EB, f), midpoint_tau(ED, f)
            r = float(np.mean((sEB <= tB) & (sED <= tD)))
            ret.append(r)
            tot[f][0] += int(round(r * len(canon))); tot[f][1] += len(canon)
        print(line); print('   retention  :  ' + '  '.join(f'{r:5.2f}' for r in ret))
        print('   gain r/f^2 :  ' + '  '.join(f'{r / (f * f):5.1f}' for r, f in zip(ret, fracs)), flush=True)
    print('pooled over all solution orbits above:')
    for f in fracs:
        k, n = tot[f]
        print(f'   keep {f:4.2f} of each list: {k}/{n} = {k / max(n, 1):.2f} of the solution orbits retained; a random cut would keep {f * f:.4f}; gain {k / max(n, 1) / (f * f):.1f}x; join work x{f * f:.4f} (pairs)')
    return tot


def checkc(tmpdir):
    """The C engine's energy cut must reproduce the numpy cut exactly (B and D lists and the solution subset)."""
    ok = True
    for (v, u, rows, f) in ((49, 18, (11, 5, 5), 0.3), (43, 6, (7, 1, 11), 0.2), (31, 5, (5, 5, 7), 0.5)):
        st = Structure(v, [u])
        ref = T.reference_solutions(st, rows)
        mA, YA, FA = psd_pool_exact(st, st.A, rows[0], 4 * v)
        mB, YB, FB = psd_pool_exact(st, st.B, rows[1], 2 * v)
        mD, YD, FD = psd_pool_exact(st, st.B, rows[2], 4 * v)
        em = energy_model(st, moments_of(FA), moments_of(FB), moments_of(FD))
        EB, ED = energies(FB, *em['B']), energies(FD, *em['D'])
        tB, tD = midpoint_tau(EB, f), midpoint_tau(ED, f)
        passB = {int(m) for m, e in zip(mB, EB) if e <= tB}
        passD = {int(m) for m, e in zip(mD, ED) if e <= tD}
        want = {s for s in ref if s[1] in passB and s[2] in passD}
        pf = Path(tmpdir) / 'low.txt'
        P.write_problem(st, rows, pf, canonB=False, canonD=False, psd=(True, True, True), seed=3,
                        energy=dict(B=(tB, *em['B']), D=(tD, *em['D'])))
        r = P.run_p3(pf, threads=2, ns=2)
        got = set(r['hits'])
        sizes_ok = r['lists'][1] == len(passB) and r['lists'][2] == len(passD) and r['lists'][0] == len(mA)
        print(f'v={v} rows={rows} f={f}: lists C {r["lists"]} vs numpy ({len(mA)}, {len(passB)}, {len(passD)}) {"ok" if sizes_ok else "MISMATCH"}; '
              f'hits C {len(got)} vs numpy-subset {len(want)} of {len(ref)} {"ok" if got == want else "MISMATCH"}', flush=True)
        ok &= sizes_ok and got == want
    print('ALL OK' if ok else 'FAILED')
    return ok


def plan(v, gen, rows, samples, fB, fD, out, seed=1):
    """Sampled-moment low-energy problem (large structures): moments of the PSD-passing members of uniform samples,
    tau from the empirical energy quantile of those samples."""
    st = Structure(v, [gen])
    rng = np.random.default_rng(seed)
    def pool(blk, r, bound, n):
        Y = blk.sample(r, n, rng)
        F = blk.values(Y)
        ok = (st.psd(F) <= bound + 1e-6).all(1)
        return F[ok]
    FA = pool(st.A, rows[0], 4 * v, min(samples, 200000))
    FB = pool(st.B, rows[1], 2 * v, samples)
    FD = pool(st.B, rows[2], 4 * v, samples)
    print(f'PSD-passing samples: A {len(FA)}, B {len(FB)}, D {len(FD)}', flush=True)
    em = energy_model(st, moments_of(FA), moments_of(FB), moments_of(FD))
    EB, ED = energies(FB, *em['B']), energies(FD, *em['D'])
    tB, tD = float(np.quantile(EB, fB)), float(np.quantile(ED, fD))
    P.write_problem(st, rows, out, canonB=True, canonD=True, psd=(True, True, True), seed=v,
                    energy=dict(B=(tB, *em['B']), D=(tD, *em['D'])))
    print(f'wrote {out}: tau_B {tB:.2f} (keeps {fB}), tau_D {tD:.2f} (keeps {fD}) of the PSD-passing pools')
    return em, tB, tD


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('cmd', choices=['study', 'checkc', 'plan'])
    ap.add_argument('v', nargs='?', type=int); ap.add_argument('gen', nargs='?', type=int)
    ap.add_argument('--rows', default=None); ap.add_argument('--samples', type=int, default=400000)
    ap.add_argument('--fB', type=float, default=0.05); ap.add_argument('--fD', type=float, default=0.05)
    ap.add_argument('--out', default='problem_lowe.txt')
    a = ap.parse_args()
    if a.cmd == 'study':
        study()
    elif a.cmd == 'checkc':
        sys.exit(0 if checkc(tempfile.mkdtemp(prefix='p3low_')) else 1)
    else:
        plan(a.v, a.gen, tuple(int(x) for x in a.rows.split(',')), a.samples, a.fB, a.fD, a.out)
