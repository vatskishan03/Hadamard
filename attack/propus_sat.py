"""Propus structures as CP-SAT models (OR-tools): one Boolean per orbit (1 = entry -1), the row sums as linear
constraints, and every class equation PAF_A + 2 PAF_B + PAF_D = 0 at a representative shift s as
    sum_X w_X sum_{orbit pairs} m_X(o1, o2; s) [x_o1 XOR x_o2] = 2v        (w = 1, 2, 1),
since PAF_X(s) = v - 2 #{g : X(g) != X(g + s)}. A is constant on the orbits of <H, -1>, B and D on those of H.

  python3 attack/propus_sat.py v gens ra,rb,rd [--time 600] [--workers 4] [--seed 1]

A solution is checked on Z_v (PAF sum 4v delta) and the symmetric Hadamard matrix is rebuilt and verified; it is
printed as a 'HIT verified' line for attack/hit2rec.py."""
import argparse, sys, time
import numpy as np
from pathlib import Path as _P
sys.path[:0] = [str(_P(__file__).resolve().parent), str(_P(__file__).resolve().parent.parent / 'verify')]
import propus as PR
import verify
from ortools.sat.python import cp_model


def build(v, gens, rows):
    st_ad, st_b = PR.structures(v, gens, False)
    m = cp_model.CpModel()
    X = {}
    for name, st in (('A', st_ad), ('B', st_b), ('D', st_b)):
        X[name] = [m.NewBoolVar(f'{name}{o}') for o in range(st.O)]
    sts = {'A': st_ad, 'B': st_b, 'D': st_b}
    for (name, r) in zip('ABD', rows):
        st = sts[name]
        m.Add(sum(int(st.sizes[o]) * (1 - 2 * X[name][o]) for o in range(st.O)) == r)   # row sum
    xor_cache = {}

    def xor(name, o1, o2):
        k = (name, min(o1, o2), max(o1, o2))
        if k not in xor_cache:
            z = m.NewBoolVar(f'x{name}_{k[1]}_{k[2]}')
            m.AddBoolXOr([X[name][o1], X[name][o2], z.Not()])
            xor_cache[k] = z
        return xor_cache[k]

    for c, s in enumerate(st_ad.reps):
        terms = []
        for name, w in (('A', 1), ('B', 2), ('D', 1)):
            orb = sts[name].orb
            cnt = {}
            for g in range(v):
                o1, o2 = int(orb[g]), int(orb[(g + s) % v])
                if o1 == o2: continue
                k = (min(o1, o2), max(o1, o2))
                cnt[k] = cnt.get(k, 0) + 1
            for (o1, o2), mlt in cnt.items():
                terms.append(w * mlt * xor(name, o1, o2))
        m.Add(sum(terms) == 2 * v)
    return m, X, sts


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('v', type=int); ap.add_argument('gens'); ap.add_argument('rows')
    ap.add_argument('--time', type=float, default=600); ap.add_argument('--workers', type=int, default=4)
    ap.add_argument('--seed', type=int, default=1)
    a = ap.parse_args()
    v, gens, rows = a.v, [int(x) for x in a.gens.split(',')], [int(x) for x in a.rows.split(',')]
    t = time.time()
    m, X, sts = build(v, gens, rows)
    print(f'v={v} H={gens} rows={rows}: orbits A {sts["A"].O} B {sts["B"].O}; '
          f'{len(m.Proto().variables)} variables, {len(m.Proto().constraints)} constraints [{time.time() - t:.1f}s]', flush=True)
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = a.time
    solver.parameters.num_search_workers = a.workers
    solver.parameters.random_seed = a.seed
    st = solver.Solve(m)
    print('status', solver.StatusName(st), f'[{time.time() - t:.1f}s]', flush=True)
    if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        seq = {}
        for name in 'ABD':
            vals = np.array([1 - 2 * solver.Value(x) for x in X[name]], dtype=np.int64)
            seq[name] = vals[sts[name].orb]
        tot = verify.paf_sum([seq['A'], seq['B'], seq['B'], seq['D']])
        ok = tot[0] == 4 * v and not tot[1:].any()
        Hm = PR.symmetric_matrix_v3(seq['A'], seq['B'], seq['D'])
        ok2 = np.array_equal(Hm, Hm.T) and np.array_equal(Hm @ Hm.T, 4 * v * np.eye(4 * v, dtype=Hm.dtype))
        print(f'HIT {"verified" if ok and ok2 else "WRONG"} A ' + ''.join('+' if x > 0 else '-' for x in seq['A']) +
              ' B ' + ''.join('+' if x > 0 else '-' for x in seq['B']) + ' D ' + ''.join('+' if x > 0 else '-' for x in seq['D']), flush=True)


if __name__ == '__main__':
    main()
