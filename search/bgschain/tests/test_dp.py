"""DP (orbit vs grouped), fibre enumeration and class values on real structures with a planted sign vector: every fibre contains the planted vector."""
import sys, os, json, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import numpy as np
import hlib

def run(cfgname, towers, seeds=3):
    cfg = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'configs', cfgname)))
    P = hlib.Problem(cfg); ok = True
    rng = np.random.default_rng(11)
    for tower in towers:
        levels = [hlib.Level(P, k) for k in tower[:-1]]; new = hlib.Level(P, tower[-1]) if tower[-1] is not None else None
        for s in range(seeds):
            # planted sign vector with the pattern's coset sums is not needed: use its own coset sums as y
            for _ in range(1000):
                Y = 2 * rng.integers(0, 2, P.O) - 1
                mask = int(sum(1 << o for o in range(P.O) if Y[o] > 0))
                yv = P.Qlevel.of_masks(np.array([mask], dtype=np.uint64), P.O)[0]
                break
            fixed = [(l, l.of_masks(np.array([mask], dtype=np.uint64), P.O)[0]) for l in levels]
            t = time.time()
            a = P.dp(yv, fixed, new, strategy='orbit', max_states=40_000_000)
            ta = time.time() - t; t = time.time()
            b = P.dp(yv, fixed, new, strategy='grouped', max_states=40_000_000)
            tb = time.time() - t
            if new is None: same = np.allclose(a[1], b[1])
            else: same = np.array_equal(a[0], b[0]) and np.allclose(a[1], b[1])
            contains = True
            if new is not None:
                tv = new.of_masks(np.array([mask], dtype=np.uint64), P.O)[0].astype(np.int16)
                contains = bool((a[0] == tv[None, :]).all(axis=1).any())
            cnt = float(a[1].sum())
            # enumeration of the fibre (all levels fixed) against the DP count
            en = None
            if cnt < 5e6:
                fx2 = fixed + ([(new, new.of_masks(np.array([mask], dtype=np.uint64), P.O)[0])] if new is not None else [])
                m = P.enumerate_fibre(yv, fx2)
                c2 = P.fibre_count(yv, fx2)
                en = (len(m) == round(c2)) and bool((m == np.uint64(mask)).any())
                # class values of the enumerated masks vs direct sums
                if len(m):
                    V = P.class_values(m[:50]); seqs = [P.sequence(int(x)) for x in m[:50]]
                    for j, x in enumerate(seqs[:20]):
                        for c, d in enumerate(P.reps):
                            sh = P.G.translate(d)
                            if int(x @ x[sh]) != V[j, c]: en = False
            print(f"{cfgname} tower {tower} seed {s}: states {len(a[0]) if a[0] is not None else 1} peak {a[2]}/{b[2]} orbit {ta:.2f}s grouped {tb:.2f}s equal {same} contains-planted {contains} enum-ok {en}", flush=True)
            ok &= same and contains and (en is None or en)
    return ok

if __name__ == '__main__':
    ok = True
    ok &= run('n1039.json', [[[3], [0, 1, 3]], [[3], [0, 1, 3], [2]], [[3], [0, 1, 3], [2], [0, 1, 2]], [[2], [3]], [[3], [2]]], seeds=2)
    ok &= run('n3964.json', [[[2], [0, 1, 2]], [[2], [0, 1, 2], [0, 1, 3]], [[2], [0, 1, 2], [0, 1, 3], [2, 3]]], seeds=2)
    print('DP/enumeration tests:', 'ALL OK' if ok else 'FAILURES')
