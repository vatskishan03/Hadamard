"""Exhaustive search of a Goethals-Seidel multiplier structure with small pools (all patterns)."""
import os, struct, subprocess, sys, tempfile, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "verify"))
import spence_lib as L, gs_lib as GS, count_exact as X
import verify

HERE = os.path.dirname(os.path.abspath(__file__))

def run(n, gens_ints, max_pairs=4e10, K=5, verbose=True):
    G = L.AbGroup([n])
    st = L.Structure(G, [G.perm_from_matrix(np.array([[g]])) for g in gens_ints])
    pm = n % 4; target = -pm
    pats = [q for q in GS.patterns(n) if all(st.pool(r) > 0 for r in q)]
    if verbose: print(f"n={n} gens={gens_ints} |orbits|={st.O} classes={st.C} patterns={len(pats)}", flush=True)
    cache, found = {}, []
    def vecs(r):
        if r not in cache:
            Y = X.pool_vectors(st, r)
            V = np.rint((st.class_values(Y) - pm) / 4).astype(np.int16)
            cache[r] = (Y.astype(np.int8), V)
        return cache[r]
    for q in pats:
        sizes = [st.pool(r) for r in q]
        # pair the lists to balance work: (largest, smallest) and the middle two
        order = sorted(range(4), key=lambda i: sizes[i])
        pairing = [order[3], order[0], order[2], order[1]]
        work = sizes[pairing[0]] * sizes[pairing[1]] + sizes[pairing[2]] * sizes[pairing[3]]
        if work > max_pairs:
            if verbose: print(f"  pattern {q}: pools {[int(s) for s in sizes]} work {work:.2e} > cap, skipped", flush=True)
            continue
        lists = [vecs(q[i]) for i in pairing]
        fd, path = tempfile.mkstemp(suffix=".bin", dir=HERE); os.close(fd)
        try:
            with open(path, "wb") as f:
                f.write(struct.pack("iiiii", *[len(l[1]) for l in lists], st.C))
                f.write(np.full(st.C, target, dtype=np.int16).tobytes())
                for l in lists: f.write(np.ascontiguousarray(l[1]).tobytes())
            t = time.time()
            out = subprocess.run([os.path.join(HERE, "find4b" if os.path.exists(os.path.join(HERE, "find4b")) else "find4"), path, str(K)], capture_output=True, text=True, check=True).stdout
        finally:
            os.remove(path)
        total = int(out.strip().split("TOTAL")[-1])
        if verbose: print(f"  pattern {q}: pools {[int(s) for s in sizes]} work {work:.2e}: {total} solutions [{time.time()-t:.1f}s]", flush=True)
        for line in out.splitlines():
            if line.startswith("SOL"):
                idx = list(map(int, line.split()[1:]))
                seqs = [None] * 4
                for slot, li in enumerate(pairing):
                    y = lists[slot][0][idx[slot]].astype(np.int64)
                    seqs[li] = y[st.orb]
                tot = verify.paf_sum(seqs)
                ok = tot[0] == 4 * n and not tot[1:].any()
                if ok:
                    H = verify.goethals_seidel(*map(verify.circulant, seqs))
                    Hf = H.astype(np.float64)
                    had = np.array_equal(Hf @ Hf.T, 4 * n * np.eye(4 * n))
                    found.append((q, seqs, had))
                    if verbose: print(f"    VERIFIED quadruple: PAF ok, H({4*n}) Hadamard = {had}", flush=True)
    return found, st

if __name__ == "__main__":
    n = int(sys.argv[1]); gens = list(map(int, sys.argv[2:]))
    found, st = run(n, gens)
    print(f"done: {len(found)} verified solutions")
