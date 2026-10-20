"""Exhaustive Spence sweep: G = Z_2 x Z_m (m = (n-1)/2, n = 3 mod 4 open orders), H <= Z_m^*,
slice-admissible structures with small pools; find all (2,0,0,0) quadruples with find4 and
verify by building the bordered matrix (verify.build_spence) and checking H H^T = nI."""
import os, re, struct, subprocess, sys, tempfile, time, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'verify'))
import spence_lib as L, count_exact as X, spence_slice_check as SC
import verify

H4 = [[1,1,1,1],[1,-1,1,-1],[1,1,-1,-1],[1,-1,-1,1]]
P = (np.ones((4,4), int) - 2*np.eye(4, dtype=int)).tolist()

def search(n, H, cap=3e10, K=3):
    m = (n - 1) // 2
    G = L.AbGroup([2, m])
    st = L.Structure(G, [G.perm_from_matrix(np.array([[1, 0], [0, h]])) for h in H])
    nA, nB = st.pool(2), st.pool(0)
    if nA == 0 or nB == 0: return None, st, 0
    work = nA * nB + nB * nB
    if work > cap: return None, st, work
    pm = st.G.v % 4; target = (-4 - 4 * pm) // 4
    A, B = X.pool_vectors(st, 2), X.pool_vectors(st, 0)
    VA = np.rint((st.class_values(A) - pm) / 4).astype(np.int16)
    VB = np.rint((st.class_values(B) - pm) / 4).astype(np.int16)
    fd, path = tempfile.mkstemp(suffix=".bin", dir=HERE); os.close(fd)
    try:
        with open(path, "wb") as f:
            f.write(struct.pack("iiiii", len(VA), len(VB), len(VB), len(VB), st.C))
            f.write(np.full(st.C, target, dtype=np.int16).tobytes())
            for V in (VA, VB, VB, VB): f.write(np.ascontiguousarray(V).tobytes())
        out = subprocess.run([os.path.join(HERE, "find4b"), path, str(K)], capture_output=True, text=True, check=True).stdout
    finally:
        os.remove(path)
    found = []
    for line in out.splitlines():
        if line.startswith("SOL"):
            ia, ib, ic, id_ = map(int, line.split()[1:])
            seqs = [Y[i].astype(np.int64)[st.orb] for Y, i in ((A, ia), (B, ib), (B, ic), (B, id_))]
            cert = {"group": [2, m], "sequences": [''.join('+' if x > 0 else '-' for x in s) for s in seqs], "P": P, "Q": H4}
            Hm, _ = verify.build_spence(cert)
            N = Hm.shape[0]; Hf = Hm.astype(np.float64)
            found.append((cert, bool(np.isin(Hm, (-1, 1)).all() and np.array_equal(Hf @ Hf.T, N * np.eye(N)))))
    return found, st, work

if __name__ == "__main__":
    fr = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'FRONTIER.md')).read()
    ns = [int(n) for _, n in re.findall(r'^\|\s*(\d+)\s*\|\s*(\d+)\s*\|', fr, re.M)]
    MAXO, CAP = int(sys.argv[1]), float(sys.argv[2])
    log = open(os.path.join(HERE, 'spence_sweep_hits.json'), 'a')
    for n in ns:
        if n % 4 != 3: continue
        m = (n - 1) // 2; t0 = time.time(); tried = 0
        for H in sorted(SC.subgroups(m), key=len, reverse=True):
            sizes = SC.orbit_sizes(m, H)
            if 2 * len(sizes) > MAXO: continue
            ok, wit = SC.admissible(n, sizes)
            if not ok: continue
            found, st, work = search(n, H, cap=CAP)
            if found is None: continue
            tried += 1
            for cert, had in found:
                print(f"*** n={n} |H|={len(H)} orbits={st.O}: Spence solution, H({4*n}) Hadamard={had}", flush=True)
                log.write(json.dumps(dict(n=n, H=sorted(H), cert=cert, hadamard=had)) + "\n"); log.flush()
        print(f"n={n}: {tried} admissible small structures searched ({time.time()-t0:.0f}s)", flush=True)
