"""Write the orbit-tabu input for GS over Z_v with multiplier generators; verify a solution line."""
import sys, struct, numpy as np
sys.path.insert(0, '/home/user/Hadamard/attack'); sys.path.insert(0, '/home/user/Hadamard/verify')
import spence_lib as L
def build(v, gens):
    G = L.AbGroup([v]); st = L.Structure(G, [G.perm_from_matrix(np.array([[g]])) for g in gens])
    S = st.S.astype(np.int64)                                  # (C, k, k)
    M = S + S.transpose(0, 2, 1)
    for c in range(st.C): np.fill_diagonal(M[c], 0)
    Dg = np.array([int(np.trace(S[c])) for c in range(st.C)])
    return st, M, Dg
def write(path, st, M, Dg):
    with open(path, 'wb') as f:
        f.write(struct.pack('ii', st.O, st.C))
        f.write(np.asarray(st.csize, dtype=np.int32).tobytes())
        f.write(M.astype(np.int32).tobytes()); f.write(Dg.astype(np.int32).tobytes())
if __name__ == '__main__':
    v = int(sys.argv[1]); gens = [int(x) for x in sys.argv[2].split(',')]; out = sys.argv[3]
    st, M, Dg = build(v, gens); write(out, st, M, Dg)
    print(f"v={v} gens={gens} orbits={st.O} classes={st.C} sizes={sorted(set(st.sizes.tolist()))}")
