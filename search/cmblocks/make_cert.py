"""Write a certificate (verify/verify.py construction quaternary-miyamoto-blocks, without its sha256) for q = 2m + 1,
m = (p + 1)/4 with K from Paley's quaternary matrix (p = 3 mod 8 prime), or q = p + 2 with K = C + iI (p = 1 mod 4,
source "conf"), and T-matrices of order t:  python3 make_cert.py p t OUT.json [paley|conf]"""
import json, sys, numpy as np
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
import cmb, paley_circ as P
def exps(Z):
    return [[int({1: 0, 1j: 1, -1: 2, -1j: 3}[complex(round(z.real), round(z.imag))]) for z in row] for row in Z]
def digits(row_r, row_i):
    e = np.where(row_r == 1, 0, np.where(row_i == 1, 1, np.where(row_r == -1, 2, 3)))
    return "".join(map(str, e))
def make(qP, t, out, source="paley"):
    if source == "paley":
        Kr, Ki, m = P.block_circulant_K(qP, check=(qP < 400))
        src = f"Paley quaternary Hadamard matrix of order (p+1)/2, p = {qP}, in 2x2 block circulant form (paley_circ.py)"
    else:
        import conf_circ
        Kr, Ki, m = conf_circ.block_circulant_K_conf(qP)
        src = f"K = C + iI, C Paley's symmetric conference matrix of order p+1, p = {qP}, in 2x2 block circulant form (conf_circ.py)"
    q = 2 * m + 1
    A, B = cmb.border_pair()
    T = cmb.t_matrices(t)
    cert = {"order": 8 * q * t, "construction": "quaternary-miyamoto-blocks", "q": q, "g": cmb.prim_root(q),
            "K": {f"{a}{b}": digits(Kr[a * m, b * m:(b + 1) * m], Ki[a * m, b * m:(b + 1) * m]) for a in range(2) for b in range(2)},
            "K_source": src, "A": exps(A), "B": exps(B),
            "T": ["".join({1: "+", -1: "-", 0: "0"}[int(x)] for x in r) for r in T]}
    json.dump(cert, open(out, "w"), indent=1)
    print("wrote", out)
if __name__ == "__main__":
    make(int(sys.argv[1]), int(sys.argv[2]), sys.argv[3], sys.argv[4] if len(sys.argv) > 4 else "paley")
