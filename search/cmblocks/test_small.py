import numpy as np, sys
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
import cmb, paley_circ as P
for qP in (11, 19, 43):
    Kr, Ki, m = P.block_circulant_K(qP)
    q = 2 * m + 1
    assert q % 4 == 3
    print("K from Paley", qP, "order", 2 * m, "-> Miyamoto q =", q, "prime?", P.is_prime(q))
    if not P.is_prime(q): continue
    g = cmb.prim_root(q)
    N = cmb.build_blocks(q, g, Kr, Ki)
    cmb.check_blocks(N)
