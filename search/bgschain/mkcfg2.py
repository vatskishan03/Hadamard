"""General builder: G = Z4^3 x (odd group given by moduli), H = <(T3, Hodd)>, S = 2Z4^3 x odd, k = 8, t_in = -28, t_out = 4."""
import numpy as np
T3 = [[2, 2, 1], [1, 0, 0], [2, 3, 3]]

def build(odd_mods, Hodd, A, a, tower='Z4,2Z4', split='Odd', bdepth=1, name='z4cube', **extra):
    r = 3 + len(odd_mods)
    H = [[0] * r for _ in range(r)]
    for i in range(3):
        for j in range(3): H[i][j] = T3[i][j]
    for i in range(len(odd_mods)):
        for j in range(len(odd_mods)): H[3 + i][3 + j] = int(Hodd[i][j]) if hasattr(Hodd, '__getitem__') and hasattr(Hodd[0], '__getitem__') else (int(Hodd) if i == j else 0)
    e = lambda i, c=1: [c if t == i else 0 for t in range(r)]
    KER = {'Z4': [e(i) for i in range(3)], '2Z4': [e(i, 2) for i in range(3)], 'Odd': [e(3 + j) for j in range(len(odd_mods))],
           'S': [e(i, 2) for i in range(3)] + [e(3 + j) for j in range(len(odd_mods))]}
    cfg = dict(name=name, mods=[4, 4, 4] + list(odd_mods), S_gens=KER['S'], gens=[H], t_in=-28, t_out=4,
               ys=[[A[i]] + [a[i]] * 7 for i in range(4)],
               levels=[{'ker': KER[w], 'name': w} for w in tower.split(',') if w], split=[{'ker': KER[w], 'name': w} for w in split.split(',') if w],
               bdepth=bdepth, sym_swap=True)
    cfg.update(extra)
    return cfg
