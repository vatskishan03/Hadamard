"""Row-sum (DC) feasibility of the Paley shapes Z_m x F_q with H = {1} x squares of F_q* (H_K trivial, G = Z_{m q}, the generator by the Chinese remainder
theorem): the rigid candidates 473 = 11 * 43, 517 = 11 * 47, 1139 = 17 * 67 and 1207 = 17 * 71 of NOTES.md section 88 have a row-sum pattern
(the quaternary class has no "2m - 1 is not a sum of two squares" obstruction), as have 299, 589, 413 and 575.
usage: python3 paley_dc.py      prints v, m, q, the number of orbits n and 'DC feasible' (dcnp.plain_dc) for the nine shapes in the list below"""
import sys
from sympy import primitive_root
from qdrv import shape
from dcnp import plain_dc
def crt_gen(v,m,q):
    g=primitive_root(q); sq=g*g%q
    for h in range(2,v):
        if h%m==1 and h%q==sq: return h
for (m,q) in [(11,43),(11,47),(17,67),(17,71),(13,23),(19,31),(7,59),(25,23),(5,7)]:
    v=m*q
    try:
        h=crt_gen(v,m,q)
        d=shape(v,[h])
        print(v,m,q,'n',d['n'],'DC feasible',plain_dc(v,d['sz'],d['neg'],d['n']),flush=True)
    except Exception as ex: print(v,m,q,'err',ex)
