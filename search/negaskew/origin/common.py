import sys, os, json, itertools, math, cmath
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from qshape import kshape
U = [1, 1j, -1, -1j]

def load_solutions(fn):
    sols = []
    for line in open(fn):
        if not line.startswith('SOL'): continue
        blocks = line.split()[1:]
        sol = []
        for b in blocks:
            x0, x1, x2 = b.split('|')
            sol.append((tuple(map(int, x0)), tuple(map(int, x1)), tuple(map(int, x2))))
        sols.append(tuple(sol))
    return sols

class Shape:
    def __init__(self, m, hK):
        self.m = m
        self.sh = kshape([m], hK)
        self.orb = self.sh['orb']
        self.nu = self.sh['nu']
        self.ell = self.sh['ell']
    def expand(self, s):
        return np.array([U[s[self.orb[k]]] for k in range(self.m)])

def paf(x):
    m = len(x)
    return np.array([np.sum(x * np.conj(np.roll(x, -k))) for k in range(m)])
