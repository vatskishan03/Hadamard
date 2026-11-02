"""Gaussian estimate E (expected number of solutions) for 239/<10>/(15,1,27) with the repository's own estimator
(search/propus_rank.estimate, imported read-only), two seeds; and the orbit count E/|Gamma|."""
import os, sys, math
os.environ['HADAMARD_WORK'] = str(__import__('pathlib').Path(__file__).resolve().parent / 'work' / 'hw')
_R = __import__('pathlib').Path(__file__).resolve().parent.parent     # search/ when this directory is search/propus3
sys.path[:0] = [str(_R), str(_R.parent / 'verify')]
import numpy as np
import propus as PR, propus_rank as PRk
st_ad, st_b = PR.structures(239, [10], False)
for seed in (1, 2):
    rng = np.random.default_rng(seed)
    lE, rk = PRk.estimate(st_ad, st_b, st_b, (15, 1, 27), 15000, rng, st_ad.csize)
    print(f'seed {seed}: log10 E = {lE:.2f}  (E = {10**lE:.1f} solutions), rank {rk}; |Gamma| = 34*2 = 68 -> {10**lE/68:.2f} expected orbits; P(>=1 orbit) = {1-math.exp(-10**lE/68):.2f} if Poisson')
