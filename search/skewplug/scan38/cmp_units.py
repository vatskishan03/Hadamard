import sys, glob, math, numpy as np
"""the check of the unit reduction of the skew block in exh38.py (EX_UNITS): the classes of the level-1 instances of one slot with and without it generate the same unit orbits.  Needs the dumps of two runs of the slot (2,2,4,7)/(2,1,3,4) of 385 (K = Z35, q = 11, H_K = <11>) with EX_DUMP=dump and EX_UNITS=0 (tag tstA) and EX_UNITS=1 (tag tstB).  usage (in search/skewplug): python3 scan38/cmp_units.py"""
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from sk import SShape
sh=SShape([35],11,[[11]]); n=35
_Wh = np.random.default_rng(12345).integers(1, 2 ** 62, size=n, dtype=np.uint64) * np.uint64(2) + np.uint64(1)
_fix = [int(P[0]) for P in sh.kCD.pos if len(P) == 1]; _neg = sh.negperm
_permsB = [np.asarray(sh.shift[t], dtype=np.int64)[_neg] for t in _fix] + [np.asarray(sh.shift[t], dtype=np.int64) for t in _fix]
_permsA = [np.arange(n), np.asarray(_neg, dtype=np.int64)]
_MUL = [np.uint64(0x9E3779B97F4A7C15), np.uint64(0xC2B2AE3D27D4EB4F), np.uint64(0x165667B19E3779F9), np.uint64(0x27D4EB2F165667C5)]
def _canon(Sx, Ex, perms, negs):
    best = None
    for sg in ((1, -1) if negs else (1,)):
        v = ((sg * Sx.astype(np.int64) + 2) * 3 + (sg * Ex.astype(np.int64) + 1) + 1).astype(np.uint64)
        for P in perms:
            h = (v[:, P] * _Wh).sum(axis=1, dtype=np.uint64)
            best = h if best is None else np.minimum(best, h)
    return best
def inst_keys(S, E):
    key = np.zeros(len(S), dtype=np.uint64)
    for X in range(4):
        key = key + _canon(S[:, X, :], E[:, X, :], _permsA if X == 0 else _permsB, X > 0) * _MUL[X]
    return key
units=[(u*np.arange(n))%n for u in range(1,n) if math.gcd(u,n)==1]
def load(pref):
    Ss=[];Es=[]
    for f in sorted(glob.glob(f'dump/dump_{pref}_*.npz')):
        z=np.load(f); Ss.append(z['S']); Es.append(z['E'])
    return np.concatenate(Ss),np.concatenate(Es)
def unitkey(S,E):
    best=None
    for P in units:
        k=inst_keys(S[:,:,P],E[:,:,P]); best=k if best is None else np.minimum(best,k)
    return best
SA,EA=load('tstA'); SB,EB=load('tstB')
kA=inst_keys(SA,EA); kB=inst_keys(SB,EB)
print('A: raw',len(SA),'classes',len(set(kA.tolist())),' B(units): raw',len(SB),'classes',len(set(kB.tolist())))
uA=set(unitkey(SA,EA).tolist()); uB=set(unitkey(SB,EB).tolist())
print('unit-orbit classes: A',len(uA),'B',len(uB),'equal',uA==uB,'A-B',len(uA-uB),'B-A',len(uB-uA))
