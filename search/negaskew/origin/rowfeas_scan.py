import sys, os
sys.path.insert(0,'.')
from rowfeas import *
from oo_shapes import shapes, combos
from sympy import isprime
sh=shapes(13,int(sys.argv[2]),int(sys.argv[1]),5)
sh.sort(key=lambda s:(s['nu'],s['v']))
for s in sh:
    m,q,h,nu=s['m'],s['q'],s['hsize'],s['nu']
    if not isprime(m): continue
    p=(nu-1)//2
    ok0=ok1=False; both=[]
    for cb in combos(m,q):
        r0=feasible(m,q,h,p,cb,0); r1=feasible(m,q,h,p,cb,1)
        if all(r0): ok0=True
        if all(r1): ok1=True
        both.append((cb,r0,r1))
    print('nu',nu,'m',m,'q',q,'h',h,'p',p,'rowsum-feasible z=0:',ok0,' z=1:',ok1, ' #combos',len(both))
