# Independent DC (row-sum / Z_2-compression) census for the plain multiplier structures at p = 599, 719, 739, 839, 1031, 1087, for ALL subgroup orders d | p-1 (odd and even;
# ../census/c8p_basic.py lists the odd orders only).  Real Goethals-Seidel over Z_2p and quaternary (negacyclic) over Z_p.  Independent code (search/NOTES.md section 92).
import math
from sympy import divisors
def real_live(p,d):
    c=(p-1)//d
    vals=sorted({u+d*t for u in (1,-1) for t in range(-c,c+1,2)})
    B=8*p
    pr={}
    for s0 in vals:
        if s0*s0>B: continue
        for s1 in vals:
            q=s0*s0+s1*s1
            if q<=B: pr.setdefault((q,s0*s1),0); pr[(q,s0*s1)]+=1
    keys=list(pr)
    two=set()
    for a in keys:
        for b in keys:
            if a[0]+b[0]<=B: two.add((a[0]+b[0],a[1]+b[1]))
    return any(((B-x,-y) in two) for (x,y) in two)
def quat_live(p,d):
    c=(p-1)//d; B=4*p
    R2=set()
    for (ur,ui) in ((1,0),(-1,0),(0,1),(0,-1)):
        for a in range(-c,c+1):
            for b in range(-c,c+1):
                if abs(a)+abs(b)>c or (a+b-c)%2: continue
                r=(ur+d*a)**2+(ui+d*b)**2
                if r<=B: R2.add(r)
    two={a+b for a in R2 for b in R2 if a+b<=B}
    return any((B-s) in two for s in two)
for p in (599,719,739,839,1031,1087):
    print('p=%d  p-1=%d'%(p,p-1))
    for d in divisors(p-1):
        if d==1: continue
        c=(p-1)//d
        r=real_live(p,d); q=quat_live(p,d)
        print('   |H|=%4d (%s)  orbits: real Z_2p %4d, quaternary Z_p %4d ;  real DC %s, quaternary DC %s'%(d,'even' if d%2==0 else 'odd ',2+2*c,1+c,'LIVE' if r else 'dead','LIVE' if q else 'dead'))
