# Independent census of the plain cyclic Spence shapes at v = 2p-1: all subgroups H of Z_v^* (own enumeration: joins of cyclic subgroups), orbit counts, row-sum feasibility (+1)
# (independent code of the judge of search/NOTES.md section 92).  usage: plain_shapes.py OMAX v [v ...]
import math, sys, itertools, collections, json
import numpy as np
def units(v): return [u for u in range(1,v) if math.gcd(u,v)==1]
def cyc(g,v):
    S=[1]; x=g
    while x!=1: S.append(x); x=x*g%v
    return frozenset(S)
def mul_closure(A,B,v):
    return frozenset(a*b%v for a in A for b in B)
def all_subgroups(v):
    U=units(v)
    C={cyc(g,v) for g in U}
    subs=set(C); frontier=set(C)
    while frontier:
        new=set()
        for A in frontier:
            for B in C:
                if not B<=A:
                    J=mul_closure(A,B,v)     # product of two subgroups of an abelian group is a subgroup
                    if J not in subs: new.add(J)
        subs|=new; frontier=new
    return subs
def orbit_sizes(v,H):
    seen=[False]*v; sizes=[]; lab=[-1]*v
    for g in range(v):
        if not seen[g]:
            o=set(g*h%v for h in H)
            for x in o: seen[x]=True; lab[x]=len(sizes)
            sizes.append(len(o))
    return sizes,lab
def rowsum_count(sizes):
    # number of sign vectors x with sum size*x = +1 (DP over sums)
    dp={0:1}
    for s in sizes:
        nd=collections.defaultdict(int)
        for t,c in dp.items(): nd[t+s]+=c; nd[t-s]+=c
        dp=nd
    return dp.get(1,0)
if __name__=="__main__":
    Omax=int(sys.argv[1]); out={}
    for v in map(int,sys.argv[2:]):
        subs=all_subgroups(v); rows=[]
        for H in subs:
            sizes,lab=orbit_sizes(v,H); O=len(sizes)
            if 2<=O<=Omax: rows.append((O,len(H),rowsum_count(sizes),min(H-{1}) if len(H)>1 else 1))
        rows.sort()
        print('v=%d: %d subgroups of Z_v^*, %d shapes with 2<=O<=%d; with L>0 (a +1 row-sum vector exists): %d'%(v,len(subs),len(rows),Omax,sum(1 for r in rows if r[2]>0)),flush=True)
        for r in rows:
            if r[2]>0: print('    O=%d |H|=%d L=%d'%(r[0],r[1],r[2]))
