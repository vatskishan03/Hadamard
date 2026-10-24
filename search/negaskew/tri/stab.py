"""Stabilisers in the group Z_e x tilde^3 x S_3 x (A -> -A) (PROOFS.md P5) of a sample of at most 2500 solutions (seed 2), by type: shift != 0, sign change of A, number of tildes, permutation != identity.
usage: python3 stab.py v k FILE        (FILE: SOL lines of ../qskew or of tri2)"""
import sys, os, itertools; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from trilib import *
v,k=int(sys.argv[1]),int(sys.argv[2]); fn=sys.argv[3]
d,sols=load_sols(fn,v,k); e=d['e']; h=e//2
def shift(b,a): return (b[0],tuple(b[1][(j+a)%e] for j in range(e)))
def tilde(b): return (b[0], tuple((-b[1][(j+h)%e])%4 for j in range(e)))
def negA(b): return ((3 if b[0]==0 else 0), tuple((c+2)%4 for c in b[1]))
def image(s,a,na,tl,perm):
    A=shift(s[0],a); A=negA(A) if na else A
    fr=[shift(x,a) for x in s[1:]]
    fr=[tilde(x) if t else x for x,t in zip(fr,tl)]
    fr=[fr[p] for p in perm]
    return (A,)+tuple(fr)
# find solutions with nontrivial stabilizer and describe the stabilizer elements
import random; random.seed(2); sols=random.sample(sols,min(2500,len(sols)))
cnt=collections.Counter(); ex={}
seen=set()
for s in sols:
    if s in seen: continue
    stab=[]
    for a in range(e):
        for na in (0,1):
            for tl in itertools.product((0,1),repeat=3):
                for perm in itertools.permutations(range(3)):
                    if (a,na,tl,perm)==(0,0,(0,0,0),(0,1,2)): continue
                    if image(s,a,na,tl,perm)==s: stab.append((a,na,tl,perm))
    # (orbit marking skipped: report only first few)
    if stab:
        key=tuple(sorted(set((a!=0,na,sum(tl),perm!=(0,1,2)) for a,na,tl,perm in stab)))
        cnt[key]+=1; ex.setdefault(key,(s,stab))
    if sum(cnt.values())>400: break
for kk,c in cnt.most_common(): print(c,'solutions with stabilizer types (shift!=0, negA, #tilde, perm!=id):',kk,'example stab elems',ex[kk][1][:3])
print('scanned until',sum(cnt.values()),'stabilized solutions')
