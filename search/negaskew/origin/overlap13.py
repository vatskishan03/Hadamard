import sys
sys.path.insert(0,'.')
from common import *
def key(sol, mode):
    bl=[]
    for b in sol:
        if mode=='x0x1x2': bl.append((b[0][1:],b[1][1:],b[2][1:]))
        elif mode=='x0': bl.append(b[0][1:])
        elif mode=='x0x1': bl.append((b[0][1:],b[1][1:]))
    return tuple(sorted(bl))
files={'299':'runs/299_all.out','247':'runs/13_19_all.out','351':'runs/13_27_all.out'}
S={}
for nm,fn in files.items():
    sols=load_solutions(fn)
    S[nm]=sols
    print(nm,len(sols))
for mode in ('x0x1x2','x0x1','x0'):
    K={nm:set(key(s,mode) for s in sols) for nm,sols in S.items()}
    print(mode,{nm:len(k) for nm,k in K.items()},'  overlaps 299&247',len(K['299']&K['247']),'299&351',len(K['299']&K['351']),'247&351',len(K['247']&K['351']))
# block-level overlap: individual (x0,x1,x2) off-origin blocks
B={nm:set((b[0][1:],b[1][1:],b[2][1:]) for s in sols for b in s) for nm,sols in S.items()}
print('distinct off-origin blocks',{nm:len(b) for nm,b in B.items()},'overlap 299&351',len(B['299']&B['351']),'299&247',len(B['299']&B['247']))
