"""Validation of the symmetry reductions of pal_lists.py: runs pal_lists on a shape with all classes and compares the total count,
weighted by 2 for every halved block (sigma or nu) that is not fixed, with the full count of qskew.
usage: python3 validate_sym.py v K q "[gens]" tag expected
Example: python3 validate_sym.py 77 7 11 "[36,23]" v77s 3840     (Z_7 x F_11, H_K of order 3; 5 s; the count of qskew2 is 3840)
Prints 'weighted total N expected M OK' (or MISMATCH).  Needs ./qjoin."""
import sys, json, os, subprocess
import numpy as np
sys.path.insert(0,'.')
from pal_lists import *
v,K,q=int(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3]); gens=json.loads(sys.argv[4]); tag=sys.argv[5]
expect=int(sys.argv[6])
env=dict(os.environ); env['PAL_MAXPRINT']='10000000'
df=f"runs/pal_{tag}.done.jsonl"
if os.path.exists(df): os.remove(df)
r=subprocess.run([sys.executable,'pal_lists.py',str(v),str(K),str(q),json.dumps(gens),'5',tag],capture_output=True,text=True,env=env)
ps=PalShape(v,K,q,gens); n=ps.d['n']
def pack(s): 
    x=0
    for i,c in enumerate(s): x|=int(c)<<(2*i)
    return np.array([x],dtype=np.uint64)
def nu_all(arrs): return [swap_ex(ps,a) for a in arrs]
tot=0
for l in open(df):
    rec=json.loads(l)
    cnt=0
    for sol in rec.get('sols',[]):
        parts=sol.split()[1:]; X=[pack(p) for p in parts]
        w=1
        for kind,j in rec['halved']:
            if kind=='sigma':
                if sigma_ex(ps,X[j])[0]!=X[j][0]: w*=2
            else:
                Y=nu_all(X)
                if any(Y[t][0]!=X[t][0] for t in range(4)): w*=2
        cnt+=w
    tot+=cnt*rec['mult']
print('weighted total',tot,'expected',expect,'OK' if tot==expect else 'MISMATCH')
