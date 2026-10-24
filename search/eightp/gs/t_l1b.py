import numpy as np, time
from gsk import *
from gl1 import *
sh=GShape([5],3); l1=L1(sh)
print(sh.nu, sh.feasible_totals(), sh.Utarget, sh.repsPM)
t0=time.time()
lists={n:enum_weight(sh,n) for n in range(0,6)}
print({n:(None if v is None else len(v[0])) for n,v in lists.items()}, time.time()-t0)
al=(0,0,0,0)
for al in [(0,0,0,0),(1,1,1,1),(0,2,3,3)]:
    if lists[al[0]] is None: continue
    specs=[BlockSpec(al[b],S=lists[al[b]][0],E=lists[al[b]][1]) for b in range(4)]
    t0=time.time(); cnt,S,E,st=l1.run(specs); print(al,cnt,st,time.time()-t0,flush=True)
