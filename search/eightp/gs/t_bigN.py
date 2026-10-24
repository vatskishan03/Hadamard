# large |K| (> 64, > 126 channels truncated) : exhaustive engine vs G-level brute force on shapes with very few H-orbits
import sys,time; sys.path.insert(0,'.')
from gsk import *; from gs_exh import exhaust; from gbf import count_G
for (mode,K,r,hg,T) in [('gs',[127],3,[([9],1)],None),('gs',[127],3,[([3],-1)],None),('gs',[127],5,[([9],1)],None),('gs',[127],7,[([9],1)],None),('gs',[73],3,[([2],1)],None),('gs',[85],3,[([4],1)],None),('gs',[129],3,[([4],1)],None)]:
    t0=time.time()
    try: gc,info=count_G(K,r,[(u,p) for u,p in hg],mode,T)
    except Exception as ex: print(mode,K,r,hg,'G-level skipped',repr(ex)[:80],flush=True); continue
    sh=GShape(K,r,[(u,p) for u,p in hg],mode,T)
    tot,sols=exhaust(sh,tmax=150,max_inst=800000)
    print(mode,K,r,hg,'nu',sh.nu,'G-level',gc,'K-level',tot['E_G'],'exhaustive',tot['exhaustive'],'MATCH' if (tot['exhaustive'] and gc==tot['E_G']) else '??',f'{time.time()-t0:.0f}s',flush=True)
