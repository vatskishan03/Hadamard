# prime-power r (GF(9)) : exhaustive K-level engine vs independent G-level brute force
import sys,time; sys.path.insert(0,'.')
from gsk import *; from gs_exh import exhaust; from gbf import count_G
for (mode,K,r,hg,T) in [('spence',[2],9,[],None),('gs',[2],9,[],None),('gs',[3],9,[],None),('bord',[2],9,[],[(0,),(1,)]),('spence',[2],9,[([1],-1)],None),('gs',[5],9,[([4],1)],None),('gs',[3],25,[],None)]:
    t0=time.time()
    try: gc,info=count_G(K,r,[(u,p) for u,p in hg],mode,T)
    except Exception as ex: print(mode,K,r,'G-level skipped',repr(ex)[:60],flush=True); continue
    sh=GShape(K,r,[(u,p) for u,p in hg],mode,T)
    tot,sols=exhaust(sh,tmax=100,max_inst=800000)
    print(mode,K,r,hg,T,'G-level',gc,'K-level',tot['E_G'],'exhaustive',tot['exhaustive'],'MATCH' if (tot['exhaustive'] and gc==tot['E_G']) else '??',f'{time.time()-t0:.0f}s',flush=True)
