import sys; sys.path.insert(0,'.')
from gs_samp import *
def go(mode,K,r,hg,T,allocs,LA,LC,secs=60,maxrounds=8):
    sh=make_shape(mode,K,r,hg,T)
    al=parse_allocs(allocs,LA,LC)
    print('==',mode,K,r,hg,T,'feasible E',sh.feasible_totals(),flush=True)
    res=run(sh,al,secs,777,f'runs/t2_{mode}',log=lambda *a: print(*a,flush=True),maxrounds=maxrounds)
    print('result',res[0],res[1],flush=True)
go('spence',[4],3,[],None,'2:2:2:2/2:0:0:0,3:2:2:2/2:0:0:0,2:2:2:3/2:0:0:0',300,3000)
go('spence',[2,2],5,[],None,'2:2:2:2/2:0:0:0,2:2:2:3/2:0:0:0',300,3000)
