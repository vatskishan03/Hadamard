import sys; sys.path.insert(0,'.')
from gs_samp import *
def go(mode,K,r,hg,T,allocs,LA,LC,secs=60):
    sh=make_shape(mode,K,r,hg,T)
    al=parse_allocs(allocs,LA,LC)
    print('==',mode,K,r,hg,T,'feasible E',sh.feasible_totals(),flush=True)
    res=run(sh,al,secs,12345,f'runs/t_{mode}',log=lambda *a: print(*a,flush=True),maxrounds=8)
    print('result',res[0],res[1],flush=True)
go('gs',[7],3,[[[2],1]],None,'1:3:3:3,3:3:3:3,4:3:3:3',300,3000)
go('spence',[4],3,[],None,'3.0:3.0:3.0:3.0/2:0:0:0,4.0:3.0:3.0:3.0/2:0:0:0,2.0:2.0:2.0:2.0/2:0:0:0',300,3000)
go('bord',[6],5,[],[[0],[2],[4]],'1.0:1.0:1.0:1.0/2:2:2:0,2.0:2.0:1.0:1.0/2:2:2:0,0:2:2:2/2:2:2:0',300,3000)
