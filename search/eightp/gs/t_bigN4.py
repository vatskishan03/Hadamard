# positive control at |K| = 283 (the n = 1981 census shape's K and H): exhaustive K-level count for r = 3, 5, 9, 11, 13 (solutions are checked exactly by check_klevel inside exhaust)
import sys,time; sys.path.insert(0,'.')
from gsk import *; from gs_exh import exhaust; from gl2io import channels
for r in (3,5,9,11,13,7):
    t0=time.time()
    sh=GShape([283],r,[([4],1)],'gs')
    tot,sols=exhaust(sh,tmax=150,max_inst=500000,want_sols=5)
    print('gs K=[283] H=<4> r=',r,'nu',sh.nu,'channels',len(channels(sh)),'allocs',tot['allocs'],'n1',tot['n1'],'E_G',tot['E_G'],'exhaustive',tot['exhaustive'],'skipped',len(tot['skipped']),f'{time.time()-t0:.0f}s',flush=True)
