# POSITIVE CONTROL at large |K|: K = Z_p (p prime, 65 < p < 330), H = squares of Z_p* (3 orbits), r in 3..13, mode gs: G-level count vs K-level exhaustive; looks for NONZERO counts
import sys,time; sys.path.insert(0,'.')
from gsk import *; from gs_exh import exhaust; from gbf import count_G
import sympy
found=0
for p in sympy.primerange(67,330):
    g=sympy.primitive_root(p)
    for r in (3,5,7,9,11,13):
        for tw in (1,):
            hg=[([g*g%p],1)]
            t0=time.time()
            try: sh=GShape([p],r,hg,'gs')
            except AssertionError: continue
            try: gc,info=count_G([p],r,[(u,ph) for u,ph in hg],'gs',None)
            except Exception as ex: continue
            if gc==0: continue
            tot,sols=exhaust(sh,tmax=120,max_inst=800000)
            print('gs',[p],r,hg,'nu',sh.nu,'G-level',gc,'K-level',tot['E_G'],'exhaustive',tot['exhaustive'],'MATCH' if (tot['exhaustive'] and gc==tot['E_G']) else '??',f'{time.time()-t0:.0f}s',flush=True)
            found+=1
            if found>=6: sys.exit(0)
print('scan finished; nonzero found:',found,flush=True)
