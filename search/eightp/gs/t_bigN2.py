# scan large cyclic K (65..130) with big H (few orbits) for NONZERO exhaustive counts and compare with the G-level brute force
import sys,time,itertools; sys.path.insert(0,'.')
from gsk import *; from gs_exh import exhaust; from gbf import count_G
from pk import signed_subgroups
found=0
for mode in ('spence','gs'):
  for N in [66,70,72,78,84,90,102,110,120,65,69,75,77,85,91,95,99,105,117]:
    if mode=='spence' and N%2: continue
    for r in (3,5,7):
        try: subs=signed_subgroups([N],maxgens=2)
        except Exception: continue
        for (G_,gens) in subs:
            hg=[(tuple(u),int(s)) for (u,s) in gens]
            try: sh=GShape([N],r,hg,mode,None)
            except AssertionError: continue
            if sh.nu>5 or sh.nu<3: continue
            t0=time.time()
            try: gc,info=count_G([N],r,hg,mode,None)
            except Exception: continue
            if gc==0: continue
            tot,sols=exhaust(sh,tmax=100,max_inst=500000)
            print(mode,[N],r,hg,'nu',sh.nu,'G-level',gc,'K-level',tot['E_G'],'exhaustive',tot['exhaustive'],'MATCH' if (tot['exhaustive'] and gc==tot['E_G']) else '??',f'{time.time()-t0:.0f}s',flush=True)
            found+=1
            if found>=6: sys.exit(0)
