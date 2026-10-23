"""Driver of oskew: writes the shape file of (v, H = <gens>) (H of odd order, or at least -1 not in H) and runs the orbit-level exhaustive join of the octic
skew class (16v, README.md).  usage: python3 odrv.py v "[gens]"          e.g.  python3 odrv.py 61 "[12]"   -> normalised count 15360 (H of order 15, 5 orbits)
Needs ./oskew (cc -O2 -o oskew oskew.c).  The count is normalised by the 8 units w^k of each of B, C, D (multiply by 8^3 for all ordered)."""
import numpy as np, subprocess, sys, json, os, time
from dcfeas import orbit_data
def oshape(v,gens,fname):
    H,orb,reps,sz,neg=orbit_data(v,gens); n=len(reps)
    N=np.zeros((n,n,n),dtype=np.int64)
    for c in range(1,n):
        s=reps[c]
        for g in range(v): N[c,orb[g],orb[(g+s)%v]]+=1
    with open(fname,'w') as f:
        f.write(f"{v} {n}\n"+" ".join(map(str,sz))+"\n"+" ".join(map(str,neg))+"\n")
        for c in range(1,n):
            for a in range(n): f.write(" ".join(map(str,N[c,a]))+"\n")
    return dict(v=v,n=n,H=H,orb=orb,reps=reps,sz=sz,neg=neg)
def run(v,gens,maxprint=3,firstonly=0,tag='t',timeout=None):
    os.makedirs('shapes',exist_ok=True); fn=f"shapes/{tag}_{v}_{'_'.join(map(str,gens))}.txt"
    d=oshape(v,gens,fn); t0=time.time()
    r=subprocess.run(['./oskew',fn,str(maxprint),str(firstonly)],capture_output=True,text=True,timeout=timeout)
    return d,r.stdout,r.stderr,time.time()-t0
if __name__=="__main__":
    v=int(sys.argv[1]); gens=json.loads(sys.argv[2]) if len(sys.argv)>2 else [1]
    d,out,err,t=run(v,gens,maxprint=3)
    print(err.strip()); print(out.strip()); print('time %.1f s'%t)
