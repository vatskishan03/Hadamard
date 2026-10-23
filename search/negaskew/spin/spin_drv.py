"""Driver of the quaternary SPIN structure (C = B o mu, D = B o mu^2 with mu^3 in H, A twisted-skew, B free) over Z_v with the
multiplier group H = <gens> (odd order): writes the shape file for qspin / qspin_fast, runs the engine, and checks a solution
at the group level and as a real skew matrix (v <= 60).
usage: python3 spin_drv.py v "[gens]" mu      e.g.  python3 spin_drv.py 7 "[]" 2     (COUNT 1820, order 56)
                                                    python3 spin_drv.py 13 "[]" 3    (COUNT 110448, order 104)
Needs ./qspin_fast (cc -O3 -march=haswell -o qspin_fast qspin_fast.c); SPIN_ENGINE=./qspin selects the plain version.
The shape file is shapes/<tag>_<v>_<mu>.txt: v n / sizes / neg / orbit permutation of mu (twice) / the tensor N."""
import numpy as np, subprocess, sys, json, os
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..'))
from qdrv import *
def spin_shape(v,gens,mu,fname):
    d=shape(v,gens)
    orb=d['orb']; reps=d['reps']; n=d['n']
    op=[orb[(mu*reps[a])%v] for a in range(n)]
    assert sorted(op)==list(range(n)) and op[0]==0
    with open(fname,'w') as f:
        f.write(f"{v} {n}\n"+" ".join(map(str,d['sz']))+"\n"+" ".join(map(str,d['neg']))+"\n"+" ".join(map(str,op))+"\n"+" ".join(map(str,op))+"\n")
        N=d['N']
        for c in range(1,n):
            for a in range(n): f.write(" ".join(map(str,N[c,a]))+"\n")
    return d,op
def run_spin(v,gens,mu,maxprint=3,timeout=None,tag='sp'):
    os.makedirs('shapes',exist_ok=True)
    fn=f"shapes/{tag}_{v}_{mu}.txt"; d,op=spin_shape(v,gens,mu,fn)
    r=subprocess.run([os.environ.get('SPIN_ENGINE','./qspin_fast'),fn,str(maxprint)],capture_output=True,text=True,timeout=timeout)
    return d,op,r.stdout,r.stderr
def g_blocks(d,op,v,mu,A,B):
    ea=expand(d,A); eb=expand(d,B)
    ec=[eb[(mu*g)%v] for g in range(v)]; ed=[eb[(mu*mu*g)%v] for g in range(v)]
    return [ea,eb,ec,ed]
def check_G(v,E):
    tot=0
    for e in E:
        u=np.array([1j**x for x in e]); tot=tot+np.array([np.sum(u[(np.arange(v)+s)%v]*np.conj(u)) for s in range(v)])
    return tot[0],np.max(np.abs(tot[1:]))
if __name__=="__main__":
    v=int(sys.argv[1]); gens=json.loads(sys.argv[2]); mu=int(sys.argv[3])
    d,op,out,err=run_spin(v,gens,mu)
    print(err.strip()); print(out.strip())
    for line in out.split('\n'):
        if line.startswith('SOL'):
            A=[int(c) for c in line.split()[1][2:]]; B=[int(c) for c in line.split()[2][2:]]
            E=g_blocks(d,op,v,mu,A,B); print('G-level',check_G(v,E))
            if v<=60:
                H=real_matrix(v,E); print('real order',H.shape[0],check_real(H))
