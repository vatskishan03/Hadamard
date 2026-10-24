"""EXACT count of Lock-consistent local states per +-pair of orbits for an origin-defect pattern.
 block A (skew): u(-O) = i conj u(O); free block: u(O), u(-O) arbitrary units.  Phi_X(k) = e_X conj u_X(k) + conj(e_X) u_X(-k) ;  Lock_z : sum_X Phi_X(k) = -z .
"""
import itertools, sys
U=[1,1j,-1,-1j]
def states(eps_list, kinds, z):
    """eps_list: origin e of each defect block; kinds: 'A' or 'F' (free). returns list of states"""
    opts=[]
    for kd in kinds:
        if kd=='A': opts.append([(a,1j*a.conjugate()) for a in U])
        else: opts.append([(a,b) for a in U for b in U])
    res=[]
    for st in itertools.product(*opts):
        tot=0
        for e,(a,b) in zip(eps_list,st): tot+=e*a.conjugate()+e.conjugate()*b
        if abs(tot+z)<1e-9: res.append(st)
    return res
def run(name, eps_choices, kinds, zs=(0,1)):
    out={}
    for z in zs:
        n=0
        for eps in itertools.product(*eps_choices): n+=len(states(list(eps),kinds,z))
        out[z]=(n, sum(1 for _ in itertools.product(*eps_choices)))
    return out
h=0.5
Meps=[(1+1j)*h,(1-1j)*h,(-1+1j)*h,(-1-1j)*h]   # all four half-defect e values
Feps=[1,-1,1j,-1j]
# For counting we fix the origin data of each block (e, cost class) and report the count for each eps
def report(name, spec, zs=(0,1)):
    """spec: list of (kind, list of e choices) for the defect blocks, kind in 'A','F' """
    kinds=[k for k,_ in spec]; ch=[c for _,c in spec]
    print(name)
    for z in zs:
        for eps in itertools.product(*ch):
            n=len(states(list(eps),kinds,z))
            if n: print('   z=%d eps=%s : %d local states'%(z,[complex(round(e.real,2),round(e.imag,2)) for e in eps],n))
# patterns (free block origin x0=1; costs fix e): F cost q-3: e=-1 ; F cost q-1: e=+-i ; F cost q+1: e=1 ; M cost (q+1)/2: e in {(1+i)/2,(1-i)/2} ; M cost (q-3)/2: e in {(-1+i)/2,(-1-i)/2}
# block A: x0 in {1,-i}: cost (q+1)/2 => x0=1,e=(1+i)/2 or x0=-i, e = (1-i)/2*? : Re(x0 conj e)=1/2: x0=1:e=(1+-i)/2 ; skew requires e = x2/x1-related: e = (x1+x2)/2 with x2 = i conj x1.
