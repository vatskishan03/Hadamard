"""CP-SAT model (OR-Tools) of the quaternary skew orbit-level system of README.md: one 4-valued variable per orbit and
block, A twisted-skew, B, C, D normalised, the PAF equations over all classes as linear constraints on indicator variables.
usage: python3 cpsat.py v "[gens]" seconds [workers]       e.g.  python3 cpsat.py 191 "[180]" 240
Prints the solver status and, if feasible, a SOL line (four orbit-exponent strings, as qskew).  It does not find the
known solution of v = 191 (n = 11 orbits) within 240 s on one worker (NOTES.md section 88)."""
import sys,json,time
from ortools.sat.python import cp_model
from qdrv import shape
def solve(v,gens,tlim,workers=1,seed=1):
    d=shape(v,gens); n=d['n']; N=d['N']; neg=d['neg']
    m=cp_model.CpModel()
    y=[[[m.NewBoolVar(f"y{X}_{a}_{k}") for k in range(4)] for a in range(n)] for X in range(4)]
    for X in range(4):
        for a in range(n): m.AddExactlyOne(y[X][a])
        # normalisation / skew structure
    for X in range(1,4):
        m.Add(y[X][0][0]==1)
    # A: u0 in {1,-i}: exps {0,3}; pair relation e'=1-e mod 4
    m.Add(y[0][0][1]+y[0][0][2]==0)
    for a in range(1,n):
        b=neg[a]
        if a<b:
            for k in range(4): m.Add(y[0][b][(1-k)%4]==y[0][a][k])
    # delta indicators
    tre=[[0]*0 for _ in range(n)]
    expr_re={c:[] for c in range(1,n)}; expr_im={c:[] for c in range(1,n)}; const_re={c:0 for c in range(1,n)}
    for X in range(4):
        for a in range(n):
            for b in range(a+1,n):
                # only if some class uses pair (a,b)
                need=any(N[c,a,b]!=0 or N[c,b,a]!=0 for c in range(1,n))
                if not need: continue
                w=[m.NewBoolVar(f"w{X}_{a}_{b}_{dl}") for dl in range(4)]
                m.AddExactlyOne(w)
                for k in range(4):
                    for l in range(4):
                        m.AddBoolOr([y[X][a][k].Not(),y[X][b][l].Not(),w[(k-l)%4]])
                for c in range(1,n):
                    s=int(N[c,a,b]+N[c,b,a]); t_=int(N[c,a,b]-N[c,b,a])
                    if s: expr_re[c]+= [s*w[0],-s*w[2]]
                    if t_: expr_im[c]+= [t_*w[1],-t_*w[3]]
    # diagonal terms are constants N[c][a][a]*1 (real)
    for c in range(1,n):
        const=0
        for a in range(n): const+=4*int(N[c,a,a])
        m.Add(sum(expr_re[c])+const==0)
        if expr_im[c]: m.Add(sum(expr_im[c])==0)
    s=cp_model.CpSolver(); s.parameters.max_time_in_seconds=tlim; s.parameters.num_workers=workers; s.parameters.random_seed=seed
    t0=time.time(); r=s.Solve(m)
    print('status',s.StatusName(r),round(time.time()-t0,1),'s')
    if r in (cp_model.OPTIMAL,cp_model.FEASIBLE):
        E=[[ [k for k in range(4) if s.Value(y[X][a][k])][0] for a in range(n)] for X in range(4)]
        return d,E
    return d,None
if __name__=="__main__":
    v=int(sys.argv[1]); gens=json.loads(sys.argv[2]); tl=float(sys.argv[3]); w=int(sys.argv[4]) if len(sys.argv)>4 else 1
    d,E=solve(v,gens,tl,w)
    if E:
        print("SOL",*["".join(map(str,e)) for e in E])
