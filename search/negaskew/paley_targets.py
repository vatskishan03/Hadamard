"""Paley plug-in shapes G = Z_K x F_q (q = 3 mod 4 prime >= 7, K odd) for the census-undecided v, with the minimum number of positions
at which the defect weights can add up to 4|K| (the defect identity (q+1) e+ + (q-3) e- = 4|K| of NOTES.md sections 88 and 89).
For a skew block A the quaternary class always has the forced half defect (q+1)/2 or (q-3)/2 at k = 0; the other positions weigh
2q+2, 2q-2, 2q-6 (a pair k, -k in A) or q+1, q-3 (one position of B, C, D).  The rigid shapes are those with N = 2 positions:
q = 4|K| - 1 or 4|K| + 3, i.e. v = 473, 517, 1139 and 1207; with one position only (N = 1) the shape is empty (413 = 7 * 59: only
A's own half defect fits).
usage: python3 paley_targets.py          the shapes with K <= 45 of the Table 5 entries m >= 4 (or undecided beyond 999) not in skew/
       python3 paley_targets.py 413 473 517 1139 1207       the rows of these v, whether or not skew/ settles them
Rows: (v, Table 5 entry m, q, K, N real class, N quaternary class including A's half defect); writes runs/paley_targets.json."""
import json,os,sys
from sympy import factorint
HERE=os.path.dirname(os.path.abspath(__file__))
T5={int(k):v for k,v in json.load(open(os.path.join(HERE,'t5.json'))).items()}
cen=[int(l) for l in open(os.path.join(HERE,'..','data','skew_census_undecided.txt')) if l.strip() and not l.startswith('#')]
found=set(int(d)//4 for d in os.listdir(os.path.join(HERE,'..','..','skew')) if d.isdigit()); found.add(919)
args=[int(a) for a in sys.argv[1:]]
und=args if args else [v for v in cen if v not in found]
def pp(q): return len(factorint(q))==1
def mincount(weights,total,cap):
    INF=10**9; dp=[INF]*(total+1); dp[0]=0
    for s in range(1,total+1):
        for w in weights:
            if w>0 and s>=w and dp[s-w]+1<dp[s]: dp[s]=dp[s-w]+1
    return dp[total] if total>=0 else INF
rows=[]
for v in und:
    m5=T5.get(v,2) if v<=999 else None
    for q in range(7,v):
        if v%q or not pp(q) or q%4!=3: continue
        K=v//q
        if K%2==0: continue
        # real
        rr=mincount([q+1,q-3],4*K,4*K)
        # quaternary doubled weights: others (any block, any position incl A pairs as 2w): single positions
        W=[2*q+2,2*q-2,2*q-6,q+1,q-3]
        best=10**9
        for wa in (q+1,q-3):
            rem=8*K-wa
            if rem<0: continue
            c=mincount(W,rem,4*K)
            best=min(best,c+1)
        rows.append((v,m5,q,K,rr if rr<10**8 else None,best if best<10**8 else None))
os.makedirs('runs',exist_ok=True)
json.dump(rows,open('runs/paley_targets.json','w'))
sel=rows if args else [r for r in rows if r[5] is not None and r[3]<=45 and (r[1] is None or r[1]>=4)]
sel.sort(key=lambda r:(r[3],r[2]))
print('v m5 q K realN quatN(incl A0 forced half-defect)')
for r in sel: print(r)
