"""G-level count of the normalised quaternary skew quadruples of a plug-in shape Z_m x F_q with H = H_K x squares, by qskew over Z_v (the
reference of qp.c's counts, NOTES.md section 89).
usage: python3 g_count.py m q "[H_K generators]" [stride]       e.g.  python3 g_count.py 7 11 "[2]"      -> COUNT 3840
Writes shapes/g_<v>_<gens>.txt, runs ./qskew with mode 1 (quaternary), stride (default 1) = full count; a stride k counts every k-th
left tuple only.  Compare with  python3 qshape.py "[m]" q "[H_K]" shapes/p.txt && ./qp shapes/p.txt maxprint=0  (the count in the STATS line)."""
import sys, json, subprocess, os
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from qshape import cyc_gens
import qdrv
os.makedirs('shapes',exist_ok=True)
m=int(sys.argv[1]); q=int(sys.argv[2]); hK=json.loads(sys.argv[3]); stride=sys.argv[4] if len(sys.argv)>4 else '1'
gens=cyc_gens(m,q,hK); v=m*q
d=qdrv.shape(v,gens,f"shapes/g_{v}_{'_'.join(map(str,gens))}.txt")
print('v',v,'gens',gens,'n orbits',d['n'], flush=True)
r=subprocess.run(['./qskew',f"shapes/g_{v}_{'_'.join(map(str,gens))}.txt",'1','0','0',stride],capture_output=True,text=True)
print(r.stderr.strip()); print(r.stdout.strip())
