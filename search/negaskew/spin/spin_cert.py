"""Certificate of a spin solution: expands a solution (A, B) of qspin / qspin_fast to the four orbit-exponent strings
(C = B o mu, D = B o mu^2) and calls ../certify.py (G-level check, real matrix, out/cert_<8v>.json in ..).
usage: python3 spin_cert.py v "[gens]" mu Astring Bstring"""
import sys,json,subprocess,os
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..'))
from spin_drv import *
v=int(sys.argv[1]); gens=json.loads(sys.argv[2]); mu=int(sys.argv[3]); A=sys.argv[4]; B=sys.argv[5]
d=shape(v,gens); orb=d['orb']; reps=d['reps']; n=d['n']
op=[orb[(mu*reps[a])%v] for a in range(n)]
e=[int(c) for c in B]
C=[e[op[a]] for a in range(n)]; D=[e[op[op[a]]] for a in range(n)]
strs=[A,B,"".join(map(str,C)),"".join(map(str,D))]
print(" ".join(strs))
subprocess.run([sys.executable,os.path.join(HERE,'..','certify.py'),str(v),json.dumps(gens)," ".join(strs)],cwd=os.path.join(HERE,'..'))
