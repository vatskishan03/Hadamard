"""The published solution of order 2216 = 8 * 277 (certificate of skew/2216, orbit exponents in the qskew order) in the encoding of tri2: converts it to coset order, reduces it by the shifts and the tilde exactly as tri2 with red=1 does (A: least encoding over the e shifts; free blocks: least of block and tilde, sorted) and writes the candidate (A, B big, {C, D}) to runs/cand277.txt.  With a FILE (a tri2 run with red=1 at v = 277, e.g. data/c277_red.out) it also prints the matching SOL line (the check of search/NOTES.md section 92 that the published solution is among the 646 reduced solutions).
usage: python3 find277.py [FILE]"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from trilib import *
v,k=277,23
d=setup(v,k); e=d['e']; h=e//2; cos=d['cos']
orbj=[]; seen=set()
for x in range(1,v):
    j=cos[x]
    if j not in seen: seen.add(j); orbj.append(j)
# quaternary_orbit_exponents of skew/2216/certificate.json
cert="3313220001321 0301322112300 0023212321000 0033103132201".split()
blocks=[]
for p in cert:
    ex=[int(c) for c in p]; c=[0]*e
    for o in range(1,e+1): c[orbj[o-1]]=ex[o]
    blocks.append((ex[0],tuple(c)))
P=sum(paf(d,b) for b in blocks); print('PAF sum max |.| off origin',np.abs(P[1:]).max(),' P[0]',P[0])
print('row sums |R|^2',[int(round(abs(spectrum(d,b)[0])**2)) for b in blocks], 'S',[sum(I4[c] for c in b[1]) for b in blocks])
def shift(b,a): return (b[0],tuple(b[1][(j+a)%e] for j in range(e)))
def tilde(b): return (b[0], tuple((-b[1][(j+h)%e])%4 for j in range(e)))
def enc(b):  # tri2 encoding: ex = sum x[o]<<2o, o=0..e ; x[0]=u0, x[1+j]=c_j
    ex=b[0]
    for j in range(e): ex|=b[1][j]<<(2*(1+j))
    return ex
A=blocks[0]; fr=blocks[1:]
big=[b for b in fr if abs(sum(I4[c] for c in b[1]))>1e-9][0]; sm=[b for b in fr if b is not big]
if sum(I4[c] for c in big[1])!=(1+1j): big=tilde(big)
repA=min(((shift(A,a)) for a in range(e)),key=enc)
# shift of A by a is accompanied by the same shift of all blocks
best=None
for a in range(e):
    Ab=shift(A,a)
    if enc(Ab)!=enc(repA): continue
    Bb=shift(big,a); 
    Cs=[min((shift(x,a) for x in (c,tilde(c))),key=enc) for c in sm]
    cand=(enc(Ab),enc(Bb),tuple(sorted(enc(c) for c in Cs)))
    print('candidate (shift a=%d)'%a,cand)
    best=cand
HERE=os.path.dirname(os.path.abspath(__file__)); os.makedirs(os.path.join(HERE,'runs'),exist_ok=True)
open(os.path.join(HERE,'runs','cand277.txt'),'w').write(repr(best))
print('candidate',best)
if len(sys.argv)>1:     # new in the port: look the candidate up in a reduced tri2 run (SOL lines: A B {C, D} in tri2's digit strings)
    dig=lambda x:''.join(str((x>>(2*a))&3) for a in range(e+1))
    cA,cB,cCD=dig(best[0]),dig(best[1]),set(dig(x) for x in best[2]); found=0
    for l in open(sys.argv[1]):
        if l.startswith('SOL'):
            p=l.split()
            if p[4]==cA and p[5]==cB and set(p[6:8])==cCD: found+=1; print('found in',sys.argv[1],':',l.strip())
    print('matches in the log:',found)

