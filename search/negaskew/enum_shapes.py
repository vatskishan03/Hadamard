import json,os,sys
from math import gcd
from sympy import totient
T5={int(k):v for k,v in json.load(open('t5.json')).items()}
cen=[int(l) for l in open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','data','skew_census_undecided.txt')) if l.strip() and not l.startswith('#')]
found=set(int(d)//4 for d in os.listdir(os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','..','skew')) if d.isdigit()); found.add(919)
und=[v for v in cen if v not in found]
def order(h,v):
    k=1;x=h
    while x!=1: x=x*h%v;k+=1
    return k
def gen_sub(v,gens):
    H={1};fr=[1]
    while fr:
        x=fr.pop()
        for g in gens:
            y=x*g%v
            if y not in H: H.add(y);fr.append(y)
    return frozenset(H)
def norb(v,H):
    seen=bytearray(v);c=0
    for g in range(v):
        if not seen[g]:
            c+=1
            for h in H: seen[g*h%v]=1
    return c
def shapes(v,nmax):
    units=[h for h in range(1,v) if gcd(h,v)==1]
    cyc={}
    for h in units:
        o=order(h,v)
        if o%2==1 and o>1:
            S=gen_sub(v,[h]); cyc[S]=h
    subs={S:[h] for S,h in cyc.items()}
    # pairwise joins
    items=list(cyc.items())
    for i in range(len(items)):
        for j in range(i+1,len(items)):
            S=gen_sub(v,[items[i][1],items[j][1]])
            if S not in subs and all(x%2==1 or True for x in [1]) and len(S)%2==1: subs[S]=[items[i][1],items[j][1]]
    out=[]
    for S,g in subs.items():
        n=norb(v,S)
        if n<=nmax: out.append((n,len(S),g))
    out.sort()
    return out
if __name__=="__main__":
    nmax=int(sys.argv[1]); vmax=int(sys.argv[2])
    res={}
    for v in und:
        if v>vmax: continue
        m=T5.get(v,2) if v<=999 else None
        if v<=999 and m<4: continue
        sh=shapes(v,nmax)
        if sh: res[v]=(m,sh)
    json.dump({str(k):v for k,v in res.items()},open(f'enum_n{nmax}_v{vmax}.json','w'))
    cnt=sum(len(x[1]) for x in res.values())
    print(len(res),'v with shapes',cnt,'shapes')
    for v,(m,sh) in sorted(res.items()): print(v,m,[(n,o) for n,o,g in sh])
