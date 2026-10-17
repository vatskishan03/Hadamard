"""Sweep the undecided skew census values over cyclotomic products (e >= 2) with an extra multiplier group K_m on Z_m
(search/cyclo_product.py): stage 1 exactly; if it has classes and the lists are small, join every class (stop at the
first hit of each class, --sub 8000). One JSON line per structure on stdout; hits go to $HADAMARD_WORK/skew_hits.jsonl.

  python3 search/skew_kmult_sweep.py [MAXMEM [DONE_LOG]]      MAXMEM: members per block (default 2e7)

The run of 2026-10-13 (MAXMEM 2e7) found v = 1417 (skew order 5668), see NOTES section 16."""
import sys, json, math, time, io, contextlib
from sympy import isprime, divisors
from pathlib import Path as _P; sys.path.insert(0, str(_P(__file__).resolve().parent))
import cyclo_product as CP
L=[int(l) for l in open(_P(__file__).resolve().parent / "data" / "skew_census_undecided.txt") if l.strip() and not l.startswith('#')]
solved={209,253,259,275,285,287,309,341,343,345,353,377,391,403,451,465,469,513,601,613,703,763,889,903,933,1501}
R=[v for v in L if v not in solved]
maxmem=float(sys.argv[1]) if len(sys.argv)>1 else 2e7
done=set()
if len(sys.argv)>2:
    for l in open(sys.argv[2]):
        try: r=json.loads(l); done.add((r['v'],r['q'],r['e'],r['k']))
        except Exception: pass
def subgroup_gens(m):
    """one generator per cyclic subgroup of Z_m^* (order > 1) not containing -1"""
    seen=set(); out=[]
    for k in range(2,m):
        if math.gcd(k,m)!=1: continue
        K={1}; x=k
        while x!=1: K.add(x); x=x*k%m
        K=frozenset(K)
        if K in seen or (m-1) in K: continue
        seen.add(K); out.append(k)
    return out
jobs=[]
for v in R:
    for q in divisors(v):
        if q<3 or not isprime(q): continue
        m=v//q
        if m%q==0 or m<3: continue
        for e in divisors(q-1):
            if e<2 or e%2 or ((q-1)//e)%2==0: continue
            F=(q-1)//e
            if F<3 or 4*F-4>4*m*e or 32*(F-1)<4*m*e: continue
            for k in [1]+subgroup_gens(m):
                Km={1}; x=k
                while x!=1: Km.add(x); x=x*k%m
                norb=len({min(a*t%m for t in Km) for a in range(m)})
                mem=math.comb(e,e//2)**norb
                if norb>22 or mem>maxmem*50: continue
                if (v,q,e,k) in done: continue
                jobs.append((mem,v,q,e,k))
jobs.sort()
print(len(jobs),'jobs',flush=True)
for _,v,q,e,k in jobs:
    t=time.time()
    try:
        S=CP.Structure(v,q,e,kgen=k)
        # estimate per-block member count
        nkr=len(S.kreps); pat=math.comb(e,e//2)
        mem=pat**nkr
        if mem>maxmem*50: print(json.dumps(dict(v=v,q=q,e=e,k=k,skip='members',log10=round(math.log10(mem),1))),flush=True); continue
        quads=CP.stage1(S)
        rec=dict(v=v,q=q,e=e,k=k,K=len(S.K),classes=len(quads),sec=round(time.time()-t))
        if quads and mem<=maxmem:
            hits=0
            for c in range(len(quads)):
                buf=io.StringIO()
                with contextlib.redirect_stdout(buf):
                    found=CP.run(v,q,e,False,c,1,8000,-1,2,False,k)
                if found: hits+=1; rec['hit_config']=c; break
            rec['searched']=True; rec['hits']=hits
        print(json.dumps(rec),flush=True)
        if rec.get('hits'): print('HIT',v,flush=True)
    except MemoryError:
        import gc; gc.collect(); print('{"v": %d, "q": %d, "e": %d, "k": %d, "error": "MemoryError"}' % (v,q,e,k),flush=True)
    except Exception as ex:
        print(json.dumps(dict(v=v,q=q,e=e,k=k,error=repr(ex)[:150])),flush=True)
