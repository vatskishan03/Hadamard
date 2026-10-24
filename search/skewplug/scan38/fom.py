"""figures of merit of slot files: the sum of 10^logE over the 40 best slots, the E-weighted mean of 10^logrho.  usage (in search/skewplug): python3 scan38/fom.py est/NAME.ext.jsonl ..."""
import json,sys,math
for nm in sys.argv[1:]:
    rows=[json.loads(l) for l in open(nm)]
    rows=[r for r in rows if 'logE' in r and r['logE']>-50]
    rows.sort(key=lambda r:-r['logE'])
    top=rows[:40]
    E=[10**r['logE'] for r in top]; rho=[10**r['logrho'] for r in top]
    fom=sum(e*r for e,r in zip(E,rho))/sum(E)
    print(nm,len(rows),'top40 sumE %.3g'%sum(E),'E-weighted rho 10^%.2f'%math.log10(fom),'best',top[0]['logE'],top[0]['logrho'],'logN1 range',min(r['logN1'] for r in top),max(r['logN1'] for r in top))
    # best rho slots with logE>-3
    rr=sorted([r for r in rows if r['logE']>-2],key=lambda r:-r['logrho'])[:3]
    print('  best rho w/ logE>-2:',[(r['al'],r['logE'],r['logrho']) for r in rr])
