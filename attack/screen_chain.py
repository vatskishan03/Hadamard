"""Run the exact compression-chain search (chain.py) on the survivors of census_screen.py, best model estimate first.

  python3 attack/screen_chain.py [BUDGET] [MAXO] [MAXBRANCHES]   -> $HADAMARD_WORK/screen_chain[bK].jsonl

The model estimate only sets the order (for composite v it is unreliable); each structure gets the full chain
search with the given time budget. A v is skipped once it has a verified hit."""
import json, sys, time
from pathlib import Path as _P; sys.path[:0] = [str(_P(__file__).resolve().parent)]
import workdir as W
import chain

budget = float(sys.argv[1]) if len(sys.argv) > 1 else 1800
maxo = int(sys.argv[2]) if len(sys.argv) > 2 else 46
maxbranches = int(sys.argv[3]) if len(sys.argv) > 3 else None     # breadth-first pass: heaviest branches only
tag = f'b{maxbranches}' if maxbranches else ''
out = W.WORK / f'screen_chain{tag}.jsonl'
done = set()
if out.exists():
    for l in open(out):
        r = json.loads(l); done.add((r['v'], tuple(r['gens']), tuple(r['pattern'])))
import os
solved = {json.loads(l)['n'] for l in open(W.hits())} if os.path.exists(W.hits()) else set()
while True:
    recs = [json.loads(l) for l in open(W.WORK / 'census_screen.jsonl')]
    todo = [r for r in recs if r['alive'] and r['O'] <= maxo and r['v'] not in solved
            and (r['v'], tuple(r['gens']), tuple(r['pattern'][1:] if len(r['pattern']) == 4 else r['pattern'])) not in done]
    if not todo: break
    r = max(todo, key=lambda r: r['log10E_gen'])
    pat = tuple(r['pattern'])
    t = time.time()
    try:
        res = chain.run(r['v'], r['gens'], pat, exhaustive=False, verbose=True, budget=budget, maxbranches=maxbranches)
    except Exception as e:
        res = dict(status=f'error: {e!r}')
    rec = dict(v=r['v'], gens=r['gens'], H=r['H'], O=r['O'], pattern=list(pat), log10E_gen=r['log10E_gen'],
               status=res.get('status'), branches=res.get('branches'), truncated=res.get('truncated'), seconds=round(time.time() - t))
    with open(out, 'a') as f: f.write(json.dumps(rec) + '\n')
    print(json.dumps(rec), flush=True)
    done.add((r['v'], tuple(r['gens']), pat))
    if res.get('found'): solved.add(r['v'])
