"""Run the exact compression-chain search (chain.py) on the survivors of census_screen.py, best model estimate first.

  python3 search/screen_chain.py [BUDGET] [MAXO] [MAXBRANCHES] [--maxfibre N] [--skip v,v,...] [--order model|orbits]
      -> $HADAMARD_WORK/screen_chain[bK|x].jsonl

MAXBRANCHES limits each structure to its heaviest level-1 branches (breadth-first pass, file tag bK); 0 searches
every branch (exhaustive pass, file tag x). The model estimate only sets the order (for composite v it is
unreliable); each structure gets the full chain search with the given time budget. `--order orbits` takes the
structures with the fewest orbits first instead. A v is skipped once it has a verified hit."""
import argparse, json, os, sys, time
from pathlib import Path as _P; sys.path[:0] = [str(_P(__file__).resolve().parent)]
import workdir as W
import chain

ap = argparse.ArgumentParser()
ap.add_argument('budget', type=float, nargs='?', default=1800)
ap.add_argument('maxo', type=int, nargs='?', default=46)
ap.add_argument('maxbranches', type=int, nargs='?', default=None)
ap.add_argument('--maxfibre', type=float, default=3e6)
ap.add_argument('--skip', default='')
ap.add_argument('--order', choices=('model', 'orbits'), default='model')
a = ap.parse_args()
maxbranches = a.maxbranches or None
tag = f'b{a.maxbranches}' if a.maxbranches else ('x' if a.maxbranches == 0 else '')
out = W.WORK / f'screen_chain{tag}.jsonl'
skip = {int(x) for x in a.skip.split(',') if x}
done = set()
if out.exists():
    for l in open(out):
        r = json.loads(l); done.add((r['v'], tuple(r['gens']), tuple(r['pattern'])))
solved = {json.loads(l)['n'] for l in open(W.hits())} if os.path.exists(W.hits()) else set()
while True:
    recs = [json.loads(l) for l in open(W.WORK / 'census_screen.jsonl')]
    todo = [r for r in recs if r['alive'] and r['O'] <= a.maxo and r['v'] not in solved and r['v'] not in skip
            and (r['v'], tuple(r['gens']), tuple(r['pattern'][1:] if len(r['pattern']) == 4 else r['pattern'])) not in done]
    if not todo: break
    r = max(todo, key=lambda r: r['log10E_gen'] if a.order == 'model' else (-r['O'], r['log10E_gen']))
    pat = tuple(r['pattern'])
    t = time.time()
    try:
        res = chain.run(r['v'], r['gens'], pat, exhaustive=False, verbose=True, budget=a.budget, maxbranches=maxbranches,
                        join_timeout=900 if maxbranches else None, maxfibre=int(a.maxfibre))
    except Exception as e:
        res = dict(status=f'error: {e!r}')
    rec = dict(v=r['v'], gens=r['gens'], H=r['H'], O=r['O'], pattern=list(pat), log10E_gen=r['log10E_gen'],
               status=res.get('status'), branches=res.get('branches'), truncated=res.get('truncated'),
               solutions=res.get('solutions'), seconds=round(time.time() - t))
    with open(out, 'a') as f: f.write(json.dumps(rec) + '\n')
    print(json.dumps(rec), flush=True)
    done.add((r['v'], tuple(r['gens']), pat))
    if res.get('found'): solved.add(r['v'])
