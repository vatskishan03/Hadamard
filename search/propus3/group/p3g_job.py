"""Resumable exhaustive job: one group structure (v, H given by ambient generators), all feasible row patterns, p3 engine.
  python3 p3g_job.py --id ID --v 343 --hgens "2,1" [--rows a,b,d ...] [--threads T] [--ns 3] [--shards N] [--workdir DIR] [--countlists]
hgens: ';'-separated generators, each = exponents of S_i (Singer / unit, one per prime of v) then of F_i (Frobenius).
Writes DIR/HIT_<id>.json ONLY for a hit that passed the exact group-PAF check and the matrix check H = H^T, H H^T = 4 v I.
Completed (pattern, shard) pairs are recorded in DIR/done_<id>.txt; p3 checkpoints inside a unit (-ckpt/-resume)."""
import sys, os, json, time, argparse
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import numpy as np
import p3g, p3run as P
from p3g import GStruct

ap = argparse.ArgumentParser()
ap.add_argument('--id', required=True); ap.add_argument('--v', type=int, required=True); ap.add_argument('--hgens', required=True)
ap.add_argument('--rows', action='append', default=None); ap.add_argument('--threads', type=int, default=int(os.environ.get('THREADS', 4)))
ap.add_argument('--ns', type=int, default=3); ap.add_argument('--shards', type=int, default=1)
ap.add_argument('--workdir', default=os.path.join(HERE, 'work')); ap.add_argument('--countlists', action='store_true')
ap.add_argument('--slice', default=None, help='BENCHMARK ONLY s/n'); ap.add_argument('--dmult', action='store_true', help='spend Mult on D (p3 -dperms), B only reversal-canonical')
a = ap.parse_args()
os.makedirs(a.workdir, exist_ok=True)
gens = [tuple(int(x) for x in g.split(',')) for g in a.hgens.split(';')]
st = GStruct(a.v, hgens=gens)
pats = st.feasible_patterns()
if a.rows:
    want = {tuple(int(x) for x in r.split(',')) for r in a.rows}
    pats = [p for p in pats if p['rows'] in want]
print(f'{a.id}: v={a.v} |H|={st.Hord} OA={st.A.O} OB={st.B.O} C={st.C} |mult|={len(st.mult_reps)} |fix|={len(st.fix)} |Gamma_amb|={st.gamma_order()} patterns={[p["rows"] for p in pats]}', flush=True)
donef = os.path.join(a.workdir, f'done_{a.id}.txt')
done = set(open(donef).read().split('\n')) if os.path.exists(donef) else set()
binary = os.path.join(HERE, 'p3')
hitfile = os.path.join(a.workdir, f'hits_{a.id}.txt')
verified = []
for p in pats:
    rows = p['rows']; tag = '_'.join(map(str, rows))
    pf = os.path.join(a.workdir, f'prob_{a.id}_{tag}.txt')
    dpf = pf + '.dperms'
    if not os.path.exists(pf):
        if a.dmult: p3g.write_problem_dmult(st, rows, pf, dpf, seed=a.v)
        else: P.write_problem(st, rows, pf, canonB=True, canonD=True, psd=(True, True, True), seed=a.v)
    bfile = os.path.join(a.workdir, f'B_{a.id}_{tag}.bin')
    for sh in range(a.shards):
        unit = f'{tag}:{sh}/{a.shards}'
        if unit in done: continue
        ck = os.path.join(a.workdir, f'ck_{a.id}_{tag}_{sh}_{a.shards}.txt')
        t = time.time()
        extra = ['-slice', a.slice] if a.slice else []
        if a.countlists: extra += ['-countlists']
        bvalid = os.path.exists(bfile + '.ok')
        blist = None
        if a.dmult: extra += ['-dperms', dpf]
        if a.shards > 1 and not a.dmult and not a.slice and not a.countlists:       # the B list is the same in every shard: build it once, reuse it
            if bvalid: blist = bfile
            else: extra += ['-bsave', bfile]
        r = P.run_p3(pf, binary=binary, threads=a.threads, ns=a.ns, shard=(sh, a.shards) if a.shards > 1 else None, ckpt=ck,
                     resume=os.path.exists(ck) and os.path.getsize(ck) > 0, hits=hitfile, extra=extra, blist=blist)
        if a.shards > 1 and not a.dmult and not bvalid and not a.slice and not a.countlists:
            open(bfile + '.ok', 'w').close()
        print(f'  rows {rows} shard {sh}/{a.shards}: lists(A,B,D) {r["lists"]} done {r["done"]} hits {len(r["hits"])} [{time.time()-t:.0f}s]', flush=True)
        if a.slice or a.countlists:
            continue
        for h in r['hits']:
            A, B, D = st.masks_to_seqs(*h)
            if st.exact_check(A, B, D) and st.matrix_check(A, B, D):
                rec = dict(id=a.id, v=a.v, group=[int(x) for x in st.G.mods], hgens=[list(g) for g in gens], rows=[int(A.sum()), int(B.sum()), int(D.sum())],
                           seqs=[''.join('+' if x > 0 else '-' for x in s) for s in (A, B, D)], order=4 * a.v, verified='group PAF + symmetric Hadamard matrix check')
                verified.append(rec)
                with open(os.path.join(a.workdir, f'HIT_{a.id}.json'), 'w') as f: json.dump(rec, f)
                print('  VERIFIED HIT', flush=True)
            else:
                print('  engine hit FAILED exact verification', flush=True)
        with open(donef, 'a') as f: f.write(unit + '\n')
print(f'{a.id}: finished; verified hits {len(verified)}', flush=True)
