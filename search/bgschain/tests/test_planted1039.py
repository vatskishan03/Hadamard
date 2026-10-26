"""End-to-end check at full 1039 scale: plant a quadruple with the real coset sums (target = its own PAF sum), fix its compressions on the first
three levels (Z37, Q x Z37, Z7) and let the engine run level 4 (Q x Z7: DP + join) and the leaves; the planted quadruple must be found."""
import json, sys, os, time, tempfile, shutil
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import numpy as np
import hlib, hchain
base = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'configs', 'n1039_prod.json')))
P0 = hlib.Problem(base); rng = np.random.default_rng(77); G = P0.G; v = G.v

def sample(ytarget):
    for attempt in range(2000):
        Y = 2 * rng.integers(0, 2, P0.O) - 1
        for it in range(400):
            c = np.array([(Y[P0.orb] * (P0.cos == cc)).sum() for cc in range(4)]); err = np.abs(c - ytarget).sum()
            if err == 0: return Y
            o = rng.integers(P0.O); Y2 = Y.copy(); Y2[o] = -Y2[o]
            c2 = np.array([(Y2[P0.orb] * (P0.cos == cc)).sum() for cc in range(4)])
            if np.abs(c2 - ytarget).sum() <= err: Y = Y2

masks = []; seqs = []
for i in range(4):
    Y = sample(np.array(base['ys'][i])); masks.append(int(sum(1 << o for o in range(P0.O) if Y[o] > 0))); seqs.append(Y[P0.orb])
shift = np.array([G.translate(s) for s in range(v)]); tgt = np.zeros(v, dtype=np.int64)
for x in seqs: tgt += (x[shift] * x[None, :]).sum(1)
cfg = dict(base); cfg.pop('normalizer'); cfg.pop('comment', None); cfg.update(name='n1039_planted', tgt=tgt.tolist())
tmpd = tempfile.mkdtemp(prefix='hchain_planted_')
ch = hchain.Chain(cfg, tmpd, verbose=True); P = ch.P
fixed = [(P.levels[l], [P.levels[l].of_masks(np.array([m], dtype=np.uint64), P.O)[0] for m in masks]) for l in range(3)]
st = hchain.Stats(); ctx = dict(mode='run', found=[], complete=True, stop=False, w=0.0, mult=1)
t = time.time(); ch.node(3, fixed, [False] * 4, st, ctx, [0])
found = [tuple(m) for m in ctx['found']]
print('time %.1fs' % (time.time() - t), {k: v for k, v in st.items() if k in ('sols_L4', 'leaves', 'solutions_unordered', 'solutions_ordered')})
print('planted quadruple found:', tuple(masks) in found, '| complete:', ctx['complete'])
shutil.rmtree(tmpd, ignore_errors=True)
