"""nm_5912.py -- near-miss statistics of level 2 at the 5912 plug-in (calibration of the completion probability, README "Order 5912").
Level 1 as in prod_5912.py (rng = default_rng([SEED, round]), token = round mod 4, no deduplication), then level 2 with gl2n_* (gl2.c built with -DNEARMISS) on chunks of 100 instances:
OMIT randomly chosen constraints of the 14 of level 2 (7 I-lags, 7 Q-lags) are dropped from the join per chunk and the PSD filters are off, so every t-quadruple that satisfies the other 14 - OMIT
is counted and its residuals on the dropped ones are recorded (a quadruple with all residuals 0 is a solution: it is saved as OUTPREFIX.SOL_<time>.json).
usage: nm_5912.py OUTPREFIX SEED OMIT ROUNDS [MAXINST]      one JSON line per instance in OUTPREFIX.jsonl: seed, rd, maxinst, al, omit, i, nm (events x 2: gauge t -> -t), compl, NM (residuals);
MAXINST level-1 instances per round (default 1000).  nm_cmp.py compares the counts with the Gaussian-lattice prediction."""
import sys, os, json, time, subprocess, tempfile
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from gs_samp import make_shape, parse_allocs
from gl1 import L1
from gl2io import header_text, instance_text, TMP
from prod_5912 import SHAPE, TOKENS, LA, LC, level1

CHUNK = 100


def run_gl2n(sh, insts, omit, seed):
    """gl2n on a list of instances: ({inst: events}, {inst: completions}, {inst: [residual lists]}, [(inst, [t0..t3])])"""
    hdr, nch = header_text(sh)
    exe = os.path.join(HERE, 'gl2n_32' if nch <= 32 else ('gl2n_64' if nch <= 64 else 'gl2n_128'))
    fd, path = tempfile.mkstemp(suffix='.gl2in', dir=TMP)
    with os.fdopen(fd, 'w') as f:
        f.write(hdr); f.write(f'{len(insts)}\n')
        for blocks in insts: f.write(instance_text(sh, blocks))
    env = dict(os.environ, GL2_OMIT=','.join(map(str, omit)), GL2_NOPSD='1')
    try: res = subprocess.run([exe, path, '5', '1e9', str(seed), '1500000', '0'], capture_output=True, text=True, env=env)
    finally: os.unlink(path)
    if res.returncode: raise RuntimeError(f'gl2n exit {res.returncode}: {res.stderr[:200]}')
    out = res.stdout
    nm, compl, NM, sols = {}, {}, {}, []
    for line in out.splitlines():
        p = line.split()
        if p[0] == 'INST': compl[int(p[1])] = int(next((x.split('=')[1] for x in p if x.startswith('completions=')), 0))
        elif p[0] == 'NMTOT': nm[int(p[1])] = int(p[2])
        elif p[0] == 'NM': NM.setdefault(int(p[1]), []).append([int(x) for x in p[2:]])
        elif p[0] == 'SOL': segs = line.split('|'); sols.append((int(segs[0].split()[1]), [[int(x) for x in sg.split()] for sg in segs[1:5]]))
    return nm, compl, NM, sols


def main():
    outp = sys.argv[1]; seed = int(sys.argv[2]); osz = int(sys.argv[3]); rounds = int(sys.argv[4])
    maxinst = int(sys.argv[5]) if len(sys.argv) > 5 else 1000
    sh = make_shape(*SHAPE); l1 = L1(sh); allocs = parse_allocs(TOKENS, LA, LC); pcache = {}
    nk = len([k for k in sh.repsPM if k != 0]) + len(sh.repsQ)
    fo = open(outp + '.jsonl', 'a'); t0 = time.time()
    for rd in range(rounds):
        rng = np.random.default_rng([seed, rd])
        al, cnt, S, E = level1(sh, l1, allocs, pcache, rng, rd, maxinst)
        for c0 in range(0, len(S), CHUNK):
            sl = range(c0, min(len(S), c0 + CHUNK))
            omit = sorted(rng.choice(nk, osz, replace=False).tolist())
            insts = [[(S[i, b], E[i, b]) for b in range(4)] for i in sl]
            nm, compl, NM, sols = run_gl2n(sh, insts, omit, int(rng.integers(1, 1 << 30)))
            for j, i in enumerate(sl):
                fo.write(json.dumps(dict(seed=seed, rd=rd, maxinst=maxinst, al=list(al), omit=omit, i=int(i), nm=nm.get(j), compl=compl.get(j), NM=NM.get(j, []))) + '\n')
            fo.flush()
            print(f'[{time.time()-t0:.0f}s] round {rd} {al} chunk {c0} omit {omit} events {sum(nm.values()) // 2} completions {sum(compl.values())}', flush=True)
            if sols:
                json.dump([(sl[j], ts) for j, ts in sols], open(f'{outp}.SOL_{int(time.time())}.json', 'w')); print('*** SOLUTION FOUND (level 2) in the near-miss run', flush=True)


if __name__ == '__main__': main()
