"""Low-energy four-list k-tree on a prime field with GENERATED lists (lg_pt) and twin-deduplicated lists, plugged into kt4.

  lg_lowe.py v e cfg --LA 3e5 --LB 1e8 --LC 6e7 --LD 6e7 --nB 1.5e7 --nC 1e7 --nD 1e7 [--secs-gen 3600] [--threads 4]
             [--T 0.42] [--nsig 2000] [--join-core-h 2] [--K 1] [--seed 1] [--workdir DIR] [--keep-lists]

L* are ball sizes in members (twins counted), n* the numbers of DISTINCT (reversal-canonical) members to collect by the chain
(stopping at the time cap); block A (2^(e/2) members) is enumerated exactly below its threshold (Gray code).
Writes HIT_<tag>.json (only on a verified solution) and appends the hit to $HADAMARD_WORK/skew_hits.jsonl."""
import argparse, json, math, os, struct, subprocess, sys, tempfile, time, resource
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import numpy as np
import lg_lib as LG, gen_run as GR
import cp_lowe as CL, cyclo_product as CP, workdir as W
from scipy import stats

def seq_from_mask(S, mask, f):
    """+-1 sequence of length v: f at 0, the sign of class c (bit c of the class mask set: -1) on the elements of class c."""
    x = 1 - 2 * np.array([(mask >> c) & 1 for c in range(S.e)], dtype=np.int64)
    cls = np.asarray(S.cls)
    seq = np.empty(S.v, dtype=np.int64); seq[0] = f
    seq[1:] = x[cls[1:]]
    return seq


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('v', type=int); ap.add_argument('e', type=int); ap.add_argument('cfg', type=int)
    ap.add_argument('--LA', type=float, default=3e5); ap.add_argument('--LB', type=float, default=1e8); ap.add_argument('--LC', type=float, default=6e7); ap.add_argument('--LD', type=float, default=6e7)
    ap.add_argument('--nB', type=float, default=1.5e7); ap.add_argument('--nC', type=float, default=1e7); ap.add_argument('--nD', type=float, default=1e7)
    ap.add_argument('--secs-gen', type=float, default=3600.0); ap.add_argument('--threads', type=int, default=4); ap.add_argument('--T', type=float, default=0.42)
    ap.add_argument('--nsig', type=int, default=2000); ap.add_argument('--mcpairs', type=int, default=4000000); ap.add_argument('--maxpairs', type=float, default=3e8)
    ap.add_argument('--join-core-h', type=float, default=0.0, help='stop kt4 after this many core-hours of wall*threads (0: no cap)')
    ap.add_argument('--K', type=int, default=1); ap.add_argument('--seed', type=int, default=1); ap.add_argument('--samples', type=int, default=400000)
    ap.add_argument('--workdir', default=''); ap.add_argument('--keep-lists', action='store_true'); ap.add_argument('--S', default='')
    ap.add_argument('--kB', type=float, default=0); ap.add_argument('--kC', type=float, default=0); ap.add_argument('--kD', type=float, default=0)
    ap.add_argument('--dry-curve', default='', help='with --dry: comma-separated numbers of cells; prints predicted pairs and expected key-quadruples for each')
    ap.add_argument('--plant', default='', help='validation: add the members of the first solution of this hits file to the lists and join at its split value only')
    ap.add_argument('--dry', action='store_true', help='stop before kt4 (print the estimate and the predicted pairs of the chosen cells)')
    ap.add_argument('--share-bc', action='store_true', help='use the list of block B also as block C (same row sum: same structure)')
    ap.add_argument('--resume', action='store_true', help='reuse the block lists already present in the work directory')
    ap.add_argument('--lists-from', default='', help='directory of lists written by a previous run (skips generation)')
    a = ap.parse_args()
    t0 = time.time(); log = lambda s: print('%s [%.0fs]' % (s, time.time() - t0), flush=True)
    wd = a.workdir or os.path.join(os.environ['HADAMARD_WORK'], 'lg_%d_%d_%d' % (a.v, a.e, a.cfg)); os.makedirs(wd, exist_ok=True)
    S, y, N = LG.make_structure(a.v, a.e, a.cfg); e = a.e; D = S.nb
    log('v=%d e=%d F=%d config %d y=%s D=%d' % (a.v, e, S.F, a.cfg, y, D))
    import pickle
    mp = os.path.join(wd, 'mom_%d_%d.pkl' % (a.samples, a.seed))
    if os.path.exists(mp): mom = pickle.load(open(mp, 'rb'))
    else:
        mom = LG.moments(S, N, y, nsamp=a.samples, seed=a.seed); pickle.dump(mom, open(mp, 'wb'))
    Lball = [a.LA, a.LB, a.LC, a.LD]; ntarget = [None, a.nB, a.nC, a.nD]
    lists = [None] * 4; masks = [None] * 4; cpu0 = resource.getrusage(resource.RUSAGE_CHILDREN)
    lg = os.path.join(HERE, 'lg_pt')
    for k in range(4):
        pth = os.path.join(wd, 'L%d.npz' % k)
        src = os.path.join(a.lists_from, 'L%d.npz' % k) if a.lists_from else pth
        if os.path.exists(src) and (a.lists_from or a.resume):
            z = np.load(src); masks[k], lists[k] = z['m'], z['k']; log('block %d: %d members loaded from %s' % (k, len(masks[k]), src)); continue
        if k == 2 and a.share_bc and mom[1]['f'] == mom[2]['f'] and mom[1]['s'] == mom[2]['s']:
            masks[2], lists[2] = masks[1], lists[1]; log('block 2: shares the list of block 1'); continue
        ck, Wk, Z = LG.energy_setup(mom, k)
        delta = Wk @ (mom[k]['mu'] - ck); nc = float(delta @ delta)
        size = 2.0 ** (e // 2) if k == 0 else math.comb(e, (e + mom[k]['s']) // 2)
        p = min(1.0, Lball[k] / size)
        tau = 1e30 if p >= 1 else GR.tau_for(p, Z, nc)
        inp = os.path.join(wd, 'in%d.bin' % k); outp = os.path.join(wd, 'out%d.bin' % k)
        LG.export_block(inp, S, N, mom, k, ck=ck, W=Wk)
        if k == 0:
            maxd = int(min(size / 2, 4.0 * p * size / 2)) + 1000
            cmd = [lg, inp, outp, repr(tau), repr(tau), '1e9', '1000000000', '1', str(a.seed), '1', '0.4', '0.4', '1000', str(maxd), os.path.join(wd, 'stat0'), '3']
        else:
            tgt = int(ntarget[k]); maxd = int(1.05 * tgt) + 1000
            cmd = [lg, inp, outp, repr(tau), repr(tau), str(a.secs_gen), str(tgt), str(a.threads), str(a.seed * 10 + k), '1', str(a.T), str(a.T), '1000', str(maxd), os.path.join(wd, 'stat%d' % k), '0']
        log('block %d: size %.3g ball %.3g (p %.3g) rank %d nc %.2f tau %.3f -> collect %s' % (k, size, Lball[k], p, Z, nc, tau, 'all' if k == 0 else int(ntarget[k])))
        r = subprocess.run(cmd, capture_output=True, text=True); log('  ' + r.stdout.strip())
        if r.returncode: print(r.stderr[-2000:]); sys.exit(1)
        rec = LG.read_records(outp, D); os.remove(outp)
        masks[k] = rec['mask'].copy(); lists[k] = rec['key'].copy(); del rec
        np.savez(pth, m=masks[k], k=lists[k])
    for k in (1, 2, 3):                                                    # keep the lowest-energy members when asked
        keepn = [a.kB, a.kC, a.kD][k - 1]
        if keepn and keepn < len(masks[k]):
            ck, Wk, Z = LG.energy_setup(mom, k)
            en = np.empty(len(masks[k]))
            for i in range(0, len(en), 1 << 20): en[i:i + (1 << 20)] = (((lists[k][i:i + (1 << 20)].astype(float) - ck) @ Wk.T) ** 2).sum(1)
            idx = np.argpartition(en, int(keepn))[:int(keepn)]
            log('block %d: kept the %d lowest of %d: energy cut %.3f (max %.3f)' % (k, len(idx), len(en), en[idx].max(), en.max()))
            masks[k], lists[k] = masks[k][idx], lists[k][idx]
    plant_keys = None
    if a.plant:
        rec0 = json.loads(open(a.plant).readline()); plant_keys = []
        sq0 = [np.array([1 if c == '+' else -1 for c in s_]) for s_ in rec0['seqs']]
        for k in range(4):
            x = np.array([sq0[k][S.reps[1 + c]] for c in range(e)])
            cm = int(LG.canon(np.array([sum(1 << c for c in range(e) if x[c] < 0)], dtype=np.uint64), e)[0])
            paf = CL.paf_vector(S, [sq0[k]]); key = np.array([paf[S.reps[1 + c]] for c in range(D)], dtype=np.int16)
            plant_keys.append(key.astype(np.int64))
            masks[k] = np.concatenate([masks[k], np.array([cm], dtype=np.uint64)]); lists[k] = np.concatenate([lists[k], key[None, :]])
        log('planted the four members of the known solution (sum of keys %s)' % (sum(plant_keys),))
    ru = resource.getrusage(resource.RUSAGE_CHILDREN); log('generation CPU %.0f s' % (ru.ru_utime + ru.ru_stime - cpu0.ru_utime - cpu0.ru_stime))
    sizes = [len(m) for m in masks]
    rng = np.random.default_rng(a.seed)
    samples = [l if len(l) <= 2000000 else l[rng.choice(len(l), 2000000, replace=False)] for l in lists]
    est, rank, off = CL.gauss_estimate(samples, np.zeros(D), sizes=sizes)
    log('Gaussian estimate of key-quadruples inside the lists %s: 10^%.2f (rank %d, null offset %.2f)' % (sizes, est, rank, off))
    var = sum(CL.colvar(l) for l in lists); absmax = np.max([np.maximum(l.max(0).astype(np.int64), -l.min(0).astype(np.int64)) for l in lists], axis=0)
    fits = absmax <= 1000
    S_idx = [int(x) for x in a.S.split(',')] if a.S else sorted([int(c) for c in np.argsort(-var) if fits[c]][:5])
    sigmas, info = CL.choose_sigmas(samples, S_idx, np.zeros(D, np.int64), a.nsig, a.mcpairs, rng, sizes=sizes)
    log('split values: %d chosen from %d; share of solutions %.3g; predicted pairs AB %.3g CD %.3g' % (len(sigmas), info['candidates'], info['share'], info['ab'], info['cd']))
    if plant_keys is not None:
        sigmas = np.array([(plant_keys[0] + plant_keys[1])[S_idx]], dtype=np.int64); log('join only at the split value of the planted solution: %s' % sigmas[0])
    if a.dry and a.dry_curve:
        for ns in [int(x) for x in a.dry_curve.split(',')]:
            sg, inf = CL.choose_sigmas(samples, S_idx, np.zeros(D, np.int64), ns, a.mcpairs, np.random.default_rng(a.seed), sizes=sizes)
            log('CURVE cells %d: pairs AB %.3g CD %.3g total %.3g; share %.4g; expected key-quadruples %.4g' % (ns, inf['ab'], inf['cd'], inf['ab'] + inf['cd'], inf['share'], 10 ** est * inf['share']))
    log('expected key-quadruples in the %d chosen cells: %.3g (10^%.2f x share %.3g)' % (len(sigmas), 10 ** est * info['share'], est, info['share']))
    if a.dry: return
    path = os.path.join(wd, 'lists.kt4')
    with open(path, 'wb') as f:
        f.write(struct.pack('5i', *sizes, D)); f.write(np.zeros(D, np.int16).tobytes())
        for l in lists: f.write(np.ascontiguousarray(l, dtype=np.int16).tobytes())
    del samples
    import re
    errpath = os.path.join(wd, 'kt4_cells.err'); done = 0
    if os.path.exists(errpath) and not a.plant:                           # resume: the ranking is deterministic, skip the cells already searched
        offset = 0
        for ln in open(errpath):
            if ln.startswith('LGOFFSET '): offset = int(ln.split()[1]); done = max(done, offset)
            m_ = re.match(r'kt4: sigma (\d+)/', ln)
            if m_: done = max(done, offset + int(m_.group(1)))
        if done: log('resuming the join after %d cells already searched' % done); sigmas = sigmas[done:]
    spath = path + '.sig'; np.savetxt(spath, sigmas, fmt='%d')
    O = CL.Orbits(S, 0); gens = CL.generators(S, O)
    env = dict(os.environ); env['OMP_NUM_THREADS'] = str(a.threads)
    errf = open(errpath, 'a'); errf.write('LGOFFSET %d\n' % done); errf.flush()
    p = subprocess.Popen([W.bin('kt4'), path, ','.join(map(str, S_idx)), '@' + spath, str(a.K), str(int(a.maxpairs))], stdout=subprocess.PIPE, stderr=errf, text=True, env=env)
    tj = time.time(); found = []
    import threading
    timer = None
    if a.join_core_h > 0:
        timer = threading.Timer(a.join_core_h * 3600.0 / a.threads, p.terminate); timer.daemon = True; timer.start()
    fbits = [0 if mom[k]['f'] == 1 else 1 for k in range(4)]
    for line in p.stdout:
        if line.startswith('TOTAL'): log('kt4: ' + line.strip()); continue
        if not line.startswith('SOL'): continue
        sol = list(map(int, line.split()[1:]))
        seqs = [seq_from_mask(S, int(masks[k][sol[k]]), mom[k]['f']) for k in range(4)]
        paf_ok, skew_ok, inv_ok = CL.check_hit(S, seqs, gens)
        log('HIT: autocorrelations %s skew %s invariant %s, row sums %s' % (paf_ok, skew_ok, inv_ok, [int(x.sum()) for x in seqs]))
        if paf_ok and skew_ok and inv_ok:
            rec = dict(n=a.v, gens=gens, pattern=[int(x.sum()) for x in seqs], hadamard=True, skew=True, method='lg_lowe e=%d cfg=%d seed=%d' % (e, a.cfg, a.seed),
                       seqs=[''.join('+' if z > 0 else '-' for z in x) for x in seqs])
            with open(W.hits(), 'a') as fh: fh.write(json.dumps(rec) + '\n')
            json.dump(rec, open('HIT_%d_%d_%d.json' % (a.v, a.e, a.cfg), 'w'))
            found.append(rec)
            p.terminate(); break
    p.wait(); errf.close()
    log('kt4 wall %.0f s; %d verified solutions' % (time.time() - tj, len(found)))
    if not a.keep_lists:
        for k in range(4):
            for nm in ('L%d.npz' % k, 'in%d.bin' % k):
                try: os.remove(os.path.join(wd, nm))
                except OSError: pass
    for pth in (path, spath):
        try: os.remove(pth)
        except OSError: pass

if __name__ == '__main__':
    main()
