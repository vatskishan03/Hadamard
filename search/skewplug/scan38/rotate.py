"""rotate.py: time-sliced rotation of exh38.py over the shapes of data/queue/shapes.json (env ROT_QUEUE: another list; stop: touch runs/rotate.STOP; per-shape logs runs/x<idx>_<v>.log) (resumable; every slice resumes from runs/NAME.slots.jsonl).
usage (in search/skewplug/scan38): python3 -u rotate.py DEADLINE_EPOCH PASSES_JSON IDX1,IDX2,...    PASSES_JSON = list of [MAXCLASS, base_slice_seconds]
Per pass: every unfinished shape gets one slice, length = base * clamp(yield / median yield, 0.3, 3) from the second pass on (yield = new classes per wall second of the shape's previous slices, with
completions weighted x1), i.e. the rotation follows the measured yield of new classes per second; a shape whose pass ends with ALL SLOTS FINISHED is left until the next pass (larger MAXCLASS)."""
import json, subprocess, sys, time, os, statistics
HERE = os.path.dirname(os.path.abspath(__file__)); W = os.path.abspath(os.path.join(HERE, '..'))      # search/skewplug: exh38.py runs there
Q = os.environ.get('ROT_QUEUE', os.path.join(HERE, 'data', 'queue', 'shapes.json'))
os.makedirs(W + '/runs', exist_ok=True)
deadline = float(sys.argv[1]); passes = json.loads(sys.argv[2]); idxs = [int(x) for x in sys.argv[3].split(',')]
shapes = {s['idx']: s for s in json.load(open(Q))}
def stats(name):
    p = f'{W}/runs/{name}.slots.jsonl'; nn = nk = 0; t = 0.0
    if os.path.exists(p):
        for l in open(p):
            try: r = json.loads(l)
            except Exception: continue
            nn += r.get('nnew', 0); nk += r.get('nok', 0); t += r.get('t_slot', 0)
    return nn, nk, t
yld = {}
def log(*a):
    print(time.strftime('%H:%M:%S', time.gmtime()), *a, flush=True)
for pi, (maxclass, base) in enumerate(passes):
    log(f'pass {pi}: MAXCLASS {maxclass} base slice {base}s')
    med = statistics.median(yld.values()) if yld else None
    for idx in idxs:
        if time.time() > deadline or os.path.exists(W + '/runs/rotate.STOP'): log('deadline/STOP'); sys.exit(0)
        s = shapes[idx]; name = f"x{idx}_{s['v']}"
        sl = base
        if med and idx in yld: sl = base * min(3.0, max(0.3, yld[idx] / med if med > 0 else 1.0))
        sl = min(sl, deadline - time.time())
        n0, k0, t0 = stats(name)
        env = dict(os.environ, EX_MAXCLASS=str(maxclass), EX_SLOT_SECS='90', OMP_NUM_THREADS='1', PROD_CONTINUE='1')
        cmd = ['nice', '-n', '5', 'python3', '-u', 'exh38.py', name, str(int(sl)), json.dumps(s['K']), str(s['q']), json.dumps(s['hu']), s['est'], '100000', 'rho', 'skf']
        t1 = time.time()
        with open(f'{W}/runs/{name}.log', 'a') as lf:
            pr = subprocess.run(['bash', '-c', 'ulimit -v 3500000; cd ' + W + '; exec "$@"', 'x'] + cmd, env=env, stdout=lf, stderr=subprocess.STDOUT)
        n1, k1, t_ = stats(name); dt = time.time() - t1
        yld[idx] = (n1 - n0) / max(dt, 1.0)
        fin = 'ALL SLOTS FINISHED' in open(f'{W}/runs/{name}.log').read().splitlines()[-1] if os.path.getsize(f'{W}/runs/{name}.log') else False
        log(f'{name} v={s["v"]} slice {sl:.0f}s wall {dt:.0f}s: +{n1 - n0} classes (total {n1}), +{k1 - k0} completions (total {k1}), yield {yld[idx]:.2f}/s, finished={fin}')
log('all passes done')
