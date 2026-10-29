"""queue worker (in search/skewplug): python3 scan38/sched38.py SHAPES.json QUEUE.json SECS PASSNAME [EREL] [M]
QUEUE.json = ordered list of shape indices.  Each worker claims the next unclaimed index (mkdir claims/PASSNAME/idx), runs scan38.run_shape with the CUMULATIVE budget SECS (resumable: a shape whose
hdr says scan_secs >= SECS or COMPLETE is skipped), and goes on.  Several workers (<= 3 in total) may run the same queue.  After an interrupted run: rm -rf claims/PASSNAME and relaunch."""
import sys, os, json, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..'))      # scan38/ and the modules of search/skewplug
DATA = os.path.join(HERE, 'data')
import scan38
rows = json.load(open(sys.argv[1])); queue = json.load(open(sys.argv[2])); secs = float(sys.argv[3]); pas = sys.argv[4]
erel = int(sys.argv[5]) if len(sys.argv) > 5 else 5; M = int(sys.argv[6]) if len(sys.argv) > 6 else 60; wmode = os.environ.get('WMODE', 'top')
os.makedirs(f'claims/{pas}', exist_ok=True); os.makedirs('est', exist_ok=True)
stop = f'claims/{pas}.STOP'
for i in queue:
    if os.path.exists(stop): print('STOP file', flush=True); break
    try: os.mkdir(f'claims/{pas}/{i}')
    except FileExistsError: continue
    name = scan38.shape_name(i, rows[i], wmode); hdr = f'est/{name}.hdr.json'
    if os.path.exists(hdr):
        try:
            h = json.load(open(hdr))
            if h.get('status') in ('COMPLETE', 'EXACT_DEAD_identity') or float(h.get('scan_secs', 0)) >= secs: print('skip', i, h.get('status'), flush=True); continue
        except Exception: pass
    t0 = time.time()
    try: info = scan38.run_shape(i, rows[i], secs, erel, M, wmode=wmode)
    except Exception as ex: info = dict(idx=i, error=str(ex)[:200])
    print(json.dumps(dict(idx=i, name=name, status=info.get('status'), pool=info.get('pool'), cov=info.get('coverage'), secs=round(time.time() - t0, 1), err=info.get('error'))), flush=True)
print('queue done', flush=True)
