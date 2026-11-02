"""Checkpoint/resume test: kill the join part-way, resume from the checkpoint file, and compare the union of hits with the
uninterrupted run (order 532 structure, PSD filters off so the join is long enough to interrupt: 2.1e9 pairs)."""
import subprocess, sys, tempfile, time, signal
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import p3run as P
from p3struct import Structure
st = Structure(133, [26]); rows = (11, 5, 19)
tmp = Path(tempfile.mkdtemp(prefix='p3res_'))
P.write_problem(st, rows, tmp / 'prob.txt', canonB=False, canonD=False, psd=(False, False, False), seed=9)
binary = P.build_binary()
full = P.run_p3(tmp / 'prob.txt', binary, threads=2, ns=3)
want = set(full['hits'])
print('uninterrupted:', len(want), 'hits,', full['done']['pairs'], 'pairs', flush=True)
ck, hf = tmp / 'ck.txt', tmp / 'hits.txt'
cmd = [binary, str(tmp / 'prob.txt'), '-t', '2', '-ns', '3', '-bs', '1', '-ckpt', str(ck), '-hits', str(hf)]
p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
time.sleep(1.6); p.kill(); p.wait()
done1 = len(ck.read_text().split()) - 2
h1 = {tuple(int(x, 16) for x in l.split()) for l in hf.read_text().split('\n') if l.strip()} if hf.exists() else set()
print(f'killed after ~1.6 s: {done1} batches checkpointed, {len(h1)} hits so far', flush=True)
r = P.run_p3(tmp / 'prob.txt', binary, threads=2, ns=3, ckpt=ck, resume=True, hits=hf, extra=('-bs', '1'))
print('resumed run stderr:', [l for l in r['stderr'].splitlines() if 'resume' in l or 'done:' in l], flush=True)
h2 = {tuple(int(x, 16) for x in l.split()) for l in hf.read_text().split('\n') if l.strip()}
ok = h2 == want
print(f'union of hits {len(h2)} == uninterrupted {len(want)}: {ok}')
# a mismatching option set must refuse to resume
r3 = subprocess.run([binary, str(tmp / 'prob.txt'), '-t', '2', '-ns', '2', '-bs', '1', '-ckpt', str(ck), '-resume'], capture_output=True, text=True)
refused = r3.returncode != 0 and 'does not match' in r3.stderr
print('resume with different options refused:', refused)
print('ALL OK' if ok and refused and done1 > 0 else 'FAILED')
sys.exit(0 if ok and refused and done1 > 0 else 1)
