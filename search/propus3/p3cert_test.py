"""Test of the hit -> matrix -> certificate conversion on the repository's orders 532 and 876, and through the whole
pipeline (engine hit -> certificate -> the repository's verify/verify.py, run unmodified as a subprocess)."""
import json, subprocess, sys, tempfile, time
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from p3struct import Structure
import p3run as P, p3cert as K
REPO = HERE.parent.parent          # the repository root when this directory is search/propus3
sys.path.insert(0, str(REPO / 'verify'))
import verify as RV           # read-only import of the repository's verifier (used for a matrix comparison only)

out = HERE / 'out'
ok = True
for order, v, gens in ((532, 133, [26]), (876, 219, [4])):
    st = Structure(v, gens)
    repo_dir = REPO / 'symmetric' / str(order)
    rc = json.loads((repo_dir / 'certificate.json').read_text())
    seqs = []
    for reps in rc['blocks']:
        x = np.ones(v, dtype=np.int64)
        for r in reps:
            for h in st.H: x[r * h % v] = -1
        seqs.append(x)
    # (1) regenerate the repository's certificate from its own sequences
    cert, H = K.make_certificate(st, gens, *seqs)
    same_json = cert == rc
    text = K.dump_cert(cert)
    same_text = text == (repo_dir / 'certificate.json').read_text()
    # (2) our matrix equals the repository verifier's matrix
    Hrepo, _ = RV.build_propus(rc)
    same_H = np.array_equal(H, Hrepo)
    seq_txt = '\n'.join(K.seq_line(x) for x in seqs) + '\n'
    same_seq = seq_txt == (repo_dir / 'sequences.txt').read_text()
    print(f'{order}: certificate dict identical {same_json}, certificate.json text identical {same_text}, matrix identical to verify.build_propus {same_H}, sequences.txt identical {same_seq}', flush=True)
    ok &= same_json and same_text and same_H and same_seq
    # (3) pipeline: engine hit -> certificate -> repository verify.py CLI
    pf = Path(tempfile.mkdtemp(prefix='p3cert_')) / 'prob.txt'
    rows = tuple(rc['rowsums'])
    P.write_problem(st, rows, pf, canonB=True, canonD=True)
    res = P.run_p3(pf, threads=4, ns=2)
    assert len(res['hits']) > 0
    ma, mb, md = sorted(res['hits'])[0]
    A, B, D = P.masks_to_seqs(st, ma, mb, md)
    d = out / 'symmetric' / f'{order}_from_engine'
    K.write_certificate(d, st, gens, A, B, D, matrix=False)
    r = subprocess.run([sys.executable, str(REPO / 'verify' / 'verify.py'), str(d)], capture_output=True, text=True)
    print('  repository verify.py on the engine-derived certificate:', r.stdout.strip()[:200], r.stderr.strip()[:200], flush=True)
    ok &= r.returncode == 0 and ' ok ' in r.stdout
    # (4) the same for the regenerated repository certificate (directory named symmetric/<order>)
    d2 = out / 'symmetric' / f'{order}'
    K.write_certificate(d2, st, gens, *seqs, matrix=False)
    r2 = subprocess.run([sys.executable, str(REPO / 'verify' / 'verify.py'), str(d2)], capture_output=True, text=True)
    print('  repository verify.py on the regenerated repository certificate:', r2.stdout.strip()[:200], flush=True)
    ok &= r2.returncode == 0 and ' ok ' in r2.stdout
print('ALL OK' if ok else 'FAILED')
sys.exit(0 if ok else 1)
