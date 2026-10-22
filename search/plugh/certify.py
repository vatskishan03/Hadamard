"""Certify a K-level plug-in solution of a PLAIN shape: build the G-level +-1 sequences A,B,C,D on G = K x F_q (mixed radix [K mods..., p]*e),
check (entries, A/B symmetric, sum PAF = 4|G| delta, K K* = 2|G| I exactly via buildq.check_solution), write complex-hadamard/C<2|G|>.json into the
directory out/ (or [outdir]) with its SHA-256 from verify/verify_complex.py; run verify_complex.py when outdir is the repository's complex-hadamard/.
usage: certify.py pickle_file solution_index Kmods(json) q hg(json) [outdir]"""
import sys, os, json, pickle, shutil, subprocess, hashlib
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from pk import Shape
from buildq import check_solution, glevel_seq
from idx import GF, Abel, factor_prime_power


def certify(sh, sol, outdir=None, tag=''):
    ok, (seqs, (Kr, Ki)) = check_solution(sh, sol, verbose=True)
    if not ok: return None
    q = sh.q; p, ee = factor_prime_power(q)
    mods = list(sh.Kmods) + [p] * ee
    n = int(np.prod(mods)); order = 2 * n
    strs = [''.join('+' if v == 1 else '-' for v in f) for f in seqs]
    cert = dict(order=order, N=n, construction='szollosi-ssxx-group', group=mods, sequences=strs,
                note=f'plain Paley plug-in K={sh.Kmods}, q={q}, H={sh.hgens_raw}; K-level solution certified by search/plugh/certify.py')
    repo = os.environ.get('HADAMARD_REPO', os.path.abspath(os.path.join(HERE, '..', '..')))
    sys.path.insert(0, os.path.join(repo, 'verify'))
    import importlib.util
    spec = importlib.util.spec_from_file_location('vc_' + str(order) + tag, os.path.join(repo, 'verify', 'verify_complex.py'))
    vc = importlib.util.module_from_spec(spec); spec.loader.exec_module(vc)
    Kr2, Ki2 = vc.build(cert); vc.check(Kr2, Ki2)
    sha = hashlib.sha256(vc.text(Kr2, Ki2).encode()).hexdigest()
    cert['sha256'] = sha
    outdir = outdir or os.path.join(HERE, 'out')
    os.makedirs(outdir, exist_ok=True)
    path = os.path.join(outdir, f'C{order}.json')
    json.dump(cert, open(path, 'w'))
    ok = True
    if os.path.abspath(outdir) == os.path.abspath(os.path.join(repo, 'complex-hadamard')):
        r = subprocess.run([sys.executable, os.path.join(repo, 'verify', 'verify_complex.py'), str(order)], capture_output=True, text=True, cwd=repo)
        print(r.stdout.strip(), r.stderr.strip()[:300]); ok = 'ok' in r.stdout and r.returncode == 0
    return dict(order=order, path=path, sha256=sha, verify_ok=ok, group=mods)


if __name__ == '__main__':
    pk_ = pickle.load(open(sys.argv[1], 'rb')); i = int(sys.argv[2])
    Kmods = json.loads(sys.argv[3]); q = int(sys.argv[4])
    hg = [(tuple(u) if isinstance(u, list) else u, ph) for u, ph in json.loads(sys.argv[5])]
    sh = Shape(Kmods, q, hg, mode='plain')
    sol = [tuple(np.array(a, dtype=np.int64) for a in blk) for blk in pk_[i]]
    print(certify(sh, sol, sys.argv[6] if len(sys.argv) > 6 else None))
