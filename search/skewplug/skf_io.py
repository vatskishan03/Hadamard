"""fast level 2 for the skew plug-in (cyclic K = Z_n): pl2f with block 0 = skew block (kind 0), blocks 1..3 free (kind 1).  Input format of pl2.c (lags instead of shift tables)."""
import numpy as np, subprocess, os, sys, tempfile, re
HERE = os.path.dirname(os.path.abspath(__file__))
TMP = os.environ.get('GL2TMP', os.path.join(HERE, 'tmp')); os.makedirs(TMP, exist_ok=True)
from sl2io import instance_text


def skf_header(sh):
    assert len(sh.Kmods) == 1, 'cyclic K only'
    n = sh.N; q = sh.r; eps = sh.eps; T = sh.T2
    kI = [int(k) for k in sh.repsPM if k != 0]; kQ = [int(k) for k in sh.repsQ]
    L = [f'{n} {q} {eps} {T}', str(len(kI)) + ' ' + ' '.join(map(str, kI)), str(len(kQ)) + ' ' + ' '.join(map(str, kQ))]
    L.append(str(len(sh.at_orbits)))
    for orb in sh.at_orbits: L.append(f'{len(orb)} ' + ' '.join(f'{p} {sg}' for (p, sg) in orb))
    K = sh.kCD
    L.append(str(K.nu))
    for i in range(K.nu):
        P = K.pos[i]
        L.append(f'{len(P)} ' + ' '.join(f'{int(p)} {int(K.sig[p])}' for p in P))
    # channels: orbit representatives of psi = 1..n-1 under the multiplications by H_K (x^(h psi) = x^(psi) for H_K-invariant x)
    hs = [int(u[0]) for u in sh.hunits]; seen = set(); reps = []
    for psi in range(1, n):
        if psi in seen: continue
        orb = {psi}; fr = [psi]
        while fr:
            x = fr.pop()
            for h in hs:
                y = (x * h) % n
                if y not in orb: orb.add(y); fr.append(y)
        seen |= orb; reps.append(psi)
    L.append(str(len(reps)) + ' ' + ' '.join(map(str, reps)))
    return '\n'.join(L) + '\n'


def run_skf(sh, insts, maxsol=5, tlimit=1e9, seed=1, timeout=None, hdr=None, exe=None):
    hdr = hdr or skf_header(sh)
    nreps = int(hdr.strip().split('\n')[-1].split()[0]); exe = exe or os.path.join(HERE, 'skf2_32' if nreps <= 32 else 'skf2_64')
    fd, path = tempfile.mkstemp(suffix='.skfin', dir=TMP)
    with os.fdopen(fd, 'w') as f:
        f.write(hdr); f.write(f'{len(insts)}\n')
        for blocks in insts: f.write(instance_text(sh, blocks))
    try:
        proc = subprocess.run([exe, path, str(maxsol), str(tlimit), str(seed)], capture_output=True, text=True, timeout=timeout)
        out = proc.stdout
    finally:
        os.unlink(path)
    res = {}; sols = []
    for line in out.splitlines():
        if line.startswith('INST'):
            parts = line.split(); i = int(parts[1]); d = dict(raw=line, compl=0, status='ok')
            for tok in parts[2:]:
                if tok.startswith('completions='): d['compl'] = int(tok.split('=')[1])
                if tok.startswith('t='): d['t'] = float(tok.split('=')[1])
                if tok.startswith('fib='): d['fib'] = tok[4:]
                if tok in ('TIMEOUT', 'TOOBIG', 'rowfilter'): d['status'] = tok
            if len(parts) > 2 and parts[2] in ('rowfilter', 'TOOBIG'): d['status'] = parts[2]
            res[i] = d
        elif line.startswith('SOL'):
            segs = line.split('|'); i = int(segs[0].split()[1])
            sols.append((i, [np.array([int(x) for x in sg.split()], dtype=np.int16) for sg in segs[1:5]]))
    for i in range(len(insts)):
        if i not in res:     # no INST line: the solver died (e.g. 'oom cb' under ulimit -v); the instance is NOT decided
            res[i] = dict(raw=(proc.stderr or '')[-200:], compl=0, status='L2FAIL', rc=proc.returncode)
    return res, sols
