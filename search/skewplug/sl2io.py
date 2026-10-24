"""level 2 for the skew plug-in: gl2.c patched so that block 0 (A, skew: t even on {0} and on pairs O u -O) has its own t-orbit structure (kind 0), blocks 1..3 are free (kind 1)."""
import numpy as np, subprocess, os, sys, tempfile
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from gl2io import channels
TMP = os.environ.get('GL2TMP', os.path.join(HERE, 'tmp')); os.makedirs(TMP, exist_ok=True)


def header_text(sh):
    n = sh.N; r = sh.r; eps = sh.eps; T = sh.T2
    kI = [int(k) for k in sh.repsPM if k != 0]
    kQ = [int(k) for k in sh.repsQ]
    assert len(kI) <= 128 and len(kQ) <= 128 and n <= 640
    L = [f'{n} {r} {eps} {T}', str(len(kI))]
    for k in kI: L.append(' '.join(map(str, sh.shift[k])) + ' ' + ' '.join(map(str, sh.shneg[k])))
    L.append(str(len(kQ)))
    for k in kQ: L.append(' '.join(map(str, sh.shift[k])) + ' ' + ' '.join(map(str, sh.shneg[k])))
    # kind 0 : A
    L.append(str(len(sh.at_orbits)))
    for orb in sh.at_orbits: L.append(f'{len(orb)} ' + ' '.join(f'{p} {sg}' for (p, sg) in orb))
    K = sh.kCD
    L.append(str(K.nu))
    for i in range(K.nu):
        P = K.pos[i]
        L.append(f'{len(P)} ' + ' '.join(f'{int(p)} {int(K.sig[p])}' for p in P))
    ch = channels(sh)
    L.append(str(len(ch)))
    for (a, b, c, d) in ch:
        for arr in (a, b, c, d): L.append(' '.join(f'{x:.17g}' for x in arr))
    return '\n'.join(L) + '\n', len(ch)


def dmask_free(K, e): return [int(e[int(p[0])] != 0) for p in K.pos]


def instance_text(sh, blocks):
    """blocks: 4 x (s, e); block 0 = A"""
    L = []
    for (s, e) in blocks: L.append(' '.join(str(int(x)) for x in s))
    L.append(' '.join(map(str, sh.A_dflags(blocks[0][1]))))
    for b in range(1, 4): L.append(' '.join(map(str, dmask_free(sh.kCD, blocks[b][1]))))
    return '\n'.join(L) + '\n'


def run_sl2(sh, insts, maxsol=5, tlimit=1e9, seed=1, maxfib=1500000, timeout=None, stop_first=False, hdr=None):
    hdr, nch = hdr if hdr else header_text(sh)
    exe = os.path.join(HERE, 'sl2_32' if nch <= 32 else ('sl2_64' if nch <= 64 else 'sl2_128'))
    fd, path = tempfile.mkstemp(suffix='.sl2in', dir=TMP)
    with os.fdopen(fd, 'w') as f:
        f.write(hdr); f.write(f'{len(insts)}\n')
        for blocks in insts: f.write(instance_text(sh, blocks))
    try:
        proc = subprocess.run([exe, path, str(maxsol), str(tlimit), str(seed), str(maxfib), str(int(stop_first))], capture_output=True, text=True, timeout=timeout)
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
