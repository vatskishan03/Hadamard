"""I/O helpers between pk.Shape (plain mode, cyclic K) and the C level-2 solver pl2.c"""
import numpy as np, subprocess, os, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))


def header_text(sh):
    assert len(sh.Kmods) == 1, 'cyclic K only'
    n = sh.N; q = sh.q; eps = sh.eps; T = sh.N2
    kI = [int(k) for k in sh.repsPM if k != 0]
    kQ = [int(k) for k in sh.repsQ]
    L = [f'{n} {q} {eps} {T}', f'{len(kI)} ' + ' '.join(map(str, kI)), f'{len(kQ)} ' + ' '.join(map(str, kQ))]
    for K in (sh.kAB, sh.kCD):
        L.append(str(K.nu))
        for i in range(K.nu):
            P = K.pos[i]
            L.append(f'{len(P)} ' + ' '.join(f'{int(p)} {int(K.sig[p])}' for p in P))
    return '\n'.join(L) + '\n'


def dmask(K, e):
    """defect flag per orbit from the e array (defect iff e != 0 at the orbit's first position)"""
    return [int(e[int(p[0])] != 0) for p in K.pos]


def instance_text(sh, blocks):
    """blocks: list of 4 (s, e) arrays"""
    L = []
    for (s, e) in blocks: L.append(' '.join(str(int(x)) for x in s))
    for b, (s, e) in enumerate(blocks):
        K = sh.kAB if b < 2 else sh.kCD
        L.append(' '.join(map(str, dmask(K, e))))
    return '\n'.join(L) + '\n'


def run_pl2(sh, insts, maxsol=5, tlimit=1e9, seed=1, exe=None, timeout=None):
    """insts: list of blocks (4 x (s,e)).  returns (list of per-instance dict, list of solutions (inst idx, [t0..t3]))"""
    if exe is None: exe = os.path.join(HERE, 'pl2_32' if sh.N - 1 <= 32 else 'pl2_64')
    fd, path = tempfile.mkstemp(suffix='.pl2in', dir=os.environ.get('PL2TMP', os.path.join(HERE, 'tmp')))
    with os.fdopen(fd, 'w') as f:
        f.write(header_text(sh)); f.write(f'{len(insts)}\n')
        for blocks in insts: f.write(instance_text(sh, blocks))
    try:
        out = subprocess.run([exe, path, str(maxsol), str(tlimit), str(seed)], capture_output=True, text=True, timeout=timeout).stdout
    finally:
        os.unlink(path)
    res = {}; sols = []
    for line in out.splitlines():
        if line.startswith('INST'):
            parts = line.split()
            i = int(parts[1]); d = dict(raw=line)
            for tok in parts[2:]:
                if tok.startswith('completions='): d['compl'] = int(tok.split('=')[1])
                if tok.startswith('t='): d['t'] = float(tok.split('=')[1])
            if 'compl' not in d: d['compl'] = 0
            res[i] = d
        elif line.startswith('SOL'):
            segs = line.split('|')
            i = int(segs[0].split()[1])
            ts = [np.array([int(x) for x in sg.split()], dtype=np.int16) for sg in segs[1:5]]
            sols.append((i, ts))
    return res, sols
