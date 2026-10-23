"""exh38.py: symmetry-reduced EXPLICIT-LIST sweep of the slots of a K-symmetric skew Paley shape (cyclic or general K; unit reduction of A for cyclic K).

prod_dedup.py samples lists of B, C, D per round and joins them against the exact A list; on the K-symmetric census shapes of search/NOTES.md section 91 the block spaces are small, the lists nearly cover them and every
class of level-1 instances is found ~2000 times (the exact block symmetries of prod_dedup's `inst_keys`).  Here, per (block, defect weight, m-weight, row sums) the list is generated
until it covers the block space (adaptive, EX_COV), reduced to ONE representative per orbit of the block symmetries (reflection k -> -k, negation, translation by the points of K fixed
by H_K; the level-1 system is a sum of per-block autocorrelation vectors, so each block may be replaced by any image), cached (memory + disk, shared by all jobs of the shape), and ALL
row-sum patterns of a slot are joined with the exact gl1.c join.  Further exact reductions: A is reflected/unit-reduced (the units of K = Z_n act on all four blocks at once; A is replaced by a
representative of its orbit, EX_UNITS=1), patterns that differ by exchanging two blocks of equal type are joined once.  The new classes (inst_keys of prod_dedup, same 64-bit keys, same
runs/NAME.seen) go to the exact level 2 (skf2 / sl2) in batches; a hit is certified (sk_cert) and the run continues.
usage (from this directory):  python3 -u exh38.py NAME SECS 'K' q 'hunits' SLOTFILE TOPN ORDER ENGINE
  SLOTFILE  est/xxx.jsonl (records al, ms, logN1, logE, logrho)      ORDER  rho | E | mix | file      ENGINE  skf | sl2
Environment: EX_COV (raw coverage of a generated list, default 0.9), EX_MAXDRAW (default 600000), EX_EPOCH (seed epoch, default 0), EX_MINRHO / EX_MINE (slot filters on the estimates),
EX_CACHE_MB (memory list cache, default 500), EX_L2CAP (level-2 instances per call, default 400), EX_BATCH (raw instances per level-2 batch, default 30000), EX_SLOT_SECS (time per slot
and visit, default 900; unfinished slots are resumed by the next pass), EX_UNITS (default 1), EX_SLOTS (a-b: only slots a..b of the ordered list), EX_SHARD (i/n: slots with index % n == i),
EX_PERSIST_MB (only lists up to this size, both arrays, are kept on disk, default 3), EX_JOINCAP (skip a pattern whose join costs more than this many pair evaluations, default 2e8; counted as nskip in the slot log), EX_MAXCLASS (skip slots with more estimated classes), EX_DUMP (directory: dump the instances of each slot).  runs/NAME.STOP ends the loop; runs/NAME.slots.jsonl is the per-slot log (resumable)."""
import sys, os, time, json, math, random, warnings, ctypes, collections, re
warnings.filterwarnings('ignore')
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from sk import SShape
from sk_pat import row_patterns_skew
from gl1 import L1, BlockSpec, _lib
from sl2io import header_text, run_sl2
from skf_io import skf_header, run_skf

if os.path.exists('runs/ex_env.json'):      # optional overrides of EX_* / PROD_* variables, read at every start (lets a running rotation change its parameters)
    try: os.environ.update({k: str(v) for k, v in json.load(open('runs/ex_env.json')).items()})
    except Exception: pass
name = sys.argv[1]; secs = float(sys.argv[2]); Kmods = json.loads(sys.argv[3]); q = int(sys.argv[4]); hu = json.loads(sys.argv[5]); slotfile = sys.argv[6]; topn = int(sys.argv[7])
order = sys.argv[8] if len(sys.argv) > 8 else 'rho'; engine = sys.argv[9] if len(sys.argv) > 9 else 'skf'
COV = float(os.environ.get('EX_COV', '0.9')); MAXDRAW = int(float(os.environ.get('EX_MAXDRAW', '600000'))); EPOCH = int(os.environ.get('EX_EPOCH', '0'))
CACHE_MB = float(os.environ.get('EX_CACHE_MB', '500')); L2CAP = int(os.environ.get('EX_L2CAP', '400')); MAXA = int(float(os.environ.get('PROD_MAXA', '3e6')))
BATCH = int(float(os.environ.get('EX_BATCH', '30000'))); SLOT_SECS = float(os.environ.get('EX_SLOT_SECS', '900')); UNITS = int(os.environ.get('EX_UNITS', '1'))
DEBUG = bool(os.environ.get('EX_DEBUG')); DUMP = os.environ.get('EX_DUMP')
PERSIST_MB = float(os.environ.get('EX_PERSIST_MB', '3')); JOINCAP = float(os.environ.get('EX_JOINCAP', '2e8'))
os.makedirs('runs', exist_ok=True)
sh = SShape(Kmods, q, hu); l1 = L1(sh); n = sh.N
hdr_sl2 = header_text(sh); hdr_skf = skf_header(sh) if (len(Kmods) == 1 and sh.N <= 64) else None
if hdr_skf is None: engine = 'sl2'
shapekey = re.sub(r'[^0-9A-Za-z]+', '_', 'K' + json.dumps(Kmods) + 'q' + str(q) + 'h' + json.dumps(hu)).strip('_')
LISTDIR = os.path.join('runs', 'lists', shapekey); os.makedirs(LISTDIR, exist_ok=True)

_I = ctypes.c_int; _V = ctypes.c_void_p; _L = ctypes.c_long; _U = ctypes.c_uint64
_lib.gl1_genlist.restype = _L
_lib.gl1_genlist.argtypes = [_I, _I, _I, _V, _V, _V, _V, _V, _I, _V, _I, _I, _I, _U, _U, _L, _I, _I, _V, _V]
_lib.gl1_hashlist.restype = None
_lib.gl1_hashlist.argtypes = [_I, _I, _I, _V, _U, _V, _V, _L, _V]
_lib.gl1_set_hin.restype = None
_lib.gl1_set_hin.argtypes = [_I, _V]
HSEED = 7
def hashlist(S, E):
    S = np.ascontiguousarray(S); E = np.ascontiguousarray(E); out = np.zeros(len(S), np.uint64)
    if len(S): _lib.gl1_hashlist(n, l1.r, len(sh.repsPM), l1.shift.ctypes.data, HSEED, S.ctypes.data, E.ctypes.data, len(S), out.ctypes.data)
    return out
def genlist(nw, mw, seed, j0, cnt, fs, fe):
    S = np.zeros((cnt, n), np.int8); E = np.zeros((cnt, n), np.int8)
    m = _lib.gl1_genlist(n, l1.r, l1.nu, l1.w.ctypes.data, l1.forced.ctypes.data, l1.poff.ctypes.data, l1.plist.ctypes.data, l1.dpw.ctypes.data, l1.nfree, l1.freeidx.ctypes.data, l1.fw,
                         nw, mw, seed, j0, cnt, -100000 if fs is None else fs, -100000 if fe is None else fe, S.ctypes.data, E.ctypes.data)
    return S[:m], E[:m]

# ---- symmetries: identical to prod_dedup.py (exact symmetries of the K-level system)
_Wh = np.random.default_rng(12345).integers(1, 2 ** 62, size=n, dtype=np.uint64) * np.uint64(2) + np.uint64(1)
_fix = [int(P[0]) for P in sh.kCD.pos if len(P) == 1]
_neg = sh.negperm
_permsB = [np.asarray(sh.shift[t], dtype=np.int64)[_neg] for t in _fix] + [np.asarray(sh.shift[t], dtype=np.int64) for t in _fix]
_permsA = [np.arange(n), np.asarray(_neg, dtype=np.int64)]
_MUL = [np.uint64(0x9E3779B97F4A7C15), np.uint64(0xC2B2AE3D27D4EB4F), np.uint64(0x165667B19E3779F9), np.uint64(0x27D4EB2F165667C5)]
# units of K = Z_n (act on all blocks simultaneously; they commute with the multipliers of H_K): A is replaced by an orbit representative
_permsU = None
if UNITS and len(Kmods) == 1:
    _permsU = [(u * np.arange(n)) % n for u in range(1, n) if math.gcd(u, n) == 1]
    assert any((p == np.asarray(_neg)).all() for p in _permsU)
def _canon(Sx, Ex, perms, negs):
    best = None
    for sg in ((1, -1) if negs else (1,)):
        v = ((sg * Sx.astype(np.int64) + 2) * 3 + (sg * Ex.astype(np.int64) + 1) + 1).astype(np.uint64)
        for P in perms:
            h = (v[:, P] * _Wh).sum(axis=1, dtype=np.uint64)
            best = h if best is None else np.minimum(best, h)
    return best
def _rawkey(Sx, Ex):
    v = ((Sx.astype(np.int64) + 2) * 3 + (Ex.astype(np.int64) + 1) + 1).astype(np.uint64)
    return (v * _Wh).sum(axis=1, dtype=np.uint64)
def inst_keys(S, E):
    key = np.zeros(len(S), dtype=np.uint64)
    for X in range(4):
        key = key + _canon(S[:, X, :], E[:, X, :], _permsA if X == 0 else _permsB, X > 0) * _MUL[X]
    return key

# ---- reduced block lists, cached (memory LRU + disk)
class ListCache:
    def __init__(self, mb): self.d = collections.OrderedDict(); self.bytes = 0; self.cap = mb * 1e6; self.hits = 0; self.miss = 0; self.disk = 0
    def get(self, key):
        v = self.d.get(key)
        if v is not None: self.d.move_to_end(key); self.hits += 1
        return v
    def put(self, key, v):
        self.d[key] = v; self.bytes += v[0].nbytes * 2
        while self.bytes > self.cap and len(self.d) > 1:
            k, old = self.d.popitem(last=False); self.bytes -= old[0].nbytes * 2
cache = ListCache(CACHE_MB)
def block_list(nw, mw, fs, fe):
    """(S, E, complete, ndraw): one representative per symmetry orbit of the block configs of defect weight nw (m-weight mw) and row sums (fs, fe)"""
    t_bl = time.time(); key = (nw, mw, fs, fe, EPOCH)
    v = cache.get(key)
    if v is not None: return v
    path = os.path.join(LISTDIR, f'{nw}_{mw}_{fs}_{fe}_{EPOCH}.npz')
    if os.path.exists(path):
        try:
            z = np.load(path); v = (z['S'], z['E'], bool(z['complete']), int(z['draws']), hashlist(z['S'], z['E'])); cache.put(key, v); cache.disk += 1; return v
        except Exception: pass
    cache.miss += 1
    seed = ((((nw * 1009 + mw) * 100003 + (fs + 5000)) * 100003 + (fe + 5000)) * 1009 + EPOCH) * 7919 + 17
    rawk = np.zeros(0, np.uint64); orbk = np.zeros(0, np.uint64); Ss = []; Es = []; draws = 0; j0 = 0; batch = 4096; complete = False
    while draws < MAXDRAW:
        S, E = genlist(nw, mw, seed, j0, batch, fs, fe); j0 += batch; draws += batch
        if len(S) == 0:
            if draws >= 3 * 4096: complete = True; break
            continue
        rk = _rawkey(S, E); urk, ri = np.unique(rk, return_index=True)
        old = np.isin(urk, rawk); nold = int(old.sum()); rawk = np.concatenate([rawk, urk[~old]])
        frac_old = (len(S) - (len(urk) - nold)) / len(S)        # fraction of this batch's draws whose raw config was already present (estimates the covered fraction of the space)
        ok = ~old
        if ok.any():
            ck = _canon(S[ri[ok]], E[ri[ok]], _permsB, True); uck, ci = np.unique(ck, return_index=True)
            newo = ~np.isin(uck, orbk)
            if newo.any():
                orbk = np.concatenate([orbk, uck[newo]]); sel = ri[ok][ci[newo]]; Ss.append(S[sel]); Es.append(E[sel])
        if draws >= 2 * 4096 and frac_old >= COV: complete = True; break
        batch = min(batch * 2, 131072)
    if Ss: SS = np.concatenate(Ss); EE = np.concatenate(Es)
    else: SS = np.zeros((0, n), np.int8); EE = np.zeros((0, n), np.int8)
    v = (SS, EE, complete, draws, hashlist(SS, EE)); cache.put(key, v)
    try:
        if SS.nbytes * 2 <= PERSIST_MB * 1e6:
            tmp = path + f'.{os.getpid()}.tmp.npz'; np.savez(tmp, S=SS, E=EE, complete=complete, draws=draws); os.replace(tmp, path)
    except Exception: pass
    if DEBUG: print(f'   list nw={nw} mw={mw} fs={fs} fe={fe}: draws {draws} raw distinct {len(rawk)} orbits {len(SS)} complete {complete} t {time.time() - t_bl:.2f}s', flush=True)
    return v

# ---- per-slot data (A list, patterns)
def slot_data(al, ms):
    try:
        pats = row_patterns_skew(sh, al, ms)
        A_list = sh.enum_A(al[0], ms[0], cap=MAXA) if pats else None
    except RuntimeError as ex:      # 'A list too large' (more than PROD_MAXA skew blocks): the slot is skipped and logged as skip = A_toobig
        return 'A_toobig'
    if not pats or A_list is None: return None
    SA, EA = A_list; nA0 = len(SA)
    ka = _canon(SA, EA, _permsU if _permsU is not None else _permsA, False); _, ia = np.unique(ka, return_index=True)
    SA = np.ascontiguousarray(SA[ia]); EA = np.ascontiguousarray(EA[ia])
    # patterns related by exchanging two blocks of equal type (weight, m-weight) are equivalent (the GS array allows permuting B, C, D): keep one of each
    def keep(p):
        for i in (0, 1):
            if (al[i + 1], ms[i + 1]) == (al[i + 2], ms[i + 2]) and tuple(p[i]) > tuple(p[i + 1]): return False
        return True
    pats = sorted([tuple(tuple(int(x) for x in b) for b in p) for p in pats if keep(p)])
    return dict(SA=SA, EA=EA, HA=hashlist(SA, EA), nA0=nA0, pats=pats)

def join_pattern(al, ms, sd, p, stats):
    """level-1 instances of one pattern: reduced explicit lists, cheapest pairing; returns (S4, E4) in block order A, B, C, D (or None)"""
    t1 = time.time(); lists = []
    for b in range(3):
        fs, fe = p[b]; S, E, comp, dr, H = block_list(al[b + 1], ms[b + 1], fs, fe)
        stats['allcomplete'] &= comp
        if len(S) == 0: stats['tg'] += time.time() - t1; return None
        lists.append((S, E, H))
    stats['tg'] += time.time() - t1; t1 = time.time()
    blk = [(sd['SA'], sd['EA'], sd['HA'])] + lists; sz = [len(b[0]) for b in blk]
    # the join costs |X0||X2| + |X1||X3| (specs 0,2 hashed together, 1,3 probed): choose the cheapest pairing of {A, B, C, D}
    pr = min((((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2))), key=lambda pr_: sz[pr_[0][0]] * sz[pr_[0][1]] + sz[pr_[1][0]] * sz[pr_[1][1]])
    jc = sz[pr[0][0]] * sz[pr[0][1]] + sz[pr[1][0]] * sz[pr[1][1]]
    if jc > JOINCAP: stats['nskip'] += 1; return None      # join too large for this pass (reported as nskip in the slot log; the slot is NOT exhaustive then)
    ordr = [pr[0][0], pr[1][0], pr[0][1], pr[1][1]]; inv = np.argsort(ordr)
    specs = [BlockSpec(al[k], S=blk[k][0], E=blk[k][1], cansign=False) for k in ordr]
    mo = 20000
    for i_, k_ in enumerate(ordr): _lib.gl1_set_hin(i_, blk[k_][2].ctypes.data)
    try:
        cnt, S4, E4, st = l1.run(specs, per_class=30000.0, maxout=mo, hseed=HSEED)
        if cnt > mo:
            stats['over'] += 1; cnt, S4, E4, st = l1.run(specs, per_class=30000.0, maxout=int(cnt) + 10, hseed=HSEED)
    finally:
        for i_ in range(4): _lib.gl1_set_hin(i_, None)
    stats['tj'] += time.time() - t1
    if DEBUG: print(f'  pat {p} sizes {sz} pairing {pr} cnt {cnt} join {time.time() - t1:.3f}s', flush=True)
    if not cnt: return None
    stats['nprod'] += 1; stats['nraw'] += int(cnt)
    return S4[:, inv].copy(), E4[:, inv].copy()

# ---- slots
slots = []
for l in open(slotfile):
    r = json.loads(l)
    if 'logE' in r and 'nopat' not in r and r['logE'] > -50: slots.append((r['logE'], r.get('logrho', 0), r.get('logN1', 0), tuple(r['al']), tuple(r['ms'])))
minrho = float(os.environ.get('EX_MINRHO', '-99')); mine = float(os.environ.get('EX_MINE', '-99'))
slots = [s for s in slots if s[1] >= minrho and s[0] >= mine]
MAXCLASS = float(os.environ.get('EX_MAXCLASS', '1e30'))      # skip slots whose estimated number of classes N1 / (2 (4 |Fix|)^3) exceeds this
slots = [s for s in slots if 10 ** s[2] / (2 * (4 * max(1, len(_fix))) ** 3) <= MAXCLASS]
if order == 'rho': slots.sort(key=lambda s: (-s[1], -s[0]))
elif order == 'E': slots.sort(key=lambda s: -s[0])
elif order == 'mix': slots.sort(key=lambda s: -(s[0] + s[1]))       # logE + logrho = 2 logE - logN1
slots = slots[:topn]
if os.environ.get('EX_SLOTS'):
    a, b = os.environ['EX_SLOTS'].split('-'); slots = slots[int(a):int(b)]
if os.environ.get('EX_SHARD'):
    i_, n_ = [int(x) for x in os.environ['EX_SHARD'].split('/')]; slots = [s for k, s in enumerate(slots) if k % n_ == i_]
slog = f'runs/{name}.slots.jsonl'; state = {}      # slot -> next pattern index (epoch EPOCH); finished slots have next >= 10**9
if os.path.exists(slog):
    for l in open(slog):
        try: r = json.loads(l)
        except Exception: continue
        if r.get('epoch', 0) == EPOCH:
            k = (tuple(r['al']), tuple(r['ms'])); state[k] = max(state.get(k, 0), r.get('pat_end', 0)) if not r.get('finished') else 10 ** 9
seen_path = f'runs/{name}.seen'
seen = set(np.fromfile(seen_path, dtype=np.uint64).tolist()) if os.path.exists(seen_path) else set()
print(f'{name}: {len(slots)} slots ordered by {order}, {sum(1 for v in state.values() if v >= 10 ** 9)} finished at epoch {EPOCH}, seen {len(seen)}, engine {engine}, N {n} nu {sh.kCD.nu}, units {bool(_permsU)}', flush=True)
tot = dict(visits=0, raw=0, new=0, l2=0, nok=0, tg=0.0, tj=0.0, tk=0.0, t2=0.0, hits=0)
stat_all = collections.Counter(); t_start = time.time(); nbatch_dump = 0

def run_level2(S, E, si, al, ms, rec, stat):
    """new classes -> exact level 2; returns (number completing, number run)"""
    nok = 0; ndone = 0
    for c0 in range(0, len(S), L2CAP):
        insts = [[(S[i, b], E[i, b]) for b in range(4)] for i in range(c0, min(len(S), c0 + L2CAP))]
        if engine == 'skf': res, ss = run_skf(sh, insts, maxsol=3, tlimit=60.0, hdr=hdr_skf)
        else: res, ss = run_sl2(sh, insts, maxsol=3, tlimit=60.0, maxfib=1500000, stop_first=True, hdr=hdr_sl2)
        ndone += len(insts)
        for i, d in res.items():
            stat[d['status']] = stat.get(d['status'], 0) + 1; nok += d['compl'] > 0
        for (i, ts) in ss[:3]:
            sol = [(insts[i][X][0].astype(np.int64), ts[X].astype(np.int64), insts[i][X][1].astype(np.int64)) for X in range(4)]
            assert sh.check_solution(sol), 'K-level check failed'
            from sk_cert import certify
            tag = f'{name}_s{si}_{int(time.time())}_{c0 + i}'
            res_c = certify(sh, sol, tag=tag); tot['hits'] += 1; rec.setdefault('HIT', []).append(res_c)
            json.dump(dict(blocks=[[x.tolist() for x in b] for b in sol], al=al, ms=ms, cert=res_c, Kmods=Kmods, q=q, hu=hu, engine='exh38'), open(f'runs/{name}.HIT_{tag}.json', 'w'))
            print('HIT', res_c, 'slot', al, ms, flush=True)
    return nok, ndone

def flush_batch(Sb, Eb, si, al, ms, rec, stat):
    global nbatch_dump
    if not Sb: return
    S = np.concatenate(Sb); E = np.concatenate(Eb); t_k = time.time()
    if DUMP: np.savez(os.path.join(DUMP, f'dump_{name}_{si}_{nbatch_dump}.npz'), S=S, E=E); nbatch_dump += 1
    keys_ = inst_keys(S, E); uk_, first_ = np.unique(keys_, return_index=True)
    msk_ = np.fromiter((int(k_) not in seen for k_ in uk_), dtype=bool, count=len(uk_))
    sel_ = first_[msk_]; newkeys = uk_[msk_]; S = S[sel_]; E = E[sel_]
    rec['nclass'] = rec.get('nclass', 0) + int(len(uk_)); tot['tk'] += time.time() - t_k
    t_2 = time.time(); nok, ndone = run_level2(S, E, si, al, ms, rec, stat); tot['t2'] += time.time() - t_2
    seen.update(int(k_) for k_ in newkeys[:ndone]); open(seen_path, 'ab').write(newkeys[:ndone].tobytes())
    rec['nnew'] = rec.get('nnew', 0) + len(newkeys); rec['nok'] = rec.get('nok', 0) + int(nok); tot['new'] += len(newkeys); tot['l2'] += ndone; tot['nok'] += nok

def all_finished(): return all(state.get((s[3], s[4]), 0) >= 10 ** 9 for s in slots)
passno = 0
while time.time() - t_start < secs and not os.path.exists(f'runs/{name}.STOP') and not all_finished():
    passno += 1
    for si, (lE, lrho, lN1, al, ms) in enumerate(slots):
        if time.time() - t_start > secs or os.path.exists(f'runs/{name}.STOP'): break
        k = (al, ms)
        if state.get(k, 0) >= 10 ** 9: continue
        t_slot = time.time()
        sd = slot_data(al, ms)
        rec = dict(al=list(al), ms=list(ms), epoch=EPOCH, logE=lE, logrho=lrho, logN1=lN1)
        if sd is None or sd == 'A_toobig':
            rec.update(skip=(1 if sd is None else 'A_toobig'), finished=1, pat_end=0); open(slog, 'a').write(json.dumps(rec) + '\n'); state[k] = 10 ** 9; continue
        pats = sd['pats']; p0 = state.get(k, 0); rec.update(npat=len(pats), nA=len(sd['SA']), nA0=sd['nA0'], pat_start=p0)
        stats = dict(tg=0.0, tj=0.0, nprod=0, nraw=0, over=0, nskip=0, allcomplete=True); stat = {}; Sb = []; Eb = []; nb = 0; pi = p0
        while pi < len(pats):
            r_ = join_pattern(al, ms, sd, pats[pi], stats); pi += 1
            if r_ is not None: Sb.append(r_[0]); Eb.append(r_[1]); nb += len(r_[0])
            if nb >= BATCH: flush_batch(Sb, Eb, si, al, ms, rec, stat); Sb = []; Eb = []; nb = 0
            if time.time() - t_slot > SLOT_SECS or time.time() - t_start > secs: break
        flush_batch(Sb, Eb, si, al, ms, rec, stat)
        fin = pi >= len(pats)
        rec.update(pat_end=pi, finished=int(fin), nprod=stats['nprod'], nraw=stats['nraw'], t_gen=round(stats['tg'], 2), t_join=round(stats['tj'], 2), complete=bool(stats['allcomplete']), over=stats['over'], nskip=stats['nskip'], stat=stat, t_slot=round(time.time() - t_slot, 2))
        open(slog, 'a').write(json.dumps(rec) + '\n'); state[k] = 10 ** 9 if fin else pi
        tot['visits'] += 1; tot['raw'] += stats['nraw']; tot['tg'] += stats['tg']; tot['tj'] += stats['tj']; stat_all.update(stat)
        if tot['visits'] % 5 == 0 or DEBUG:
            print(f'{name}: pass {passno} visits {tot["visits"]} (slot idx {si}) raw {tot["raw"]} new classes {tot["new"]} l2 {tot["l2"]} completions {tot["nok"]} hits {tot["hits"]} t_gen {tot["tg"]:.0f} t_join {tot["tj"]:.0f} t_key {tot["tk"]:.0f} t_l2 {tot["t2"]:.0f} cache {len(cache.d)}/{cache.bytes / 1e6:.0f}MB hit/disk/miss {cache.hits}/{cache.disk}/{cache.miss} stat {dict(stat_all)} wall {time.time() - t_start:.0f}', flush=True)
print(f'{name}: done: ' + json.dumps(tot) + ' stat ' + json.dumps(dict(stat_all)) + f' wall {time.time() - t_start:.0f}' + (' ALL SLOTS FINISHED' if all_finished() else ''), flush=True)
