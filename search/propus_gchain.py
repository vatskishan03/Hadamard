"""Bordered propus over non-cyclic abelian groups with compressions to a quotient group (larger structures than
propus_bgroup.py joins directly).

For a characteristic subgroup K of G (K = kG or the k-torsion G[k]), the compressions x'(c) = sum_{g in c} X(g)
of the three blocks to Q = G/K satisfy PAF_a' + 2 PAF_b' + PAF_d' = -4 |K| at every c != 0, and are constant on the
orbits that H (and <H, -1> for A) induces on Q. Stage 1 lists the compressions each block attains (dynamic program
over the orbits, propus_chain.Comp) and joins them exactly; stage 2 lists the members of each compatible triple by
meet in the middle (propus_chain.members) and joins them on all autocorrelation classes with join4x. Every hit is
rebuilt as the bordered matrix and checked (H = H^T, H H^T = 4(m + 1) I).

  python3 search/propus_gchain.py 8,8 [--mino 31] [--maxo 46] [--pairs] [--maxwork 3e10] [--maxlist 3e7] [--maxend 2e6]
Hits go to $HADAMARD_WORK/sym_hits.jsonl (border: true).
"""
import argparse, json, math, os, struct, subprocess, sys, tempfile, time
import numpy as np
from pathlib import Path as _P
sys.path[:0] = [str(_P(__file__).resolve().parent), str(_P(__file__).resolve().parent.parent / 'verify')]
import workdir as W
import spence_lib as L
import propus_chain as PC
import propus_group as PG
import propus_bgroup as BG


def quotients(G):
    """(name, Q, q) for the characteristic subgroups K = kG and G[k]: q[g] is the index in Q = G/K."""
    n = [int(x) for x in G.mods]; e = math.lcm(*n); out = []
    for k in range(2, e + 1):
        if e % k: continue
        mods = [math.gcd(x, k) for x in n]                                  # G / kG
        if 1 < math.prod(mods) < G.v:
            Q = L.AbGroup(mods); out.append((f"G/{k}G", Q, Q.index(G.coords % np.array(mods))))
        mods2 = [x // math.gcd(x, k) for x in n]                            # G / G[k] = kG
        if 1 < math.prod(mods2) < G.v:
            Q = L.AbGroup(mods2)
            img = (G.coords * k) // np.array([math.gcd(x, k) for x in n])
            out.append((f"G/G[{k}]", Q, Q.index(img % np.array(mods2))))
    return out


def induced(perms, q, nq):
    """Permutations of Q induced by permutations of G that preserve the cosets of K."""
    first = np.full(nq, -1, np.int64)
    for g in range(len(q) - 1, -1, -1): first[q[g]] = g
    return [q[p[first]] for p in perms]


def make_comp(st, q, Q, qperms):
    lab, k = L.closure_partition(Q.v, qperms)
    reps = [int(np.argmax(lab == o)) for o in range(k)]
    c = PC.Comp.__new__(PC.Comp)
    c.d = Q.v; c.roid = lab; c.rreps = reps; c.k = k
    cnt = np.zeros((st.O, Q.v), np.int64)
    np.add.at(cnt, (st.orb, q), 1)
    c.cnt = cnt[:, reps]
    return c


def qpaf(Q, x):
    return np.array([int(np.dot(x, x[Q.translate(c)])) for c in range(1, Q.v)], np.int64)


def search(G, name, mats, gens, st_ad, st_b, a, log):
    m = G.v; pm = m % 4; found = []
    tgt, pats, build, n4 = BG.mode(G)                 # tgt: weighted PAF target at each g != 0
    neg = G.neg()
    best = None
    for qname, Q, q in quotients(G):
        kb = len(set(L.closure_partition(Q.v, induced(gens, q, Q.v))[0].tolist()))
        if kb <= 10 and (best is None or Q.v > best[1].v): best = (qname, Q, q)
    if best is None:
        log(f"m={m} {name}: no quotient with at most 10 orbits"); return []
    qname, Q, q = best
    ca = make_comp(st_ad, q, Q, induced(gens + [neg], q, Q.v)); cb = make_comp(st_b, q, Q, induced(gens, q, Q.v))
    for ra, rb, rd in pats:
        if min(st_ad.pool(ra), st_b.pool(rb), st_b.pool(rd)) == 0: continue
        t0 = time.time()
        Da = ca.dp(st_ad, ra); Db = cb.dp(st_b, rb); Dd = cb.dp(st_b, rd)
        if Da is None or Db is None or Dd is None:
            log(f"m={m} {name} {qname} rows ({ra}, {rb}, {rd}): compression DP above its cap; skipped"); continue
        KA, KB, KD = list(Da), list(Db), list(Dd)
        PA = np.array([qpaf(Q, ca.full(k)) for k in KA]).reshape(len(KA), -1)
        PB = np.array([2 * qpaf(Q, cb.full(k)) for k in KB]).reshape(len(KB), -1)
        PD = np.array([qpaf(Q, cb.full(k)) for k in KD]).reshape(len(KD), -1)
        triples = []
        tgtQ = np.array([int(tgt[(q == c) & (np.arange(m) != 0)].sum()) for c in range(1, Q.v)], np.int64)
        for i, j, k in PC.join3(PA, PB, PD, tgtQ):
            triples.append((Da[KA[i]] * Dd[KD[k]], KA[i], KB[j], KD[k]))
        triples.sort(key=lambda x: -x[0])
        work = sum(t[0] for t in triples)
        log(f"m={m} {name} {qname} (|Q| {Q.v}) rows ({ra}, {rb}, {rd}): orbits A {st_ad.O}, B/D {st_b.O}, classes {st_ad.C}; "
            f"compressions {len(Da)}, {len(Db)}, {len(Dd)}; {len(triples)} compatible triples, pair work {work:.2e} [{time.time() - t0:.1f}s]")
        if max(st_ad.O, st_b.O) > 46: log("  too many orbits for the meet in the middle"); continue
        done = 0.0
        for i, (w_, ka_, kb_, kd_) in enumerate(triples):
            if done + w_ > a.maxwork: log(f"  work cap reached after {i} triples"); break
            Ya = PC.members(st_ad, ca, ra, ka_, a.maxlist); Yb = PC.members(st_b, cb, rb, kb_, a.maxlist)
            Yd = PC.members(st_b, cb, rd, kd_, a.maxlist)
            if Ya is None or Yb is None or Yd is None: log(f"  triple {i}: list above cap; skipped"); continue
            def vals(st, Y, s, chunk=1 << 18):
                out = np.empty((len(Y), st.C), dtype=np.int16)
                for c0 in range(0, len(Y), chunk):
                    out[c0:c0 + chunk] = s * np.rint((st.class_values(Y[c0:c0 + chunk]) - pm) / 4)
                return out
            Va, Vb, Vd = vals(st_ad, Ya, 1), vals(st_b, Yb, 2), vals(st_b, Yd, 1)
            Z = np.zeros((1, st_ad.C), np.int16)
            roles = sorted([(len(Va), 0), (len(Vd), 1), (len(Vb), 2)]); byrole = {0: Va, 1: Vd, 2: Vb}
            lists = [byrole[roles[0][1]], byrole[roles[1][1]], byrole[roles[2][1]], Z]
            slot = {r: j for j, (_, r) in enumerate(roles)}
            fd, path = tempfile.mkstemp(suffix='.bin', dir=W.tmpdir()); os.close(fd)
            try:
                with open(path, 'wb') as f:
                    f.write(struct.pack('iiiii', *[len(x) for x in lists], st_ad.C))
                    f.write(np.array([(tgt[g] - 4 * pm) // 4 for g in st_ad.reps], dtype=np.int16).tobytes())
                    for x in lists: f.write(np.ascontiguousarray(x).tobytes())
                out = subprocess.run([W.bin('join4x'), path, '1'], capture_output=True, text=True, check=True).stdout
            finally:
                os.remove(path)
            total = int(out.strip().split('TOTAL')[-1]); done += w_
            if i < 3 or i % 200 == 0 or total:
                log(f"  triple {i}: lists {len(Ya)}, {len(Yb)}, {len(Yd)}: {total} solutions [{time.time() - t0:.0f}s]")
            for line in out.splitlines():
                if not line.startswith('SOL'): continue
                ids = list(map(int, line.split()[1:]))
                A = Ya[ids[slot[0]]].astype(np.int64)[st_ad.orb]; D = Yd[ids[slot[1]]].astype(np.int64)[st_b.orb]
                B = Yb[ids[slot[2]]].astype(np.int64)[st_b.orb]
                tot = PG.group_paf(G, A) + 2 * PG.group_paf(G, B) + PG.group_paf(G, D)
                if not (tot[0] == 4 * m and np.array_equal(tot[1:], tgt[1:])): log('  WRONG autocorrelations'); continue
                Hm, Kb = build(A, B, D)
                ok = Hm is not None and np.array_equal(Hm, Hm.T) and np.array_equal(Hm @ Hm.T, n4 * np.eye(n4, dtype=Hm.dtype))
                log(f"  hit: {'' if BG.UNBORDERED else 'bordered '}symmetric Hadamard of order {n4} {ok}")
                if ok:
                    rec = dict(n=m, order=n4, border=not BG.UNBORDERED, group=[int(x) for x in G.mods], structure=name,
                               automorphisms=[x.tolist() for x in mats], freeD=True, rows=[int(A.sum()), int(B.sum()), int(D.sum())],
                               method=f'propus_gchain {qname}', K=None if Kb is None else Kb.tolist(),
                               seqs=[''.join('+' if x > 0 else '-' for x in s) for s in (A, B, D)])
                    with open(W.hits('sym_hits.jsonl'), 'a') as fh: fh.write(json.dumps(rec) + '\n')
                    found.append(rec)
            if found: return found
        log(f"  done: {min(len(triples), i + 1) if triples else 0} of {len(triples)} triples, work {done:.2e} [{time.time() - t0:.0f}s]")
    return found


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('mods'); ap.add_argument('--mino', type=int, default=31); ap.add_argument('--maxo', type=int, default=46)
    ap.add_argument('--maxwork', type=float, default=3e10); ap.add_argument('--maxlist', type=float, default=3e7)
    ap.add_argument('--maxend', type=float, default=2e6); ap.add_argument('--pairs', action='store_true')
    ap.add_argument('--unbordered', action='store_true', help='plain propus on G (order 4|G|, |G| odd)')
    a = ap.parse_args()
    BG.UNBORDERED = a.unbordered
    log = lambda s: print(s, flush=True)
    G = L.AbGroup([int(x) for x in a.mods.split(',')])
    auts = BG.automorphisms(G, a.maxend)
    if auts is None: log(f"G={a.mods}: too many endomorphisms"); return
    neg = G.neg()
    uniq = {}
    for M, p in auts:
        perms = BG.subgroup_perms(p, G.v)
        if any(np.array_equal(x, neg) for x in perms): continue
        inv = BG.invariant(G, perms)
        if inv not in uniq: uniq[inv] = ([M], [p])
    if a.pairs:
        firsts = list(uniq.values()); cyc = [(M, p) for M, p in auts]
        extra = {}
        for (M1,), (p1,) in [(x[0][:1], x[1][:1]) for x in firsts]:
            for M2, p2 in cyc:
                lab, k = L.closure_partition(G.v, [p1, p2])
                if k < a.mino or k > a.maxo: continue
                grp = {tuple(range(G.v))}; frontier = [np.arange(G.v)]; ok = True
                while frontier and ok:
                    x = frontier.pop()
                    for g in (p1, p2):
                        r = g[x]; t = tuple(r.tolist())
                        if t not in grp:
                            if np.array_equal(r, neg): ok = False; break
                            grp.add(t); frontier.append(r)
                    if len(grp) > 4096: ok = False
                if not ok: continue
                inv = BG.invariant(G, [np.array(t) for t in grp])
                if inv not in uniq and inv not in extra: extra[inv] = ([M1, M2], [p1, p2])
        uniq.update(extra)
    todo = []
    for inv, (mats, gens) in uniq.items():
        st_b = L.Structure(G, gens); st_ad = L.Structure(G, gens + [neg])
        if a.mino <= st_b.O <= a.maxo: todo.append((st_b.O, mats, gens, st_ad, st_b))
    todo.sort(key=lambda x: x[0])
    log(f"G={a.mods} (order {G.v}): {len(auts)} automorphisms, {len(todo)} structures with {a.mino}..{a.maxo} orbits")
    for _, mats, gens, st_ad, st_b in todo:
        name = 'aut' + '_'.join(''.join(str(x) for x in M.flatten()) for M in mats)
        if search(G, name, mats, gens, st_ad, st_b, a, log):
            log(f"HIT m={G.v} G={a.mods} {name}"); return
    log(f"=== G={a.mods}: none")


if __name__ == '__main__':
    main()
