"""Tabulate the logs of the exhaustive propus searches (logs/*.log, logs/*.out) into one line per (structure, row pattern).

    python3 summarize_logs.py logs > results.tsv        # per row pattern
    python3 summarize_logs.py logs --structures         # one line per structure (the table of README.md)

Reads only what the logs print: the header line of p3 (v, C, O_A, O_B, row sums, threads, canonical-form switches), COUNT lines
(leaves and kept members of a list), the LISTS line (A', B', D'), and the DONE line (hits, matches, badmatches, pairs, join
and total seconds).  The multiplier group of a structure is read from its job name (sp_v<v>_H<|H|>_g<g1>x<g2>.. = H generated
by g1, g2, ..), and |H| is recomputed here as a check; the two special structures (v = 239, run by hand, and sym343_S38, run
by p3g_job.py over the group Z_7^3) are named in SPECIAL below.  The log of sym343_S38 is not in p3's format: its lines give
the list sizes and the DONE dictionary of each row pattern.  Nothing in the logs is estimated."""
import glob, os, re, sys
from math import gcd

SPECIAL = {                                     # log stem -> (v, generators as text, |H| stated in the job's own log or name)
    'v239_run239_0_1': (239, '10', 7),
    'v239_run239_noreduce': (239, '10', 7),
}


def subgroup_order(v, gens):
    H = {1}
    while True:
        grown = H | {h * g % v for h in H for g in gens}
        if grown == H:
            return len(H)
        H = grown


def parse_p3(text, name):
    out, cur = [], None
    for line in text.splitlines():
        m = re.match(r"p3: v=(\d+) C=(\d+) OA=(\d+) OB=(\d+) rows=\((\d+),(-?\d+),(-?\d+)\) NS=(\d+) W=(\d+) threads=(\d+) "
                     r"shard (\d+)/(\d+) canonB=(\d)\(perms (\d+)\) canonD=(\d)", line)
        if m:
            cur = dict(name=name, v=int(m[1]), C=int(m[2]), OA=int(m[3]), OB=int(m[4]), rows=(int(m[5]), int(m[6]), int(m[7])),
                       NS=int(m[8]), threads=int(m[10]), canon=int(m[13]), mult=int(m[14]) + 1)
            out.append(cur)
            continue
        if cur is None:
            continue
        m = re.match(r"LISTS (\d+) (\d+) (\d+)", line)
        if m:
            cur['lists'] = tuple(int(x) for x in m.groups())
        m = re.match(r"COUNT (\w) (\d+) (\d+)", line)
        if m:
            cur['leaves_' + m[1]] = int(m[2])
        m = re.match(r"DONE hits=(\d+) matches=(\d+) badmatches=(\d+) pairs=(\d+) join_seconds=([\d.]+) total_seconds=([\d.]+)", line)
        if m:
            cur.update(hits=int(m[1]), matches=int(m[2]), bad=int(m[3]), pairs=int(m[4]), join_s=float(m[5]), total_s=float(m[6]))
    return out


def parse_sym343(text, name):
    out = []
    hdr = re.match(r"(\S+): v=(\d+) \|H\|=(\d+) OA=(\d+) OB=(\d+) C=(\d+) \|mult\|=(\d+)", text)
    for m in re.finditer(r"rows \((\d+), (\d+), (\d+)\) shard 0/1: lists\(A,B,D\) \((\d+), (\d+), (\d+)\) done (\{[^}]*\})", text):
        d = eval(m[7])
        out.append(dict(name=name, v=int(hdr[2]), C=int(hdr[6]), OA=int(hdr[4]), OB=int(hdr[5]), rows=tuple(int(m[i]) for i in (1, 2, 3)),
                        NS=3, threads=8, canon=1, mult=int(hdr[7]), lists=tuple(int(m[i]) for i in (4, 5, 6)),
                        hits=int(d['hits']), matches=int(d['matches']), bad=int(d['badmatches']), pairs=int(d['pairs']),
                        join_s=float(d['join_seconds']), total_s=float(d['total_seconds']), Hord=int(hdr[3])))
    return out


def gens_of(name):
    m = re.match(r"sp_v(\d+)_H(\d+)_g([\dx]+)", name)
    if m:
        return int(m[1]), [int(x) for x in m[3].split('x')], int(m[2])
    return None


def load(logdir):
    rows = []
    for path in sorted(glob.glob(os.path.join(logdir, '*'))):
        stem = os.path.basename(path).rsplit('.', 1)[0]
        text = open(path).read()
        if stem.startswith('sym343_S38_run_02'):
            rs = parse_sym343(text, 'sym343_S38')
            for r in rs:
                r['gens'] = 'S^38 (Singer cycle of GF(7^3), order 342)'
        elif stem.startswith('sp_v'):
            v, gens, h = gens_of(stem.replace('_ns3', ''))
            rs = parse_p3(text, re.sub(r'_ns3$', '', stem))
            for r in rs:
                r['gens'] = ','.join(map(str, gens)); r['Hord'] = subgroup_order(v, gens)
                assert r['Hord'] == h, (stem, r['Hord'], h)
                assert r['v'] == v
        elif stem in SPECIAL:
            v, g, h = SPECIAL[stem]
            rs = parse_p3(text, stem)
            for r in rs:
                r['gens'] = g; r['Hord'] = subgroup_order(v, [int(g)]); assert r['Hord'] == h
        else:
            continue
        rows += rs
    return rows


def main():
    logdir = sys.argv[1]
    rows = load(logdir)
    if '--structures' in sys.argv:
        print("name\tv\tH\t|H|\tC\tO_A\tO_B\tpatterns\tpairs\tp3_seconds\tthreads\thits")
        by = {}
        for r in rows:
            by.setdefault(r['name'], []).append(r)
        for k, rs in by.items():
            r0 = rs[0]
            print("\t".join(map(str, (k, r0['v'], r0['gens'], r0['Hord'], r0['C'], r0['OA'], r0['OB'], len(rs), sum(r['pairs'] for r in rs),
                                      round(sum(r['total_s'] for r in rs), 1), r0['threads'], sum(r['hits'] for r in rs)))))
        return
    print("name\tv\tH\t|H|\trows\tC\tO_A\tO_B\tpool_A\tpool_B\tpool_D\tA'\tB'\tD'\tpairs\tjoin_s\ttotal_s\tthreads\thits\tmatches\tbadmatches\tcanonical_reductions")
    for r in rows:
        pools = [r.get('leaves_' + x, '-') for x in 'ABD']
        print("\t".join(map(str, (r['name'], r['v'], r['gens'], r['Hord'], ','.join(map(str, r['rows'])), r['C'], r['OA'], r['OB'], *pools, *r['lists'],
                                  r['pairs'], r['join_s'], r['total_s'], r['threads'], r['hits'], r['matches'], r['bad'], r['canon']))))


if __name__ == '__main__':
    main()
