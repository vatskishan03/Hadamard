"""Estimates of the number of normalised solutions from the logs of tri2 (search/NOTES.md section 92; HEURISTIC when the log covers only some of the M partitions, EXACT when it covers all of them).
A partition with a random residue holds N/M solutions on average (N = the normalised count: the weights w of the reduced hits of red=1 sum to N exactly as the unreduced hits do), so  E = (sum of the weights) * M / (number of partitions).
The relative error is that of a Poisson count of independent clusters: sqrt(sum w_p^2) / sum w_p over the partitions p with hits (a cluster of 32 solutions, PROOFS.md P5, is one event).
usage: python3 est_logs.py M FILE [FILE ...]                  one line per log
       python3 est_logs.py --combine M FILE_UNREDUCED FILE_REDUCED     the two logs run on the same sequence of partitions (same r0, step): the unreduced log for its partitions, the reduced log for the others (no partition counted twice)"""
import sys, math
def read(fn):
    st = []
    for l in open(fn):
        if l.startswith('STAT'):
            f = dict(x.split('=') for x in l.split()[1:]); st.append((int(f['i']), int(f['hits']), float(f['wsum'])))
    return st
def report(name, M, st, lo=0, hi=None):
    st = [s for s in st if s[0] >= lo and (hi is None or s[0] < hi)]; n = len(st); W = sum(s[2] for s in st); W2 = sum(s[2] ** 2 for s in st)
    ev = sum(1 for s in st if s[2] > 0)
    E = W * M / n if n else float('nan'); rel = math.sqrt(W2) / W if W else float('nan')
    if n == M and len(set(x[0] for x in st)) == M: print('%-34s COMPLETE scan, all %d partitions: hit partitions %d, exact normalised count (sum of the weights) = %.0f' % (name, M, ev, W)); return n, W, W2
    print('%-34s partitions %6d of M = %d  hit partitions %3d  weight %8.0f   E = %.3g +- %.0f%%' % (name, n, M, ev, W, E, 100 * rel))
    return n, W, W2
if sys.argv[1] == '--combine':
    M = int(sys.argv[2]); u = read(sys.argv[3]); r = read(sys.argv[4]); nu = len(u)
    n1, W1, V1 = report(sys.argv[3] + ' (unreduced)', M, u)
    n2, W2_, V2 = report(sys.argv[4] + ' (all)', M, r)
    n3, W3, V3 = report(sys.argv[4] + ' (i >= %d)' % nu, M, r, lo=nu)
    n = n1 + n3; W = W1 + W3; E = W * M / n; rel = math.sqrt(V1 + V3) / W
    print('combined: unreduced partitions 0..%d + reduced partitions %d..%d: %d partitions, weight %.0f, E = %.3g +- %.0f%%' % (nu - 1, nu, len(r) - 1, n, W, E, 100 * rel))
else:
    M = int(sys.argv[1])
    for fn in sys.argv[2:]: report(fn, M, read(fn))
