"""count, for every subgroup H of Z_v^* with Omin <= #orbits <= Omax, whether a row sum +1 is attainable (rho = sum size*x = 1) -- by a subset-sum DP, no enumeration"""
import sys, numpy as np
import os
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import plain_spence3 as ps
def rho_feasible(size):
    # exact: number of x in {+-1}^O with sum size*x = 1 via DP on sums
    tot = int(sum(size)); off = tot
    dp = np.zeros(2 * tot + 1, dtype=np.float64); dp[off] = 1
    for s in size:
        s = int(s)
        new = np.zeros_like(dp)
        new[s:] += dp[:-s] if s > 0 else dp
        new[:-s] += dp[s:]
        dp = new
    return dp[off + 1]
v = int(sys.argv[1]); Omin = int(sys.argv[2]); Omax = int(sys.argv[3])
subs = ps.all_subgroups(v)
rows = []
for H in subs:
    lab, reps, O, size, cls = ps.structure(v, H)
    if not (Omin <= O <= Omax): continue
    L = rho_feasible(size)
    rows.append((O, len(H), len(cls), L))
rows.sort()
print('v', v, 'subgroups', len(subs), 'shapes in range', len(rows), 'with L>0:', sum(1 for r in rows if r[3] > 0))
for r in rows:
    if r[3] > 0: print('  O %d |H| %d classes %d  L %.3g' % r)
