"""Rank-one bicliques K_{d,d} in the Paley two-graph of Z_p: X+ and X- (|X+| = |X-| = d, disjoint) with chi(y - x) = u_x v_y.
0 in X- by translation; then X+ lies in T(y, e) = {x : chi(y - x) = e chi(x)} for every y in X- minus 0. DFS over (y, e)."""
import sys
p = int(sys.argv[1]); d = int(sys.argv[2]) if len(sys.argv) > 2 else (p - 1) // 4
qr = {x * x % p for x in range(1, p)}
chi = lambda x: 0 if x % p == 0 else (1 if x % p in qr else -1)
T = {}
for y in range(1, p):
    for e in (1, -1):
        T[y, e] = sum(1 << x for x in range(1, p) if x != y and chi(y - x) == e * chi(x))
full = sum(1 << x for x in range(1, p))
count = 0; sols = []
def dfs(Y, mask, start):
    global count
    if len(Y) == d - 1:
        cand = mask & ~sum(1 << y for y, e in Y)
        if bin(cand).count('1') >= d:
            count += 1
            if len(sols) < 5: sols.append(([0] + [y for y, e in Y], [x for x in range(p) if cand >> x & 1]))
        return
    for y in range(start, p):
        for e in (1, -1):
            if not (mask >> y & 1) and True:
                pass
            m2 = mask & T[y, e] & ~(1 << y)
            # remaining candidates for X+ must stay >= d after removing future Y members
            if bin(m2).count('1') >= d:
                dfs(Y + [(y, e)], m2, y + 1)
dfs([], full, 1)
print(f"p={p} d={d}: rank-one bicliques with 0 in X- (X- sorted): {count}", sols[:3])
