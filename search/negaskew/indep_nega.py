"""independent check of a 'goethals-seidel-negacyclic' certificate from the file only: four +-1 sequences of length n,
negacyclic blocks, plain reversal, Goethals-Seidel array, no complex arithmetic.   usage: python3 indep_nega.py CERT.json"""
import json, sys, numpy as np
c = json.load(open(sys.argv[1]))
S = [[1 if ch == '+' else -1 for ch in s] for s in c['sequences']]; n = len(S[0])
a = S[0]; print('order', 4 * n, 'n', n, 'a0', a[0], 'A palindromic', all(a[j] == a[n - j] for j in range(1, n)))
def nega_paf(s, k): return sum((s[i] * s[i + k] if i + k < n else -s[i] * s[i + k - n]) for i in range(n))
tot = [sum(nega_paf(s, k) for s in S) for k in range(n)]
print('sum of negaperiodic PAFs: peak', tot[0], 'max |off-peak|', max(abs(x) for x in tot[1:]))
idx = np.arange(n)
def N(s):
    s = np.array(s, dtype=np.float64); d = idx[None, :] - idx[:, None]
    return np.where(d >= 0, s[d % n], -s[d % n])
A, B, C, D = [N(s) for s in S]; R = np.eye(n)[::-1]; T = lambda X: X.T
H = np.block([[A, B @ R, C @ R, D @ R], [-B @ R, A, T(D) @ R, -T(C) @ R], [-C @ R, -T(D) @ R, A, T(B) @ R], [-D @ R, T(C) @ R, -T(B) @ R, A]]); m = H.shape[0]
print('+-1:', bool(np.all(np.abs(H) == 1)), ' H H^T = mI:', bool(np.array_equal(H @ H.T, m * np.eye(m))), ' H + H^T = 2I:', bool(np.array_equal(H + H.T, 2 * np.eye(m))))
