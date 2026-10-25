"""Đoković's two Goethals-Seidel quadruples of order 191 (arXiv:math/0703312): H = <39> of order 5, cosets
alpha_{2i} = r_i H, alpha_{2i+1} = -r_i H, S_k = union of alpha_i over J_k, a_k = -1 on S_k."""
import numpy as np
n = 191
H = [1, 39, 184, 109, 49]
reps = [1, 2, 3, 4, 6, 8, 9, 11, 12, 13, 16, 17, 18, 19, 22, 32, 36, 38, 41]
alpha = []
for r in reps:
    alpha.append(sorted(r * h % n for h in H)); alpha.append(sorted((-r * h) % n for h in H))
J = [[[1, 7, 9, 10, 11, 13, 17, 18, 25, 26, 30, 31, 33, 34, 35, 36, 37],
      [1, 4, 7, 9, 11, 12, 13, 14, 19, 21, 22, 23, 24, 25, 26, 29, 36, 37],
      [0, 3, 4, 5, 7, 8, 9, 16, 17, 19, 24, 25, 29, 30, 31, 33, 35, 37],
      [1, 3, 4, 5, 8, 11, 14, 18, 19, 20, 21, 23, 24, 25, 28, 29, 30, 32, 34, 35]],
     [[0, 1, 6, 8, 9, 11, 12, 16, 18, 20, 21, 23, 28, 31, 33, 36, 37],
      [0, 1, 3, 4, 10, 12, 13, 17, 20, 22, 24, 31, 32, 33, 34, 35, 36, 37],
      [4, 8, 9, 10, 12, 13, 14, 16, 17, 20, 21, 24, 26, 27, 29, 31, 32, 34],
      [1, 7, 9, 10, 11, 12, 14, 15, 16, 17, 20, 22, 23, 25, 28, 29, 32, 33, 34, 37]]]
def quad(which):
    seqs = []
    for Jk in J[which]:
        x = np.ones(n, np.int64)
        for i in Jk: x[alpha[i]] = -1
        seqs.append(x)
    return seqs
def slots():
    """slot index of each residue: 0 for {0}, 1 + i for coset alpha_i"""
    s = np.zeros(n, np.int64)
    for i, a in enumerate(alpha): s[a] = 1 + i
    return s
if __name__ == '__main__':
    assert sorted(x for a in alpha for x in a) == list(range(1, n))
    for w in (0, 1):
        Q = quad(w)
        F = sum(np.abs(np.fft.fft(x)) ** 2 for x in Q)
        paf = np.rint(np.fft.ifft(F).real).astype(int)
        print(w, [int((x == -1).sum()) for x in Q], paf[0], int(np.abs(paf[1:]).max()))
