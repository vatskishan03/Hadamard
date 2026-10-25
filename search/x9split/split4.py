"""The order-36w plug-in. Y_a = J - 2 L_a (order 9), L_a the subgroup of Z3^2 orthogonal to the direction a, satisfy
Y_a Y_b^T = J (a != b) and sum_a Y_a Y_a^T = 36 I: a regular 4-set of order 9 (Seberry-Whiteman 1988, Seberry-Zhang
1993). If an Hadamard matrix K of order 4w is a disjoint sum K = K0 + K1 + K2 + K3 of weighing matrices W(4w, w), then
H = sum_a K_a (x) Y_a is an Hadamard matrix of order 36w (H H^T = w I (x) sum Y_a Y_a^T + (K K^T - sum K_a K_a^T) (x) J);
H is symmetric if every K_a is. It is the Z3^2-invariant case of the AG(2,3) hyperplane-block lift of NOTES.md §97."""
import numpy as np
pts = [(x0, x1) for x0 in range(3) for x1 in range(3)]
DIRS = [(0, 1), (1, 0), (1, 1), (1, 2)]
def Y(a):
    return np.array([[1 - 2 * int((a[0] * (x[0] - y[0]) + a[1] * (x[1] - y[1])) % 3 == 0) for y in pts] for x in pts], np.int64)
def plug(parts):
    return sum(np.kron(P, Y(DIRS[a])) for a, P in enumerate(parts))
def is_weighing(P, w):
    return np.array_equal(P @ P.T, w * np.eye(P.shape[0], dtype=np.int64))
