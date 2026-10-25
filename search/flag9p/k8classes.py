"""Switching classes (with D-labels: vertices 0-3 have D = +1, 4-7 have D = -1) of the K8 templates printed by `kn 8 4 N`
(lines `SOL k_01 k_02 ... k_67`, sixth roots e^{i pi k/3}, sign + for even k). Canonical form: the triple products
s_ab s_bc s_ac under all relabelings inside the two D-classes. Writes k8_classes.pkl.

    ./kn 8 4 100000 > k8_all.txt && python3 k8classes.py k8_all.txt"""
import itertools, pickle, sys
lines = [l for l in open(sys.argv[1]) if l.startswith('SOL')]
n = 8
pairs = [(a, b) for a in range(n) for b in range(a + 1, n)]
def canon(sgn):
    tri = lambda a, b, c: sgn[tuple(sorted((a, b)))] * sgn[tuple(sorted((b, c)))] * sgn[tuple(sorted((a, c)))]
    T = {t: tri(*t) for t in itertools.combinations(range(n), 3)}
    best = None
    for p1 in itertools.permutations(range(4)):
        for p2 in itertools.permutations(range(4, 8)):
            perm = list(p1) + list(p2)
            key = tuple(T[tuple(sorted((perm[a], perm[b], perm[c])))] for a, b, c in itertools.combinations(range(n), 3))
            if best is None or key < best: best = key
    return best
classes = set()
for l in lines:
    ks = list(map(int, l.split()[1:]))
    classes.add(canon({pr: (1 if k % 2 == 0 else -1) for pr, k in zip(pairs, ks)}))
print("solutions", len(lines), "classes with D-labels", len(classes))
pickle.dump(classes, open('k8_classes.pkl', 'wb'))
