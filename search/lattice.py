"""Integer lattice tools: HNF of the lattice spanned by integer vectors, membership tests."""

def hnf_basis(vectors):
    """Row-style Hermite basis (Python ints) of the Z-span of the given integer row vectors."""
    rows = [list(map(int, v)) for v in vectors if any(int(x) for x in v)]
    if not rows:
        return []
    n = len(rows[0])
    basis = []                       # list of (pivot col, row)
    for v in rows:
        v = v[:]
        for i, (pc, b) in enumerate(basis):
            pass
        basis = _insert(basis, v, n)
    return [b for _, b in basis]

def _insert(basis, v, n):
    """Insert v into an echelon basis (sorted by pivot column) keeping it echelon (gcd steps)."""
    basis = [list(x) for x in basis]
    out = [(pc, b[:]) for pc, b in basis]
    while True:
        nz = [j for j in range(n) if v[j] != 0]
        if not nz:
            return out
        j = nz[0]
        k = next((i for i, (pc, _) in enumerate(out) if pc == j), None)
        if k is None:
            out.append((j, v))
            out.sort(key=lambda t: t[0])
            return out
        b = out[k][1]
        # extended gcd on (b[j], v[j])
        a, c = b[j], v[j]
        g, x, y = _egcd(a, c)
        newb = [x * bi + y * vi for bi, vi in zip(b, v)]
        newv = [(c // g) * bi - (a // g) * vi for bi, vi in zip(b, v)]
        if newb[j] < 0:
            newb = [-t for t in newb]
        out[k] = (j, newb)
        v = newv

def _egcd(a, b):
    if b == 0:
        return (abs(a), (1 if a >= 0 else -1), 0)
    x0, y0, x1, y1 = 1, 0, 0, 1
    while b:
        q, a, b = a // b, b, a % b
        x0, x1 = x1, x0 - q * x1
        y0, y1 = y1, y0 - q * y1
    if a < 0:
        a, x0, y0 = -a, -x0, -y0
    return a, x0, y0

def in_lattice(basis, w):
    w = list(map(int, w))
    for b in sorted(basis, key=lambda b: next(j for j, x in enumerate(b) if x)):
        j = next(i for i, x in enumerate(b) if x)
        if w[j] % b[j]:
            return False
        q = w[j] // b[j]
        w = [wi - q * bi for wi, bi in zip(w, b)]
    return not any(w)

def index_info(basis, dim):
    """(rank, product of pivots) - product of pivots = index in Z^dim when full rank."""
    piv = 1
    for b in basis:
        j = next(i for i, x in enumerate(b) if x)
        piv *= abs(b[j])
    return len(basis), piv
