import json, re, sys, math
from sympy import factorint, isprime
from pathlib import Path
ROOT = Path(__file__).resolve().parent.parent.parent
T4 = {int(a): int(b) for a, b in (l.split() for l in (ROOT / 'search/data/cati_pasechnik_table4.txt').read_text().splitlines()
                                  if l.strip() and not l.startswith('#'))}
fr = (ROOT / 'FRONTIER.md').read_text(encoding='utf-8')
FR = sorted({int(m.group(2)) for m in re.finditer(r'^\| (\d+) \| (\d+) \|', fr, re.M)})
def pp(n):
    if n < 2: return None
    f = factorint(n)
    return (list(f)[0], list(f.values())[0]) if len(f) == 1 else None
def is_pp(n): return pp(n) is not None
# T-matrix orders (Seberry-Yamada 1992 section 3: every t <= 210 except the list; t-matrices/ adds 73, 97, 103, 107, 109, 113, 131, 151, 163, 181, 197)
EXCL = {73, 79, 83, 89, 97, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 183, 191, 193, 197, 199}
REPO_T = {9, 11, 15, 17, 19, 21, 23, 25, 29, 31, 33, 35, 37, 41, 43, 47, 49, 53, 55, 59, 67, 71, 73, 97, 101, 103, 107, 109, 113, 131, 151, 163, 181, 197}
def T_order(t):
    return t == 1 or (t % 2 == 1 and ((t <= 210 and t not in EXCL) or t in REPO_T))
def golay_odd(m):
    """odd m a product of 3, 5, 11, 13 (lengths of complex Golay pairs; products by Turyn's multiplication)"""
    if m == 1: return True
    f = factorint(m)
    return all(p in (3, 5, 11, 13) for p in f)
def source(q):
    """why a quaternary Hadamard matrix of order q-1 in 2x2 block circulant form is known for the prime power q = 3 (mod 4)"""
    out = []
    if is_pp(2 * q - 3): out.append(f"Paley complex p={2*q-3}")
    if is_pp(q - 2): out.append(f"C+iI, Paley conference p={q-2}")
    m = (q - 1) // 2
    if m > 1 and golay_odd(m): out.append(f"complex Golay pair of length {m}")
    return out
if __name__ == "__main__":
    print("FRONTIER n (4n open):", len(FR))
    print("--- (A) q = 3 mod 4 prime power, 4q open (q in FRONTIER), with the sources of QH(q-1) in block circulant form")
    elig = []
    for q in FR:
        if q % 4 == 3 and is_pp(q):
            s = source(q)
            print(q, "prime power", pp(q), "->", s if s else "no source")
            if s: elig.append(q)
    print("eligible:", elig)
    print("--- (B) Table 4 entries n(m>=4), n = q t, q eligible (4q open), t a T-matrix order >= 3")
    hits = []
    for n, m in sorted(T4.items()):
        if m < 4: continue
        for q in elig:
            if n % q == 0 and n // q >= 3 and T_order(n // q):
                hits.append((n, m, q, n // q)); print(f"n = {n}({m}) = {q} * {n//q}: order 8n = {8*n}")
    print("hits:", len(hits))
    print("--- (C) beyond Table 4: n = q t > 2999, q eligible; smallest orders per q (t T-matrix order up to 210)")
    for q in elig:
        ts = [t for t in range(3, 211, 2) if T_order(t)]
        lst = [(q * t, t) for t in ts if q * t > 2999]
        print(q, "n = q t, first:", lst[:8], "  number with n <= 9999:", sum(1 for n, t in lst if n <= 9999))
