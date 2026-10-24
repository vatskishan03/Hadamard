"""recursive8p.py -- exact checks of the recursive / direct theorems of Tables 3 and 4 of the survey (Egan-Kotsireas-London-O Cathain) at the six orders 8p.
Every statement is a decision from the stated hypotheses; hypotheses on ingredients are taken from Table 3 of the survey."""
from sympy import factorint, isprime
import itertools
P = [599, 719, 739, 839, 1031, 1087]

def is_pp(n):
    f = factorint(n); return len(f) == 1
def pp_info(n):
    f = factorint(n); return dict(f)

# Golay lengths 2^a 10^b 26^c ; complex Golay: 1,3,5,11,13, even <=26 except 14, 2mn (survey), optionally mn (Turyn-product assumption, generous)
LIM = 2200
golay = set()
for a in range(0, 12):
    for b in range(0, 4):
        for c in range(0, 3):
            v = 2**a * 10**b * 26**c
            if v <= LIM: golay.add(v)
cg = {1, 3, 5, 11, 13} | {n for n in range(2, 27, 2) if n != 14}
def close(base, mult2=True, mult1=False):
    S = set(base); ch = True
    while ch:
        ch = False
        for m in list(S):
            for n in list(S):
                for v in ((2*m*n,) if mult2 else ()) + ((m*n,) if mult1 else ()):
                    if v <= LIM and v not in S: S.add(v); ch = True
    return S
cg_strict = close(cg, True, False)
cg_generous = close(cg, True, True)

print("# Table 4 / Table 3 constructions, order 8p = 4n, n = 2p")
for p in P:
    n = 2 * p; N = 8 * p
    res = []
    # Paley I: q+1 = 8p
    q = N - 1; res.append(("Paley I (q = 8p-1 prime power)", is_pp(q), pp_info(q)))
    # Paley II: 2(q+1) = 8p, q = 4p-1 must be = 1 mod 4
    q = 4 * p - 1; res.append(("Paley II (q = 4p-1 = 1 mod 4 prime power)", is_pp(q) and q % 4 == 1, (q % 4, pp_info(q))))
    # Paley blocks: 4 s t, 4s-1 = 3 mod 8 prime power, T-matrices of order t>1: s*t = 2p -> (s,t) = (p,2), (2,p)
    q = 4 * p - 1; res.append(("Paley blocks (s,t)=(p,2): 4p-1 = 3 mod 8 prime power, T(2) exists", is_pp(q) and q % 8 == 3, (q % 8, pp_info(q))))
    res.append(("Paley blocks (s,t)=(2,p): q = 7 (=7 mod 8) needs T-matrices of order p: known t <= 197 and 97, 103, 109, 151, 163, 181 only", p <= 197, p))
    # Miyamoto III: q=p, 2q-3 prime
    res.append(("Miyamoto III (8q, q = 3 mod 4, q and 2q-3 prime)", isprime(p) and isprime(2*p-3), (2*p-3, pp_info(2*p-3))))
    # Seberry-Yamada 8q given WT((q-1)/2)
    w = (p - 1) // 2
    turyn = is_pp(2*w - 1)
    res.append((f"Seberry-Yamada 8q needs WT({w}): Turyn (2w-1 = {2*w-1} prime power)?", turyn, pp_info(2*w-1)))
    res.append((f"  WT({w}) via Miyamoto II (w = 1 mod 4 prime power + symmetric BH)", w % 4 == 1 and is_pp(w), (w % 4, pp_info(w))))
    res.append((f"  WT({w}) in Table A.1 of 1992 / odd w <= 69 or listed", w <= 69 or w in (73,79,91,93,95,97,103,109,127,133,181,111), w))
    # Cooper-Wallis 4tw = 8p: t*w = 2p
    res.append(("Cooper-Wallis/Williamson: needs WT(p) (t=2), T(p) (w=2), T(2p) (w=1) or WT(2p)", False, "no WT/T of these orders known (T-matrices t<=197)"))
    # complex Golay 8(x+y), x+y=p
    ps = [(x, p - x) for x in sorted(cg_generous) if x <= p - x and (p - x) in cg_generous]
    ps2 = [(x, p - x) for x in sorted(cg_strict) if x <= p - x and (p - x) in cg_strict]
    res.append(("Craigen-Holzmann-Kharaghani 8(x+y), x+y = p, complex Golay lengths (survey: 2mn closure)", bool(ps2), ps2))
    res.append(("  same with the generous closure under mn as well", bool(ps), ps))
    # T-sequences from real Golay: length g+g' = 2p
    gp = [(x, 2*p - x) for x in sorted(golay) if x <= 2*p - x and (2*p - x) in golay]
    res.append(("T-sequences of length g+g' = 2p (Golay lengths 2^a 10^b 26^c) -> H(4*2p)", bool(gp), gp))
    # Spence 4(1+q+q^2) / Mathon / Yamada / Sylvester / Agaian
    sp = [qq for qq in range(2, 200) if 1 + qq + qq*qq == 2*p]
    res.append(("Spence 4(1+q+q^2) = 8p", bool(sp), sp))
    mt = [qq for qq in range(3, 200, 4) if is_pp(qq) and is_pp(qq+2) and 2*(qq*qq*(qq+2)+1) == N]
    res.append(("Mathon 2(q^2(q+2)+1) = 8p", bool(mt), mt))
    y = [qq for qq in [2*p-2] if is_pp(qq)]
    res.append(("Yamada 4(q+2) = 8p (q = 2p-2 prime power)", bool(y), pp_info(2*p-2)))
    res.append(("Agaian 8ab, ab = p", False, "a = 1, b = p needs H(4p): open"))
    res.append(("Miyamoto I 4q = 8p (q = 2p not = 1 mod 4 prime power); Miyamoto II 4tq: q|2p, q = 1 mod 4 prime power", False, "none"))
    res.append(("Craigen-Seberry-Zhang 16abcd: 16 does not divide 8p", False, ""))
    print(f"\n p = {p}  (8p = {N})")
    for name, ok, info in res:
        print(f"   {'AVAILABLE' if ok else 'fails    '}  {name}   [{info}]")
