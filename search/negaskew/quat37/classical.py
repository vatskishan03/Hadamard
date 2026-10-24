"""classical-route check for a skew-Hadamard order N = 8v (v odd): prints which standard constructions could give N.  usage: python3 classical.py v ..."""
import sys
from sympy import factorint, isprime, primefactors
def is_pp(n):
    if n<2: return None
    f=factorint(n)
    if len(f)==1:
        p,e=next(iter(f.items())); return (p,e)
    return None
def routes(v):
    N=8*v; out=[]
    # Paley I : N-1 = q prime power (q = 3 mod 4 automatically)
    pp=is_pp(N-1)
    if pp: out.append(f"Paley I: {N}-1 = {pp[0]}^{pp[1]}")
    # 2^u (q+1) with q prime power: u=0..3 (N/2^u integer)
    for u in range(1,4):
        pp=is_pp(N//2**u-1)
        if pp: out.append(f"2^{u}*(q+1): q = {N//2**u-1} = {pp[0]}^{pp[1]} (q mod 8 = {(N//2**u-1)%8}) [Paley skew+doubling if q=3 mod 4; Szekeres/Wallis-Whiteman/Momihara-Xiang type otherwise]")
    # Spence 4(q^2+q+1), doubled: v = q^2+q+1 (prime power q) -> 8v = 2*4(q^2+q+1); also q^2+q+1 = v directly gives 4v
    for q in range(2,int(v**0.5)+2):
        if q*q+q+1==v and is_pp(q): out.append(f"Spence: v = q^2+q+1 with q={q}, 4v then doubled (needs the Spence skew CDS)")
    # skew 4v by Paley / other: 4v-1 prime power -> 4v Paley skew, doubled
    pp=is_pp(4*v-1)
    if pp: out.append(f"Paley I at 4v then doubling: 4v-1 = {pp[0]}^{pp[1]}")
    pp=is_pp(2*v-1)
    if pp: out.append(f"2v-1 = {pp[0]}^{pp[1]} prime power (order 2v is not a Hadamard order: v odd) -- irrelevant unless 4(q+1) type: 4(q+1)=8v needs q=2v-1")
    return out
if __name__=='__main__':
    for a in sys.argv[1:]:
        v=int(a); r=routes(v)
        print(v, 8*v, ('CLASSICAL: '+' | '.join(r)) if r else 'no classical route found')
