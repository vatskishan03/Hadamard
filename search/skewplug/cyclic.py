"""Transport a 'goethals-seidel-group' certificate over Z_m x Z_q (gcd(m, q) = 1, mixed radix, index k*q + a) to the
cyclic 'goethals-seidel' form of skew/820: Z_v, v = m q, x <-> (x mod m, x mod q). Each diagonal automorphism
diag(u, w) of the certificate becomes the multiplier h = u (mod m), h = w (mod q) (without automorphisms: h = 1 mod m,
h = g^2 mod q). Each sequence is listed by the least element of every orbit of the multiplier group on which it is -1.
usage: python3 cyclic.py IN.json OUT.json"""
import json, sys
from math import gcd


def primitive_root(q):
    return next(g for g in range(2, q) if len({pow(g, i, q) for i in range(1, q)}) == q - 1)


def crt(a, m, b, q):
    return next(x for x in range(m * q) if x % m == a % m and x % q == b % q)


def to_cyclic(cert):
    m, q = cert['group']
    assert gcd(m, q) == 1
    v = m * q
    autos = cert.get('automorphisms') or [[[1, 0], [0, pow(primitive_root(q), 2, q)]]]
    for M in autos:
        assert M[0][1] == 0 and M[1][0] == 0, 'only diagonal automorphisms'
    gens = [crt(M[0][0], m, M[1][1], q) for M in autos]
    H = {1}
    while True:
        grown = H | {h * g % v for h in H for g in gens}
        if grown == H:
            break
        H = grown
    blocks = []
    for s in cert['sequences']:
        cyc = [s[(x % m) * q + (x % q)] for x in range(v)]
        reps = sorted({min(x * u % v for u in H) for x in range(v) if cyc[x] == '-'})
        assert all(cyc[r * u % v] == '-' for r in reps for u in H)
        blocks.append(reps)
    out = dict(order=cert['order'], construction='goethals-seidel', skew=cert['skew'], group=f'Z_{v}', v=v)
    if len(gens) == 1:
        out['multiplier'] = gens[0]
    else:
        out['multipliers'] = gens
    out.update(subgroup_order=len(H), rowsums=cert['rowsums'], blocks=blocks)
    return out


def dump(c):
    lines = ['{'] + [f' {json.dumps(k)}: {json.dumps(x)},' for k, x in c.items() if k != 'blocks']
    lines += [' "blocks": ['] + [f'  {json.dumps(b)}' + (',' if i < 3 else '') for i, b in enumerate(c['blocks'])] + [' ]', '}']
    return '\n'.join(lines) + '\n'


if __name__ == '__main__':
    open(sys.argv[2], 'w').write(dump(to_cyclic(json.load(open(sys.argv[1])))))
