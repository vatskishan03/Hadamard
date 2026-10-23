"""Transport a 'goethals-seidel-group' certificate over Z_m x Z_q (gcd(m, q) = 1, mixed radix, index k*q + a) to the
cyclic 'goethals-seidel' form of skew/820: Z_v, v = m q, x <-> (x mod m, x mod q); the multiplier h = 1 (mod m),
h = g^2 (mod q) generates H = {1} x (squares of F_q*); each sequence is listed by the least element of every H-orbit
on which it is -1.   usage: python3 cyclic.py IN.json OUT.json"""
import json, sys
from math import gcd


def primitive_root(q):
    return next(g for g in range(2, q) if len({pow(g, i, q) for i in range(1, q)}) == q - 1)


def to_cyclic(cert):
    m, q = cert['group']
    assert gcd(m, q) == 1
    v = m * q
    h = next(x for x in range(1, v) if x % m == 1 and x % q == pow(primitive_root(q), 2, q))
    H = sorted({pow(h, i, v) for i in range(q)})
    blocks = []
    for s in cert['sequences']:
        cyc = [s[(x % m) * q + (x % q)] for x in range(v)]
        reps = sorted({min(x * u % v for u in H) for x in range(v) if cyc[x] == '-'})
        assert all(cyc[r * u % v] == '-' for r in reps for u in H)
        blocks.append(reps)
    out = dict(order=cert['order'], construction='goethals-seidel', skew=cert['skew'], group=f'Z_{v}', v=v,
               multiplier=h, subgroup_order=len(H), rowsums=cert['rowsums'], blocks=blocks)
    return out


if __name__ == '__main__':
    c = to_cyclic(json.load(open(sys.argv[1])))
    lines = ['{'] + [f' {json.dumps(k)}: {json.dumps(x)},' for k, x in c.items() if k != 'blocks']
    lines += [' "blocks": ['] + [f'  {json.dumps(b)}' + (',' if i < 3 else '') for i, b in enumerate(c['blocks'])] + [' ]', '}']
    open(sys.argv[2], 'w').write('\n'.join(lines) + '\n')
