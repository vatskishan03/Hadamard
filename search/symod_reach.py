"""Census values in reach of the symmetric weighted plug-in (search/NOTES.md §32) and the row-sum screen of the Balonin
form.

A symmetric OD(4t; mu, 4t - mu) and a symmetric Hadamard matrix of order h = 4k with 4t = mu h give a symmetric
Hadamard matrix of order mu h (h - 1) = 4 v, v = mu k (4k - 1).  For each symmetric-undecided v of the Doković-Seberry
census (odd v < 2500, k odd, k not itself undecided) this lists the splits (k, mu) and whether the Balonin form passes
the row-sum screen: integer vectors r, rho with r_A^2 + 2 r_B^2 + r_D^2 = mu, rho_A^2 + 2 rho_B^2 + rho_D^2 = 4t - mu,
r_A rho_A + 2 r_B rho_B + r_D rho_D = 0, exactly one of r_A, rho_A odd (A symmetric of odd length t) and
|r_X| + |rho_X| <= t.

Usage: python3 search/symod_reach.py CENSUS.txt [--maxt T]   (pdftotext of the census PDF of 2025-05-28)
"""
import argparse
import math
import re

HDR = r'CENSUS OF SYMMETRIC AND SKEW HADAMARD MATRICES OF ORDER 4v FOR ODD v < 2500\d*'


def grab(txt, start, end):
    s = txt.index(start)
    e = txt.index(end, s)
    body = re.sub(HDR, ' ', txt[s + len(start):e])
    body = re.sub(r'DRAGOMIR.*?N2L', ' ', body, flags=re.S)
    body = re.sub(r'\n\s*\d\s*\n', ' ', body)
    return sorted(set(int(x) for x in re.findall(r'\b\d+\b', body)))


def reps(n):
    out = []
    for a in range(-math.isqrt(n), math.isqrt(n) + 1):
        for b in range(-math.isqrt(n // 2) - 1, math.isqrt(n // 2) + 2):
            r = n - a * a - 2 * b * b
            if r < 0:
                continue
            d = math.isqrt(r)
            if d * d == r:
                for dd in {d, -d}:
                    out.append((a, b, dd))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('census')
    ap.add_argument('--maxt', type=int, default=400)
    a = ap.parse_args()
    txt = open(a.census, encoding='utf-8', errors='replace').read()
    sym = set(grab(txt, 'Hadamard matrices of order 4v is undecided:', 'Remark: In his paper'))
    assert len(sym) == 595, len(sym)
    rows = []
    for k in range(1, 26, 2):
        if k in sym:
            continue
        base = k * (4 * k - 1)
        for mu in range(1, 2500, 2):
            v = mu * base
            if v >= 2500:
                break
            if v in sym:
                rows.append((k * mu, v, k, mu))
    rows.sort()
    print(f'{len({r[1] for r in rows})} symmetric-undecided v in reach, {len(rows)} splits')
    for t, v, k, mu in rows:
        if t > a.maxt:
            continue
        ok = set()
        for r in reps(mu):
            for p in reps(4 * t - mu):
                if (r[0] * p[0] + 2 * r[1] * p[1] + r[2] * p[2] == 0 and r[0] % 2 != p[0] % 2
                        and all(abs(r[i]) + abs(p[i]) <= t for i in range(3))):
                    ok.add((tuple(map(abs, r)), tuple(map(abs, p))))
        print(f't={t:4d} v={v:4d} h={4 * k:3d} mu={mu:4d}  ' +
              (f'{len(ok)} row-sum classes, e.g. {sorted(ok)[:2]}' if ok else 'fails the row-sum screen'))


if __name__ == '__main__':
    main()
