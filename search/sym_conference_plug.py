#!/usr/bin/env python3
"""Symmetric Hadamard matrices of order 2n(n - 1) from a symmetric conference matrix of order n (verify.py
construction `symmetric-conference-plug`; symmetric/CONFERENCE_PLUG.md).

    python3 search/sym_conference_plug.py targets          census-undecided v = n(n - 1)/2 for Paley's orders n = q + 1
    python3 search/sym_conference_plug.py certify q [...]  writes symmetric/<2q(q + 1)>/certificate.json

q is a prime power = 1 (mod 4); for q = p^k (k > 1) the first primitive monic polynomial in lexicographic order of
its coefficients (c_0, ..., c_{k-1}) defines GF(q).
"""
import hashlib, itertools, json, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "verify"))
import verify  # noqa: E402


def prime_power(q):
    p = next(d for d in range(2, q + 1) if q % d == 0)
    k = 0
    while q % p == 0:
        q //= p
        k += 1
    return (p, k) if q == 1 else None


def field_spec(q):
    p, k = prime_power(q)
    if k == 1:
        return q
    for poly in itertools.product(range(p), repeat=k):
        try:
            verify._paley_conference_field({"p": p, "poly": list(poly)})
            return {"p": p, "poly": list(poly)}
        except ValueError:
            continue
    raise ValueError(f"no primitive polynomial for {q}")


def cert(q):
    c = {"order": 2 * q * (q + 1), "construction": "symmetric-conference-plug", "symmetric": True,
         "q": field_spec(q), "sha256": ""}
    H, _ = verify.BUILDERS[c["construction"]](c)
    sha = hashlib.sha256()
    for chunk in verify.csv_chunks(H):
        sha.update(chunk)
    c["sha256"] = sha.hexdigest()
    verify.build_checked(c)
    return c


def targets():
    census = {int(x) for x in (ROOT / "search/data/symmetric_census_undecided.txt").read_text().split()
              if x.isdigit()} - {2500}
    done = {int(d.name) // 4 for d in (ROOT / "symmetric").iterdir() if d.name.isdigit()}
    for q in range(5, 80, 4):
        pk = prime_power(q)
        if pk and q * (q + 1) // 2 in census:
            v = q * (q + 1) // 2
            print(q, v, 4 * v, "in symmetric/" if v in done else "open")


def main():
    if sys.argv[1] == "targets":
        targets()
    elif sys.argv[1] == "certify":
        for q in map(int, sys.argv[2:]):
            c = cert(q)
            d = ROOT / "symmetric" / str(c["order"])
            d.mkdir(exist_ok=True)
            (d / "certificate.json").write_text(json.dumps(c, indent=1) + "\n")
            print(f"symmetric/{c['order']}: {c['sha256']}")


if __name__ == "__main__":
    main()
