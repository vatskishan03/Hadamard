#!/usr/bin/env python3
"""Symmetric Hadamard matrices of order 4q, q = 1 (mod 4) a prime, from a periodic Golay pair of length (q - 1)/2
(verify.py construction `symmetric-involution`; symmetric/PERIODIC_GOLAY.md, search/NOTES.md section 65).

    python3 search/sym_involution.py targets          census-undecided primes q = 5 (mod 8) and the input K each needs
    python3 search/sym_involution.py certify q [...]  writes symmetric/<4q>/certificate.json

A directory that already holds a certificate is never overwritten.

The pairs: length 50 over Z_2 x Z_5^2, from an exhaustive enumeration over GF(25) of ternary quadruples constant
on the cosets of the subgroup of index 8 (program not included); length 74 over Z_74, Dokovic and Kotsireas (Numer.
Algorithms 69 (2015) 523-530, arXiv:1310.5773), parameter set (74; 36, 31; 30), the sets as unions of orbits of
H = {1, 47, 63}.
"""
import hashlib, json, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "verify"))
sys.path.insert(0, str(ROOT / "search"))
import verify  # noqa: E402
import sym_propus_miyamoto as SPM  # noqa: E402


def _signs(row):
    return "".join("+" if x > 0 else "-" for x in row)


def _orbit_seq(v, H, reps):
    a = [1] * v
    for j in reps:
        for h in H:
            a[j * h % v] = -1
    return a


PAIRS = {
    101: dict(group="2x5^2", t=25,
              a="++-+-+-+++--++++++-+--+-++++-+++---++------+-++-++",
              b="+--++-++++-+++++---+++-+-----+-++++-+-++----+++-+-",
              source="enumeration over GF(25), index 8 (search/NOTES.md section 65; program not included)"),
    149: dict(group="74", t=37,
              a=_signs(_orbit_seq(74, [1, 47, 63], [1, 4, 6, 7, 9, 12, 22, 23, 28, 29, 34, 42])),
              b=_signs(_orbit_seq(74, [1, 47, 63], [1, 2, 4, 6, 9, 12, 17, 21, 22, 37, 55])),
              source="Dokovic-Kotsireas, arXiv:1310.5773, (74; 36, 31; 30)"),
}


def cert(q):
    p = PAIRS[q]
    c = {"order": 4 * q, "construction": "symmetric-involution", "symmetric": True, "q": q,
         "group": p["group"], "t": p["t"],
         "a": [1 if ch == "+" else -1 for ch in p["a"]], "b": [1 if ch == "+" else -1 for ch in p["b"]],
         "sha256": ""}
    H, _ = verify.BUILDERS[c["construction"]](c)
    sha = hashlib.sha256()
    for chunk in verify.csv_chunks(H):
        sha.update(chunk)
    c["sha256"] = sha.hexdigest()
    verify.build_checked(c)
    return c


def targets():
    done = {int(d.name) // 4 for d in (ROOT / "symmetric").iterdir() if d.name.isdigit()}

    def two_squares(n):
        return any(int((n - a * a) ** 0.5) ** 2 == n - a * a for a in range(int(n ** 0.5) + 1))

    for q in sorted(set(SPM.undecided()) - done):
        if q % 8 == 5 and SPM.is_prime(q):
            w = (q - 1) // 4
            need = (f"K = [[P, Q], [Q, P]] of order {q - 1}, e.g. a periodic Golay pair of length {(q - 1) // 2}"
                    if two_squares(w) else "no input K exists: (q - 1)/4 is not a sum of two squares")
            print(q, 4 * q, need, "(available)" if q in PAIRS else "")


def certify(qs):
    for q in qs:
        if q not in PAIRS:
            raise SystemExit(f"no periodic Golay pair for q = {q} (available: {', '.join(map(str, PAIRS))})")
    dirs = {q: ROOT / "symmetric" / str(4 * q) for q in qs}
    for d in dirs.values():
        if (d / "certificate.json").exists():
            raise SystemExit(f"{d} exists")
    for q, d in dirs.items():
        c = cert(q)
        d.mkdir(exist_ok=True)
        (d / "certificate.json").write_text(json.dumps(c, indent=1) + "\n")
        print(f"symmetric/{4 * q}: {c['sha256']}")


if __name__ == "__main__":
    if sys.argv[1:2] == ["targets"]:
        targets()
    elif sys.argv[1:2] == ["certify"]:
        certify([int(x) for x in sys.argv[2:]])
    else:
        raise SystemExit(__doc__)
