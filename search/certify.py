#!/usr/bin/env python3
"""Turn a match into a repository certificate.

Expands the four masks of a "SOL" line printed by match (coset numbering as in
generate.py), checks that the sequences are complementary, and writes
certificate.json, sequences.txt and H<n>.csv.gz into the order directory. The
result is then checked with verify/verify.py. An existing certificate.json is
never overwritten.

    python3 search/certify.py 1123 60 "SOL 2c76b2b30 21a71af94 a95af712 4dcf4c15" /path/to/H4492

The directory is an absolute path, or a path relative to the repository root; a bare number such as 4492 means
hadamard/4492.
"""
import argparse
import gzip
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "verify"))
import verify  # noqa: E402


def coset_numbering(p, g):
    orb = -np.ones(p, dtype=np.int64)
    k = 0
    for j in range(p):
        if orb[j] < 0:
            x = j
            while orb[x] < 0:
                orb[x] = k
                x = x * g % p
            k += 1
    return orb


def write_certificate(cert, path):
    if path.exists():
        raise SystemExit(f"{path} exists; not overwritten")
    lines = ["{"]
    keys = list(cert)
    for i, k in enumerate(keys):
        if k == "blocks":
            body = ",\n".join("  " + json.dumps(b) for b in cert[k])
            text = f'"blocks": [\n{body}\n ]'
        else:
            text = f'"{k}": {json.dumps(cert[k])}'
        lines.append(" " + text + ("," if i < len(keys) - 1 else ""))
    path.write_text("\n".join(lines + ["}"]) + "\n")


def main():
    ap = argparse.ArgumentParser(description="Write a certificate for a Goethals-Seidel solution.")
    ap.add_argument("p", type=int)
    ap.add_argument("g", type=int, help="generator of the multiplier subgroup")
    ap.add_argument("sol", help='"SOL maskA maskB maskC maskD" as printed by match')
    ap.add_argument("directory", help="order directory to write into (a bare number such as 4492 means "
                                      "hadamard/4492; other paths are absolute or relative to the repository root)")
    a = ap.parse_args()

    p, g = a.p, a.g
    masks = [int(m, 16) for m in a.sol.split()[1:5]]
    orb = coset_numbering(p, g)
    seqs = [np.array([-1 if o and (m >> (o - 1)) & 1 else 1 for o in orb], dtype=np.int64) for m in masks]
    paf = verify.paf_sum(seqs)
    if paf[0] != 4 * p or np.any(paf[1:] != 0):
        sys.exit("not a complementary quadruple")

    H = sorted({pow(g, i, p) for i in range(p)})
    blocks = []
    for x in seqs:
        neg = set(np.flatnonzero(x < 0).tolist())
        reps, seen = [], set()
        for t in sorted(neg):
            if t not in seen:
                seen |= {t * h % p for h in H}
                reps.append(t)
        blocks.append(reps)

    order = 4 * p
    out = ROOT / "hadamard" / a.directory if a.directory.isdigit() else ROOT / a.directory
    out.mkdir(parents=True, exist_ok=True)
    M = verify.goethals_seidel(*map(verify.circulant, seqs))
    data = verify.canonical_csv(M)
    cert = {"order": order, "construction": "goethals-seidel", "group": f"Z_{p}", "v": p,
            "multiplier": g, "subgroup_order": len(H), "rowsums": [int(x.sum()) for x in seqs],
            "blocks": blocks, "sha256": hashlib.sha256(data).hexdigest()}
    write_certificate(cert, out / "certificate.json")
    (out / "sequences.txt").write_text("".join("".join("+" if v > 0 else "-" for v in x) + "\n" for x in seqs))
    with gzip.open(out / f"H{order}.csv.gz", "wb", compresslevel=9) as f:
        f.write(data)
    n, digest = verify.check(out)
    print(f"order {n}: certificate written to {out}, sha256 {digest}")


if __name__ == "__main__":
    main()
