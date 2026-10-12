"""Write a repository certificate for a skew Goethals-Seidel hit recorded in skew_hits.json.

    python3 attack/skew_certify.py n [index]      # index-th verified hit for n (default 0)

Writes skew/<4n>/certificate.json, sequences.txt and H<4n>.csv.gz and checks them with verify/verify.py
(H H^T = 4n I, H + H^T = 2I, SHA-256)."""
import gzip, hashlib, json, sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "verify")); sys.path.insert(0, str(ROOT / "search"))
import verify
from certify import write_certificate


def subgroup(gens, v):
    H = {1}
    while True:
        grown = H | {h * g % v for h in H for g in gens}
        if grown == H: return sorted(H)
        H = grown


def main():
    n = int(sys.argv[1]); k = int(sys.argv[2]) if len(sys.argv) > 2 else 0
    hits = [json.loads(l) for l in open(ROOT / "attack" / "skew_hits.json")]
    hits = [h for h in hits if h["n"] == n and h["hadamard"] and h["skew"]]
    h = hits[k]
    seqs = [np.array([1 if c == '+' else -1 for c in s], dtype=np.int64) for s in h["seqs"]]
    A = seqs[0]
    assert A[0] == 1 and all(A[(-t) % n] == -A[t] for t in range(1, n)), "first block is not skew-type"
    paf = verify.paf_sum(seqs)
    assert paf[0] == 4 * n and not paf[1:].any()
    gens = h["gens"]; H = subgroup(gens, n)
    blocks = []
    for x in seqs:
        assert all(x[t * g % n] == x[t] for t in range(n) for g in gens), "block not invariant under H"
        reps, seen = [], set()
        for t in sorted(np.flatnonzero(x < 0).tolist()):
            if t not in seen:
                seen |= {t * u % n for u in H}; reps.append(t)
        blocks.append(reps)
    order = 4 * n; out = ROOT / "skew" / str(order); out.mkdir(parents=True, exist_ok=True)
    M = verify.goethals_seidel(*map(verify.circulant, seqs))
    data = verify.canonical_csv(M)
    cert = {"order": order, "construction": "goethals-seidel", "skew": True, "group": f"Z_{n}", "v": n}
    if len(gens) == 1: cert["multiplier"] = gens[0]
    else: cert["multipliers"] = gens
    cert.update({"subgroup_order": len(H), "rowsums": [int(x.sum()) for x in seqs], "blocks": blocks,
                 "sha256": hashlib.sha256(data).hexdigest()})
    write_certificate(cert, out / "certificate.json")
    (out / "sequences.txt").write_text("".join("".join("+" if v > 0 else "-" for v in x) + "\n" for x in seqs))
    with gzip.open(out / f"H{order}.csv.gz", "wb", compresslevel=9) as f: f.write(data)
    m, digest = verify.check(out)
    print(f"skew order {m}: H = <{','.join(map(str, gens))}> of order {len(H)}, rows {cert['rowsums']}, sha256 {digest}")


if __name__ == "__main__":
    main()
