"""Write a repository certificate for a symmetric propus hit recorded in $HADAMARD_WORK/sym_hits.jsonl.

    python3 search/sym_certify.py n [index]      # index-th hit for n (default 0)

The hit must come from the array with only A symmetric (`propus.py --freeD`, verify.py construction "propus").
Writes symmetric/<4n>/certificate.json, sequences.txt and H<4n>.csv.gz and checks them with verify/verify.py
(H H^T = 4n I, H = H^T, SHA-256)."""
import gzip, hashlib, json, sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parent.parent
sys.path[:0] = [str(ROOT / "verify"), str(ROOT / "search")]
import verify
import workdir as W
from certify import write_certificate


def subgroup(gens, v):
    H = {1}
    while True:
        grown = H | {h * g % v for h in H for g in gens}
        if grown == H: return sorted(H)
        H = grown


def main():
    n = int(sys.argv[1]); k = int(sys.argv[2]) if len(sys.argv) > 2 else 0
    hits = [json.loads(l) for l in open(W.hits('sym_hits.jsonl'))]
    hits = [h for h in hits if h["n"] == n and h.get("freeD")]
    h = hits[k]
    seqs = [np.array([1 if c == '+' else -1 for c in s], dtype=np.int64) for s in h["seqs"]]
    A, B, D = seqs
    assert np.array_equal(A, A[(-np.arange(n)) % n]), "first sequence is not symmetric"
    paf = verify.paf_sum([A, B, D]) + verify.paf_sum([B])
    assert paf[0] == 4 * n and not paf[1:].any()
    gens = h["gens"]; H = subgroup(gens, n)
    blocks = []
    for x in seqs:
        assert all(x[t * g % n] == x[t] for t in range(n) for g in gens), "sequence not invariant under H"
        reps, seen = [], set()
        for t in sorted(np.flatnonzero(x < 0).tolist()):
            if t not in seen:
                seen |= {t * u % n for u in H}; reps.append(t)
        blocks.append(reps)
    order = 4 * n; out = ROOT / "symmetric" / str(order); out.mkdir(parents=True, exist_ok=True)
    cert = {"order": order, "construction": "propus", "symmetric": True, "group": f"Z_{n}", "v": n,
            "multipliers": gens, "subgroup_order": len(H), "rowsums": [int(x.sum()) for x in seqs], "blocks": blocks}
    M, _ = verify.build_propus(cert)
    data = verify.canonical_csv(M)
    cert["sha256"] = hashlib.sha256(data).hexdigest()
    write_certificate(cert, out / "certificate.json")
    (out / "sequences.txt").write_text("".join("".join("+" if v > 0 else "-" for v in x) + "\n" for x in seqs))
    with gzip.open(out / f"H{order}.csv.gz", "wb", compresslevel=9) as f: f.write(data)
    m, digest = verify.check(out)
    print(f"symmetric order {m}: H = <{','.join(map(str, gens))}> of order {len(H)}, rows {cert['rowsums']}, sha256 {digest}")


if __name__ == "__main__":
    main()
