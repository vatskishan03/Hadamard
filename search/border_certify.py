"""Write a repository certificate for a bordered propus hit (propus.py / propus_chain.py / propus_bgroup.py with
--border) recorded in $HADAMARD_WORK/sym_hits.jsonl.

    python3 search/border_certify.py ORDER [index]      # index-th bordered hit of that order (default 0)

Writes symmetric/<ORDER>/certificate.json (construction "propus-bordered"), sequences.txt and H<ORDER>.csv.gz and
checks them with verify/verify.py (H = H^T, H H^T = ORDER I, SHA-256)."""
import gzip, hashlib, json, sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parent.parent
sys.path[:0] = [str(ROOT / "attack"), str(ROOT / "verify"), str(ROOT / "search")]
import verify
import workdir as W
from certify import write_certificate


def qborder_R(h, group, cs, q):
    """Hadamard R of order 4q solving the border problem for the compressions of the hit, with the propus array over
    Z_q built as verify.build_propus_qbordered builds it (reflection R[c, -c] = 1)."""
    import propus as PR
    mods = np.array(group, dtype=np.int64)
    coords = np.array(np.unravel_index(np.arange(int(np.prod(mods))), tuple(int(x) for x in mods))).T
    phi = (coords @ np.array(cs, dtype=np.int64)) % q
    xs = [np.array([1 if c == '+' else -1 for c in s], dtype=np.int64) for s in h["seqs"]]
    ys = [np.array([int(x[phi == c].sum()) for c in range(q)], dtype=np.int64) for x in xs]
    zq = np.arange(q); diff = (zq[None, :] - zq[:, None]) % q
    Rq = np.zeros((q, q), dtype=np.int64); Rq[zq, (-zq) % q] = 1
    cA, cB, cD = (y[diff] for y in ys)
    Mc = np.block([[cA, cB, cB, cD], [cB.T, cD @ Rq, -cA, -cB @ Rq], [cB.T, -cA, -cD @ Rq, cB @ Rq], [cD.T, -cB @ Rq, cB @ Rq, -cA]])
    Rh = PR.qborder_solve(Mc, q)
    if Rh is None: raise SystemExit("no border for these compressions")
    return Rh.tolist()


def main():
    order = int(sys.argv[1]); k = int(sys.argv[2]) if len(sys.argv) > 2 else 0
    hits = [json.loads(l) for l in open(W.hits('sym_hits.jsonl'))]
    hits = [h for h in hits if h.get("border") and h.get("order") == order]
    h = hits[k]
    group = h.get("group") or [h["n"]]
    q = h.get("q", 1)
    if q == 1:
        cert = {"order": order, "construction": "propus-bordered", "symmetric": True, "group": group,
                "rowsums": h["rows"], "K": h["K"], "sequences": h["seqs"]}
    else:                                               # Z_q border: solve the border problem in verify's convention
        cs = h.get("quotient") or [1]
        cert = {"order": order, "construction": "propus-qbordered", "symmetric": True, "group": group, "q": q,
                "quotient": cs, "rowsums": h["rows"], "R": qborder_R(h, group, cs, q), "sequences": h["seqs"]}
    if h.get("automorphisms"):
        cert["automorphisms"] = h["automorphisms"]
    elif h.get("gens"):
        cert["multipliers"] = h["gens"]                 # Z_m: x -> g x for each multiplier g (informative)
    M, _ = (verify.build_propus_bordered if q == 1 else verify.build_propus_qbordered)(cert)
    data = verify.canonical_csv(M)
    cert["sha256"] = hashlib.sha256(data).hexdigest()
    out = ROOT / "symmetric" / str(order); out.mkdir(parents=True, exist_ok=True)
    write_certificate(cert, out / "certificate.json")
    (out / "sequences.txt").write_text("".join(s + "\n" for s in h["seqs"]))
    with gzip.open(out / f"H{order}.csv.gz", "wb", compresslevel=9) as f: f.write(data)
    m, digest = verify.check(out)
    print(f"symmetric order {m} (Z_{q}-bordered propus over Z_{'x Z_'.join(map(str, group))}), rows {h['rows']}, sha256 {digest}")


if __name__ == "__main__":
    main()
