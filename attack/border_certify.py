"""Write a repository certificate for a bordered propus hit (propus.py / propus_chain.py / propus_bgroup.py with
--border) recorded in $HADAMARD_WORK/sym_hits.jsonl.

    python3 attack/border_certify.py ORDER [index]      # index-th bordered hit of that order (default 0)

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


def main():
    order = int(sys.argv[1]); k = int(sys.argv[2]) if len(sys.argv) > 2 else 0
    hits = [json.loads(l) for l in open(W.hits('sym_hits.jsonl'))]
    hits = [h for h in hits if h.get("border") and h.get("order") == order]
    h = hits[k]
    group = h.get("group") or [h["n"]]
    cert = {"order": order, "construction": "propus-bordered", "symmetric": True, "group": group,
            "rowsums": h["rows"], "K": h["K"], "sequences": h["seqs"]}
    if h.get("automorphisms"):
        cert["automorphisms"] = h["automorphisms"]
    elif h.get("gens"):
        cert["multipliers"] = h["gens"]                 # Z_m: x -> g x for each multiplier g (informative)
    M, _ = verify.build_propus_bordered(cert)
    data = verify.canonical_csv(M)
    cert["sha256"] = hashlib.sha256(data).hexdigest()
    out = ROOT / "symmetric" / str(order); out.mkdir(parents=True, exist_ok=True)
    write_certificate(cert, out / "certificate.json")
    (out / "sequences.txt").write_text("".join(s + "\n" for s in h["seqs"]))
    with gzip.open(out / f"H{order}.csv.gz", "wb", compresslevel=9) as f: f.write(data)
    m, digest = verify.check(out)
    print(f"symmetric order {m} (bordered propus over Z_{'x Z_'.join(map(str, group))}), rows {h['rows']}, sha256 {digest}")


if __name__ == "__main__":
    main()
