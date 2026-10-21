#!/usr/bin/env python3
"""Write complex-hadamard/C<2N>.json for the first verified hit of each N in a hit log of search/ssxx.py, after
rebuilding and checking the matrix with verify/verify_complex.py.

    python3 search/ssxx_certify.py HITS.jsonl [N ...]
"""
import hashlib, json, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "verify"))
import verify_complex as vc


def main():
    hits = [json.loads(l) for l in Path(sys.argv[1]).read_text().splitlines() if l.strip()]
    want = {int(a) for a in sys.argv[2:]}
    done = set()
    (ROOT / "complex-hadamard").mkdir(exist_ok=True)
    for h in hits:
        N = h["N"]
        if N in done or (want and N not in want):
            continue
        out = ROOT / "complex-hadamard" / f"C{2 * N}.json"
        if out.exists():
            sys.exit(f"{out.relative_to(ROOT)} exists; not overwritten")
        if "border" in h:
            cert = {"order": 2 * N, "construction": "szollosi-bordered-group" if "automorphisms" in h else
                    "szollosi-bordered", "N": N, "group": h.get("group", [N - 1]),
                    "symmetry": h.get("automorphisms", h.get("gens")), "rowsums": h["rowsums"], "border": h["border"],
                    "sequences": h["sequences"]}
        elif "group" in h:
            cert = {"order": 2 * N, "construction": "szollosi-ssxx-group", "N": N, "group": h["group"],
                    "field": h["field"], "generator": h["generator"], "subgroup_order": h["subgroup_order"],
                    "rowsums": h["rowsums"], "sequences": h["sequences"]}
        else:
            cert = {"order": 2 * N, "construction": "szollosi-ssxx", "N": N, "multipliers": h["gens"],
                    "subgroup_order": h["subgroup_order"], "rowsums": h["rowsums"], "sequences": h["sequences"]}
        Kr, Ki = vc.build(cert)
        vc.check(Kr, Ki)
        cert["sha256"] = hashlib.sha256(vc.text(Kr, Ki).encode()).hexdigest()
        out.write_text(json.dumps(cert, indent=1) + "\n")
        done.add(N)
        print(f"order {2 * N}: K K* = {2 * N} I verified; written {out.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
