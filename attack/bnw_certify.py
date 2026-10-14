#!/usr/bin/env python3
"""Write the certificate of a Hadamard matrix of order 4tw from T-matrices of order t and a bordered Williamson-type
quadruple of order w (attack/bnw.py), after rebuilding and checking it with verify/verify.py.

    python3 attack/bnw_certify.py HITS.jsonl LINE T_FILE OUTDIR      (T_FILE: t-matrices/T<t>.txt, or "1", "3")
"""
import hashlib, json, sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "verify"))
import verify

SMALL_T = {"1": ["+", "0", "0", "0"], "3": ["+00", "0+0", "00+", "000"]}


def main():
    hits, line, tfile, outdir = sys.argv[1], int(sys.argv[2]), sys.argv[3], Path(sys.argv[4])
    rec = json.loads(Path(hits).read_text().splitlines()[line])
    T = SMALL_T[tfile] if tfile in SMALL_T else [ln.strip() for ln in Path(tfile).read_text().splitlines() if ln.strip()]
    t, w = len(T[0]), rec["N"]
    cert = {"order": 4 * t * w, "construction": "cooper-wallis-bordered", "t": t, "w": w,
            "bordered_williamson_group": f"Z_{w - 1}", "bordered_williamson_multipliers": rec["gens"],
            "t_matrices": T, "bordered_williamson": rec["sequences"]}
    H, _ = verify.BUILDERS[cert["construction"]](cert)
    n = H.shape[0]
    assert n == cert["order"] and verify.is_pm1(H) and verify.gram_is_scalar(H, n), "not Hadamard"
    sha = hashlib.sha256()
    for chunk in verify.csv_chunks(H):
        sha.update(chunk)
    cert["sha256"] = sha.hexdigest()
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "certificate.json").write_text(json.dumps(cert, indent=1) + "\n")
    print(f"order {n}: H H^T = {n} I verified; sha256 {cert['sha256']}")


if __name__ == "__main__":
    main()
