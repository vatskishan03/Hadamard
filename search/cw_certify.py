#!/usr/bin/env python3
"""Write the certificate of a Hadamard matrix of order 4tw by the Cooper-Wallis construction, after rebuilding and
checking it with verify/verify.py: T-matrices of order t give OD(4t; t, t, t, t), and Williamson matrices of order w
(four symmetric circulants with sum of squares 4wI), or bordered Williamson-type matrices, are substituted for its
variables.

    python3 search/cw_certify.py T_FILE W OUTDIR
        T_FILE  t-matrices/T<t>.txt (four rows over 0, +, -)
        W       an order in search/data/williamson_sage.json (Williamson matrices of SageMath's database),
                or williamson-type/WT<w>.json (bordered Williamson-type matrices)
"""
import hashlib, json, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "verify"))
import verify


def main():
    tfile, wspec, outdir = sys.argv[1], sys.argv[2], Path(sys.argv[3])
    T = [ln.strip() for ln in Path(tfile).read_text().splitlines() if ln.strip()]
    t = len(T[0])
    if wspec.isdigit():
        W = json.loads((ROOT / "search" / "data" / "williamson_sage.json").read_text())["quadruples"][wspec]
        w = len(W[0])
        cert = {"order": 4 * t * w, "construction": "cooper-wallis", "t": t, "w": w, "t_matrices": T,
                "williamson": list(W)}
    else:
        rec = json.loads(Path(wspec).read_text())
        w = int(rec["order"])
        cert = {"order": 4 * t * w, "construction": "cooper-wallis-bordered", "t": t, "w": w,
                "bordered_williamson_group": rec["group"], "bordered_williamson_multipliers": rec["multipliers"],
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
