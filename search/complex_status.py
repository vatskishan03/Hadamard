#!/usr/bin/env python3
"""Status of complex (quaternary) Hadamard matrices of order 2N, N odd, N < NMAX: known before (with the reason, from
search/ssxx_sweep.py's literature filter), constructed in this repository (complex-hadamard/C<2N>.json), or open.

    python3 search/complex_status.py [NMAX]
"""
import json, sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ssxx_sweep as S

ROOT = Path(__file__).resolve().parent.parent


def main():
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
    ours = {json.loads(p.read_text())["N"] for p in (ROOT / "complex-hadamard").glob("C*.json")}
    opn = [N for N in range(3, nmax, 2) if S.known_reason(N) is None and N not in ours]
    print(f"N < {nmax}: {len(ours & set(range(nmax)))} orders constructed here, {len(opn)} open")
    print("constructed here:", " ".join(str(2 * n) for n in sorted(ours) if n < nmax))
    print("open (order 2N):", " ".join(str(2 * n) for n in opn))


if __name__ == "__main__":
    main()
