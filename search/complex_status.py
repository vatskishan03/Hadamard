#!/usr/bin/env python3
"""Status of complex (quaternary) Hadamard matrices of order 2N, N odd, N < NMAX: first constructed in this repository,
known (from search/ssxx_sweep.py's literature filter) and also constructed here, or open.

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
    new = sorted(n for n in ours if n < nmax and S.known_reason(n) is None)
    known = sorted(n for n in ours if n < nmax and S.known_reason(n) is not None)
    opn = [N for N in range(3, nmax, 2) if S.known_reason(N) is None and N not in ours]
    print(f"N < {nmax}: {len(new)} orders first constructed here, {len(known)} known orders also constructed here, "
          f"{len(opn)} open")
    print("first constructed here:", " ".join(str(2 * n) for n in new))
    print("known, also constructed here:", " ".join(str(2 * n) for n in known))
    print("open (order 2N):", " ".join(str(2 * n) for n in opn))


if __name__ == "__main__":
    main()
