#!/usr/bin/env python3
"""Status of complex (quaternary) Hadamard matrices of order 2N, N odd, N < NMAX: known before (with the reason, from
search/ssxx_sweep.py's literature filter), constructed in this repository (complex-hadamard/C<2N>.json), constructed
on V. Zaitsev's page of October 2026 (search/data/complex_zaitsev_2026.txt), or open. The open list is compared
with the page's own list of unknown orders.

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
    ours = {n for n in ours if n < nmax}
    odd = range(3, nmax, 2)
    known = {N for N in odd if S.known_reason(N) is not None}
    conc = {N for N in odd if S.known_reason(N) is None and S.concurrent_reason(N) is not None}
    opn = [N for N in odd if N not in known | conc | ours]
    both = sorted(ours & conc)
    print(f"N < {nmax}: {len(ours)} orders constructed here, {len(both)} of them also on Zaitsev's page; "
          f"{len(conc - ours)} more on Zaitsev's page only; {len(opn)} open")
    print("constructed here:", " ".join(str(2 * n) for n in sorted(ours)))
    print("not also on Zaitsev's page:", " ".join(str(2 * n) for n in sorted(ours - conc)) or "none")
    print("on Zaitsev's page only:", " ".join(str(2 * n) for n in sorted(conc - ours)))
    print("open (order 2N):", " ".join(str(2 * n) for n in opn))
    page = {o // 2 for o, m in S.load_zaitsev().items() if m == "unknown" and o // 2 < nmax}
    if set(opn) == page:
        print(f"the open list is the page's list of unknown orders below {2 * nmax}")
    else:
        print("differences from the page's list of unknown orders: open here only",
              sorted(2 * n for n in set(opn) - page), "| unknown on the page only", sorted(2 * n for n in page - set(opn)))


if __name__ == "__main__":
    main()
