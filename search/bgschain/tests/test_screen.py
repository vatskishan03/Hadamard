"""The exact mod-p^k screen and the score never remove a solution: engine with screen_mods/score = brute force (planted instances, k = 4 and k = 8, non-split)."""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import numpy as np
import nonsplit_lib as NL
from test_small3 import run_engine, expand
from test_nonsplit import ker_z4, lv
allok = True
SCR = dict(screen_mods=[[2, 3], [2, 5], [3, 2], [5, 1], [7, 1]], score=True, leaflog=False)
def run(cfg, levels, split, bd, label):
    global allok
    sols, n = NL.brute_all(cfg, maxn=9000); bs = set(sols)
    c = dict(cfg, **SCR)
    found, ordered, complete, ch, agg, nb = run_engine(c, levels, split, bd, True)
    es = expand(ch, found)
    ok = es == bs and complete
    print(f"{label:50s} brute {len(bs):6d} engine+screen {len(es):6d} screened {agg.get('screened', 0):.0f} obstructed {agg.get('obstructed', 0):.0f} -> {'OK' if ok else 'MISMATCH'}", flush=True)
    allok &= ok
m23 = [4, 4, 3]; S23 = [[2, 0, 0], [0, 2, 0], [0, 0, 1]]; g23 = [NL.blockdiag(NL.T2, [2])]
ks = [ker_z4(2, 1, 'Z4'), ker_z4(2, 1, 'Zp'), ker_z4(2, 1, '2Z4')]
T2m = np.array(NL.T2); u1 = (np.eye(2, dtype=int) + 2 * T2m) % 4
nu = np.zeros((3, 3), dtype=int); nu[:2, :2] = u1; nu[2, 2] = 2
for seed in (31, 32, 33, 34, 35):
    cfg, _ = NL.planted_sym(m23, S23, g23, nu.tolist(), seed)
    run(cfg, lv(ks[0], ks[2]), lv(ks[1]), 1, f"Z4^2xZ3 symseed{seed} L=2 sp=1")
    run(cfg, lv(ks[1], ks[2]), [], 1, f"Z4^2xZ3 symseed{seed} L=2 (no split)")
m3 = [4, 4, 4]; S3 = [[2, 0, 0], [0, 2, 0], [0, 0, 2]]; g3 = [NL.T3]
for seed in (1, 2, 3):
    cfg, _ = NL.planted(m3, S3, g3, seed)
    run(cfg, lv(S3), [], 1, f"Z4^3 k8 planted seed{seed}")
for seed in (21, 22, 23):                                  # Z4^2 x Z3, trivial-ish: planted with larger solution sets
    cfg, _ = NL.planted(m23, S23, g23, seed)
    run(cfg, lv(ks[0], ks[2]), [], 1, f"Z4^2xZ3 planted seed{seed}")
print('SCREEN TESTS:', 'ALL OK' if allok else 'FAILURES'); sys.exit(0 if allok else 1)
