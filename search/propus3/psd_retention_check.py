"""Sampling estimate of the PSD retention (numpy, uniform samples) against the exact retention counted by the engine,
on the orders 532 and 876 (validates the estimate used for v = 239: B 0.82%, D 77.7%, A 45.5%)."""
import sys; sys.path.insert(0, '.')
import numpy as np
from p3struct import Structure
import p3run as P
rng = np.random.default_rng(4)
for (v, gens, rows, exact) in ((133, [26], (11, 5, 19), (396 / 966, 48186 / 2855320, 1466892 / 2368236)),
                               (219, [4], (3, 17, 17), (216 / 924, 99504 / 7488432, 2033376 / 7488432))):
    st = Structure(v, gens)
    out = []
    for blk, r, b in ((st.A, rows[0], 4 * v), (st.B, rows[1], 2 * v), (st.B, rows[2], 4 * v)):
        Y = blk.sample(r, 200000, rng)
        out.append(float((st.psd(blk.values(Y)) <= b + 0.5).all(1).mean()))
    print(f'v={v}: PSD retention sampled (A,B,D) = ({out[0]:.4f}, {out[1]:.4f}, {out[2]:.4f})  exact from the engine = ({exact[0]:.4f}, {exact[1]:.4f}, {exact[2]:.4f})')
