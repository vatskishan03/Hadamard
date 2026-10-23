"""cheap EXACT necessary conditions on a level-1 instance (s,e of the four blocks) for the existence of t-signs (level 2), derived from (I) and (Q):
 (F1) (I): sum_X PAF_s(k) + q sum_X PAF_t(k) = 0 (k != 0), PAF_t(k) = sum over the c_X(k) pairs (j, j+k) of two non-defect positions of +-1:  R(k) := -sum_X PAF_s(k)/q must be an integer with
      R(k) = c(k) mod 2 and |R(k)| <= c(k)  (c(k) = sum_X c_X(k)).
 (F2) (Q): Q(k) = sum_X sum_j t(j) (s(j+k) - s(j-k)) (t = +-1 on non-defect j, 0 on defects) must vanish: parity of #{(X,j): j non-defect, s(j+k) - s(j-k) odd} must be even (k in repsQ).
 score = sum_k R(k)^2 / max(c(k),1): small score = more likely to complete (HEURISTIC ranking; the filters F1, F2 are exact)."""
import numpy as np


def prefilter(sh, S, E):
    S = S.astype(np.int64); E = E.astype(np.int64); M = S.shape[0]; q = sh.q
    nd = (E == 0).astype(np.int64)
    ok = np.ones(M, bool); score = np.zeros(M)
    for k in sh.repsPM:
        if k == 0: continue
        sh_ = sh.shift[k]
        pafs = (S * S[:, :, sh_]).sum(2).sum(1)
        c = (nd * nd[:, :, sh_]).sum(2).sum(1)
        div = (pafs % q == 0)
        R = -(pafs // q)
        ok &= div & ((R - c) % 2 == 0) & (np.abs(R) <= c)
        score += R.astype(float) ** 2 / np.maximum(c, 1)
    for k in sh.repsQ:
        shp = sh.shift[k]; shm = sh.shneg[k]
        d = (S[:, :, shp] - S[:, :, shm]) & 1
        par = (nd * d).sum(2).sum(1) & 1
        ok &= (par == 0)
    return ok, score
