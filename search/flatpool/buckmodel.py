#!/usr/bin/env python3
"""buckmodel.py rjinput.bin ERRFILE... : Gaussian prediction of solutions / pairs per processed bucket vs what rj found.
Model: key vectors of each (rotation-symmetrised) role Gaussian with covariance Sigma_i; the number of tuples with partial sum
sigma on S for pair (0,1) is L0 L1 N(sigma; mu0+mu1, S0+S1); solutions in the bucket = E_tuples * N_S(sigma; conditional law)."""
import sys, re, math, numpy as np
import est
fn = sys.argv[1]
hdr = np.fromfile(fn, dtype=np.int32, count=10); C, T = hdr[0], hdr[1]; rot = hdr[2:6]; n = hdr[6:10]
off = 40; Ks = []
with open(fn, 'rb') as f:
    f.seek(40)
    for i in range(4):
        K = np.fromfile(f, dtype=np.int8, count=int(n[i]) * C).reshape(-1, C).astype(float); Ks.append(K)
import pickle, os as _os
cf = fn + '.stats.pkl'
if _os.path.exists(cf):
    mus, Ss, Ls = pickle.load(open(cf, 'rb'))
else:
    mus, Ss, Ls = [], [], []
    for K, r in zip(Ks, rot):
        m, S, L = est.stats(K, bool(r)); mus.append(m); Ss.append(S); Ls.append(L)
    pickle.dump((mus, Ss, Ls), open(cf, 'wb'))
import os
NCc = int(os.environ.get('NC', C))
idx = list(range(NCc if NCc < C else C - 1))
mu = sum(mus); S = sum(Ss)
Si = np.linalg.inv(S[np.ix_(idx, idx)])
d = np.full(len(idx), float(T)) - mu[idx]
q = d @ Si @ d; sign, logdet = np.linalg.slogdet(S[np.ix_(idx, idx)])
lpdf = -0.5 * (len(idx) * math.log(2 * math.pi) + logdet + q)
E = math.exp(lpdf) * np.prod([float(x) for x in Ls])
print("E_tuples(model) = %.4g   (E_classes %.4g)" % (E, E / C))
SAB = (Ss[0] + Ss[1])[np.ix_(idx, idx)]; muAB = (mus[0] + mus[1])[idx]
G = SAB @ Si
cm = muAB + G @ d; cV = SAB - G @ SAB
buck = []
for fe in sys.argv[2:]:
    for line in open(fe):
        mm = re.match(r"BUCKET (\d+) sigma ([\-0-9 ]+) pairsAB (\d+) pairsCD (\d+) matches (\d+) found (\d+)", line)
        if mm: buck.append((int(mm.group(1)), [int(x) for x in mm.group(2).split()], int(mm.group(3)), int(mm.group(4)), int(mm.group(5)), int(mm.group(6))))
ns = len(buck[0][1]); sl = list(range(ns))
Vs = cV[np.ix_(sl, sl)]; Vi = np.linalg.inv(Vs); ldet = np.linalg.slogdet(Vs)[1]
SABs = SAB[np.ix_(sl, sl)]; SABi = np.linalg.inv(SABs); ldab = np.linalg.slogdet(SABs)[1]
tot_pred = 0; tot_found = 0; pr_ab = 0; act_ab = 0
for b, sig, pab, pcd, mt, fd in buck:
    x = np.array(sig, float)
    pred = E * math.exp(-0.5 * (ns * math.log(2 * math.pi) + ldet + (x - cm[sl]) @ Vi @ (x - cm[sl])))
    tot_pred += pred; tot_found += (fd if NCc >= C else mt)
    pa = float(Ls[0]) * float(Ls[1]) * math.exp(-0.5 * (ns * math.log(2 * math.pi) + ldab + (x - muAB[sl]) @ SABi @ (x - muAB[sl])))
    pr_ab += pa; act_ab += pab
print("buckets %d : predicted solutions %.3f, found %d ; predicted AB pairs %.4g actual %.4g (ratio %.3f)" % (len(buck), tot_pred, tot_found, pr_ab, act_ab, act_ab / pr_ab))
