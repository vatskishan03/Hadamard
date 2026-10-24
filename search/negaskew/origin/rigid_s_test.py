import sys, itertools, time
sys.argv=['x','11','47']
src=open('rigid_s.py').read().split("# A origin variants")[0]
exec(src)
# planted test: pick random C0 (any), D0 canonical; target = -(P_C0 + P_D0) in the sense need: P_C + P_D = P_C0 + P_D0
rng=np.random.default_rng(1)
ok=0
for trial in range(5):
    ci=rng.integers(0,seqs.shape[0]); di=rng.integers(0,Dseq.shape[0])
    T2=2*(PC[ci]+PD[di])
    keyT=lin(T2[None,:])[0]
    need=keyT-keyD
    lo=np.searchsorted(keyC_sorted,need,side='left'); hi=np.searchsorted(keyC_sorted,need,side='right')
    nz=np.nonzero(hi-lo)[0]
    found=False
    for d in nz:
        for cc in order[lo[d]:hi[d]]:
            if np.max(np.abs(PC[cc]+PD[d]-PC[ci]-PD[di]))<1e-9: found=True
    ok+=found
print('planted-solution test: found',ok,'of 5')
