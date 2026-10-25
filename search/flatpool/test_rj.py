import numpy as np, subprocess, sys
C=27; n=3000; T=-3
rng=np.random.default_rng(5)
Ks=[np.clip(np.rint(rng.normal(0,2,size=(n,C))),-9,9).astype(np.int8) for _ in range(4)]
# plant: A[7], B[11], C[13] random; D[17] = T - sum
Ks[3][17]=(T-Ks[0][7].astype(int)-Ks[1][11]-Ks[2][13]).astype(np.int8)
print("planted D key range",Ks[3][17].min(),Ks[3][17].max())
with open('w/rjtest.bin','wb') as f:
    np.array([C,T,1,1,1,1,n,n,n,n],dtype=np.int32).tofile(f)
    for K in Ks: K.tofile(f)
