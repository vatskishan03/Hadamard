import itertools, numpy as np, sys
p=int(sys.argv[1]); n=8
qr={x*x%p for x in range(1,p)}
chi=np.array([0]+[1 if x in qr else -1 for x in range(1,p)])
target=np.array(sorted([-np.sqrt(13)]+[-np.sqrt(5)]*3+[np.sqrt(5)]*3+[np.sqrt(13)]))
cands=[]
rest=list(range(2,p))
combos=list(itertools.combinations(rest,6))
B=200000
hits=[]
for st in range(0,len(combos),B):
    sets=np.array([(0,1)+c for c in combos[st:st+B]])
    S=chi[(sets[:,None,:]-sets[:,:,None])%p].astype(float)   # S[a,b]=chi(x_b-x_a)
    ev=np.linalg.eigvalsh(-S)
    ok=np.all(np.abs(ev-target)<1e-6,axis=1)
    hits.extend(sets[ok].tolist())
print("p",p,"8-subsets through {0,1}:",len(combos),"with the template Seidel spectrum:",len(hits))
import pickle; pickle.dump(hits,open(f'k8hits_{p}.pkl','wb'))
