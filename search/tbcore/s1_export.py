"""Turn the .npz of s1_prepare.py into the binary input of join4.

usage: python3 s1_export.py LISTS.npz LISTS.bin
Format: int32 nA nB nC nD m; int16 VA, VB, VC, VD (row-major); int32 T. Run ./join4 LISTS.bin LOG2_CHUNKS (join4.c is built by
tbio.build_join4, e.g. `python3 -c "import tbio; print(tbio.build_join4('.'))"`).
"""
import sys, numpy as np
import tbio
d=np.load(sys.argv[1]); out=sys.argv[2]
VA,VB,VC,VD,T=[d[k] for k in ("VA","VB","VC","VD","T")]
tbio.write_binary_lists(out, VA, VB, VC, VD, T)
print("exported",len(VA),len(VB),len(VC),len(VD),VA.shape[1])
