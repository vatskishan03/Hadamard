"""Run the C exact counter on a structure's full pools."""
import subprocess, struct, os, tempfile
import numpy as np
import count_exact as X
import workdir as W

def dump(st, path):
    v = st.G.v; pm = v % 4; target = (-4 - 4 * pm) // 4
    A, B = X.pool_vectors(st, 2), X.pool_vectors(st, 0)
    VA = np.rint((st.class_values(A) - pm) / 4).astype(np.int16)
    VB = np.rint((st.class_values(B) - pm) / 4).astype(np.int16)
    with open(path, "wb") as f:
        f.write(struct.pack("iii", len(A), len(B), st.C))
        f.write(np.full(st.C, target, dtype=np.int16).tobytes())
        f.write(np.ascontiguousarray(VA).tobytes()); f.write(np.ascontiguousarray(VB).tobytes())
    return len(A), len(B)

def count(st):
    fd, path = tempfile.mkstemp(suffix=".bin", dir=W.tmpdir())
    os.close(fd)
    try:
        nA, nB = dump(st, path)
        out = subprocess.run([W.bin("count4"), path],
                             capture_output=True, text=True, check=True).stdout
        return int(out.strip()), nA, nB
    finally:
        os.remove(path)
