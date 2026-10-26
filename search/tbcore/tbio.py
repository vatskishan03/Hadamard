"""Shared helpers of search/tbcore: the output directory, the tuple file, and the on-demand build of join4.c.

Nothing is written into the source tree: every program takes --out DIR (default: the current directory) and
writes its files there; the C helper is compiled into the same directory on first use.
"""
import hashlib
import json
import os
import subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
DATA = os.path.join(HERE, "data")


def outdir(path=None):
    """The output directory (created if needed); default the current directory."""
    path = os.path.abspath(path or os.getcwd())
    os.makedirs(path, exist_ok=True)
    return path


def tuples_path(out):
    return os.path.join(out, "tuples.json")


def trows_path(out):
    return os.path.join(out, "trows.npy")


def save_tuples(out, tuples):
    """tuples: list of (zA, zB, zC, zD, weight) with z* tuples of 5 ints and weight a Python int."""
    with open(tuples_path(out), "w") as f:
        json.dump([[list(t[0]), list(t[1]), list(t[2]), list(t[3]), int(t[4])] for t in tuples], f)


def load_tuples(out):
    with open(tuples_path(out)) as f:
        raw = json.load(f)
    return [(tuple(t[0]), tuple(t[1]), tuple(t[2]), tuple(t[3]), t[4]) for t in raw]


def tuples_digest(tuples):
    """SHA-256 of the canonical text of the tuple list (order matters: the index of a tuple is its rank by weight)."""
    s = json.dumps([[list(t[0]), list(t[1]), list(t[2]), list(t[3]), int(t[4])] for t in tuples],
                   sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(s.encode()).hexdigest()


def build_join4(out, name="join4"):
    """Compile join4.c into OUT/join4 (if missing or older than the source) and return the executable path."""
    src = os.path.join(HERE, name + ".c")
    exe = os.path.join(out, name)
    if not os.path.exists(exe) or os.path.getmtime(exe) < os.path.getmtime(src):
        tmp = exe + ".tmp%d" % os.getpid()
        subprocess.check_call(["gcc", "-O2", src, "-o", tmp])
        os.replace(tmp, exe)
    return exe


def write_binary_lists(path, VA, VB, VC, VD, T):
    """The input file of join4: int32 nA nB nC nD m; int16 VA, VB, VC, VD (row-major); int32 T."""
    import numpy as np
    m = VA.shape[1]
    with open(path, "wb") as f:
        f.write(np.array([len(VA), len(VB), len(VC), len(VD), m], dtype=np.int32).tobytes())
        for V in (VA, VB, VC, VD):
            f.write(np.ascontiguousarray(V, dtype=np.int16).tobytes())
        f.write(np.ascontiguousarray(T, dtype=np.int32).tobytes())

