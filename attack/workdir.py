"""Locations shared by the attack tools.

Code is found relative to this file; nothing is written inside the repository.
- WORK is $HADAMARD_WORK, or <tmpdir>/hadamard-work when it is unset. It holds hit logs, temporary list files and compiled helpers.
- bin(name) compiles attack/<name>.c into WORK/bin on first use (or when the source is newer) and returns the binary's path.
"""
import os
import subprocess
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent            # attack/
ROOT = HERE.parent                                # repository root
WORK = Path(os.environ.get('HADAMARD_WORK') or Path(tempfile.gettempdir()) / 'hadamard-work')
WORK.mkdir(parents=True, exist_ok=True)

_FLAGS = {'find4b': ['-O3', '-march=native', '-fopenmp'],
          'join3': ['-O3', '-march=native', '-fopenmp'],
          'join4x': ['-O3', '-march=native', '-fopenmp'],
          'kt4': ['-O3', '-march=native', '-fopenmp', '-lm'],
          'match4': ['-O3', '-march=native', '-fopenmp', '-lm'],
          'skew_gen': ['-O3', '-march=native', '-fopenmp'],
          'spin_ls': ['-O3', '-march=native'],
          'spin_ls_gls': ['-O3', '-march=native'],
          'spin_walk': ['-O3', '-march=native', '-fopenmp']}


def bin(name):
    src = HERE / f'{name}.c'
    out = WORK / 'bin' / name
    out.parent.mkdir(exist_ok=True)
    if not out.exists() or out.stat().st_mtime < src.stat().st_mtime:
        flags = _FLAGS.get(name, ['-O3', '-march=native'])
        libs = [f for f in flags if f.startswith('-l')]
        tmp = out.with_name(f'.{name}.{os.getpid()}')        # compile aside, then rename: a running copy keeps its inode
        subprocess.run(['cc', *[f for f in flags if not f.startswith('-l')], str(src), '-o', str(tmp), *libs], check=True)
        os.replace(tmp, out)
    return str(out)


def tmpdir():
    d = WORK / 'tmp'
    d.mkdir(exist_ok=True)
    return str(d)


def hits(name='skew_hits.jsonl'):
    return str(WORK / name)
