"""Exact quotient screen of a job (search/compress_check.py): the four blocks are compressed to Z_d for every
prime-power divisor d of the group order and the compressions are joined exactly on the PAF conditions of the level;
then every pair of coprime levels is joined jointly (a member must realise all its compressions at once). An empty
join rules the job out, whatever the sizes of its lists. Levels whose state space is too large are skipped."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import compress_check as CCK


def quotient_kill(mode, N, gens, pat, maxstates=2_000_000, maxpairs=2e7):
    """(killed, detail) for a job of qhlib's modes; detail maps each level to its number of solutions (0 = ruled out,
    None = state space too large, "large" = join too large)"""
    try:
        res = CCK.check('bsz' if mode == 'bord' else 'ssxx', N, list(gens), tuple(int(x) for x in pat),
                        maxstates=maxstates, stop=True, joint=True, maxpairs=maxpairs)
    except MemoryError:
        return False, {'error': 'memory'}
    return any(v == 0 for v in res.values()), res
