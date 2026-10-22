// l1.c : C level-1 generator + join for PLAIN plug-in shapes (cyclic K = Z_n).  A, B: exact lists (given); C, D: LC random configs each of exact defect weight
// (uniform over (s,e) configs, canonical block sign), regenerated deterministically from (seed, index); join on the U features
//   feat[k] = sum_j u(j) u(j+k), u = s + q e, k in reps  (target given), via class-sorted hash join (A,C) x (B,D).
// Output: the instances (s,e of C and D, indices a,b of A and B).
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#define MAXN 64
typedef struct { uint64_t s; } Rng;
static inline uint64_t sm64(uint64_t *x) { uint64_t z = (*x += 0x9E3779B97F4A7C15ULL); z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL; z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL; return z ^ (z >> 31); }
static inline uint64_t xs(uint64_t *s) { uint64_t x = *s; x ^= x >> 12; x ^= x << 25; x ^= x >> 27; *s = x; return x * 0x2545F4914F6CDD1DULL; }

typedef struct {
    int n, q, nU; const int *reps; int nu; const int *w, *forced, *poff, *plist; const double *dpw; int nfree; const int *freeidx; int fw;
    uint64_t coef[64];
} Ctx;

// generate config j of the list (seed, weight nw): s,e arrays (int8, length n), canonical sign
static void gen_cfg(const Ctx *c, int nw, uint64_t seed, uint64_t j, int8_t *s, int8_t *e) {
    uint64_t st = seed ^ (j * 0xD1B54A32D192ED03ULL); uint64_t r0 = sm64(&st); uint64_t rs = r0 | 1;
    int n = c->n; int stride = n + 1;
    memset(s, 0, n); memset(e, 0, n);
    uint8_t defect[MAXN]; memset(defect, 0, sizeof(defect));
    for (int i = 0; i < c->nu; i++) if (c->forced[i]) defect[i] = 1;
    int rem = nw - c->fw;
    for (int ii = 0; ii < c->nfree; ii++) {
        int i = c->freeidx[ii]; int wi = c->w[i];
        double skip = 2.0 * c->dpw[(ii + 1) * stride + rem];
        double take = rem >= wi ? 4.0 * c->dpw[(ii + 1) * stride + rem - wi] : 0.0;
        double p = take / (take + skip + 1e-300);
        double u = (double)(xs(&rs) >> 11) * (1.0 / 9007199254740992.0);
        if (u < p) { defect[i] = 1; rem -= wi; }
    }
    int sg = 0;
    for (int i = 0; i < c->nu; i++) {
        uint64_t b = xs(&rs); int be = (b & 1) ? 1 : -1, bt = (b & 2) ? 1 : -1;
        int sv, ev;
        if (defect[i]) { ev = be; sv = (bt == 1) ? 0 : -2 * be; } else { ev = 0; sv = be; }
        for (int p = c->poff[i]; p < c->poff[i + 1]; p++) { s[c->plist[p]] = (int8_t)sv; e[c->plist[p]] = (int8_t)ev; }
        if (sg == 0) { int v = sv != 0 ? sv : ev; sg = v; }
    }
    if (sg < 0) for (int z = 0; z < n; z++) { s[z] = -s[z]; e[z] = -e[z]; }
}

static void feats(const Ctx *c, const int8_t *s, const int8_t *e, int16_t *f) {
    int n = c->n, q = c->q; int u[2 * MAXN];
    for (int z = 0; z < n; z++) { u[z] = s[z] + q * e[z]; u[z + n] = u[z]; }
    for (int i = 0; i < c->nU; i++) { int k = c->reps[i]; int a = 0; for (int z = 0; z < n; z++) a += u[z] * u[z + k]; f[i] = (int16_t)a; }
}
static inline uint64_t hsh(const Ctx *c, const int16_t *f) { uint64_t h = 0; for (int i = 0; i < c->nU; i++) h += c->coef[i] * (uint64_t)(int64_t)f[i]; return h; }

typedef struct { uint64_t key; uint32_t a, c; } ent;

long l1_run(int n, int q, int nU, const int *reps, const int16_t *target, int nu, const int *w, const int *forced, const int *poff, const int *plist,
            const double *dpw, int nfree, const int *freeidx, int fw, int nC, int nD, long LC, long LD, uint64_t seed,
            int nA, const int8_t *SA, const int8_t *EA, int nB, const int8_t *SB, const int8_t *EB, int lgM, long maxout,
            int8_t *oSC, int8_t *oEC, int8_t *oSD, int8_t *oED, uint32_t *oAB)
{
    Ctx cx = { n, q, nU, reps, nu, w, forced, poff, plist, dpw, nfree, freeidx, fw, {0} };
    uint64_t st = seed * 6364136223846793005ULL + 1; for (int i = 0; i < 64; i++) cx.coef[i] = sm64(&st) | 1;
    uint64_t M = 1ULL << lgM, mm = M - 1;
    int16_t f[64];
    uint64_t *hA = malloc(8L * nA), *hB = malloc(8L * nB);
    for (int i = 0; i < nA; i++) { feats(&cx, SA + (long)i * n, EA + (long)i * n, f); hA[i] = hsh(&cx, f); }
    for (int i = 0; i < nB; i++) { feats(&cx, SB + (long)i * n, EB + (long)i * n, f); hB[i] = hsh(&cx, f); }
    uint64_t ht = hsh(&cx, target);
    uint64_t seedC = seed * 2 + 1, seedD = seed * 2 + 2;
    // generate C, D hashes
    uint64_t *hC = malloc(8L * LC), *hD = malloc(8L * LD);
    int8_t s_[MAXN], e_[MAXN];
    for (long j = 0; j < LC; j++) { gen_cfg(&cx, nC, seedC, j, s_, e_); feats(&cx, s_, e_, f); hC[j] = hsh(&cx, f); }
    for (long j = 0; j < LD; j++) { gen_cfg(&cx, nD, seedD, j, s_, e_); feats(&cx, s_, e_, f); hD[j] = hsh(&cx, f); }
    long *bC = calloc(M + 1, sizeof(long)), *bD = calloc(M + 1, sizeof(long));
    for (long i = 0; i < LC; i++) bC[(hC[i] & mm) + 1]++;
    for (long i = 0; i < LD; i++) bD[(hD[i] & mm) + 1]++;
    for (uint64_t r = 0; r < M; r++) { bC[r + 1] += bC[r]; bD[r + 1] += bD[r]; }
    uint64_t *sC = malloc(8L * LC), *sD = malloc(8L * LD); uint32_t *iC = malloc(4L * LC), *iD = malloc(4L * LD);
    { long *fl = malloc(sizeof(long) * M); memcpy(fl, bC, sizeof(long) * M);
      for (long i = 0; i < LC; i++) { long p = fl[hC[i] & mm]++; sC[p] = hC[i]; iC[p] = (uint32_t)i; }
      memcpy(fl, bD, sizeof(long) * M);
      for (long i = 0; i < LD; i++) { long p = fl[hD[i] & mm]++; sD[p] = hD[i]; iD[p] = (uint32_t)i; } free(fl); }
    free(hC); free(hD);
    long found = 0; size_t cap = 1 << 10; ent *tab = malloc(sizeof(ent) * cap);
    size_t bmcap = 1 << 10; uint64_t *bm = malloc(8 * bmcap);
    int8_t sC_[MAXN], eC_[MAXN], sD_[MAXN], eD_[MAXN]; int16_t fa[64], fb[64], fc[64], fd[64];
    for (uint64_t r = 0; r < M; r++) {
        long npairs = 0;
        for (int a = 0; a < nA; a++) { uint64_t need = (r - hA[a]) & mm; npairs += bC[need + 1] - bC[need]; }
        if (!npairs) continue;
        size_t nc = 1; while (nc < (size_t)npairs * 2) nc <<= 1;
        if (nc > cap) { cap = nc; tab = realloc(tab, sizeof(ent) * cap); }
        uint64_t cm = nc - 1;
        for (size_t i = 0; i < nc; i++) tab[i].a = 0xffffffffu;
        // bitmap with 16 bits per pair (false-positive rate ~6%)
        int bmb = 4; while (((size_t)1 << bmb) < (size_t)npairs * 16 && bmb < 24) bmb++;
        size_t bmw = ((size_t)1 << bmb) >> 6; if (bmw < 1) bmw = 1; if (bmw > bmcap) { bmcap = bmw; bm = realloc(bm, 8 * bmcap); }
        memset(bm, 0, 8 * bmw); uint64_t bmm = ((uint64_t)1 << bmb) - 1;
        for (int a = 0; a < nA; a++) {
            uint64_t need = (r - hA[a]) & mm; uint64_t ha = hA[a];
            for (long p = bC[need]; p < bC[need + 1]; p++) {
                uint64_t key = ha + sC[p]; uint64_t i = ((key * 0x9E3779B97F4A7C15ULL) >> 24) & cm;
                { uint64_t bi = (key * 0xC2B2AE3D27D4EB4FULL) >> 40 & bmm; bm[bi >> 6] |= 1ULL << (bi & 63); }
                while (tab[i].a != 0xffffffffu) i = (i + 1) & cm;
                tab[i].key = key; tab[i].a = a; tab[i].c = iC[p];
            }
        }
        uint64_t rt = (ht - r) & mm;
        for (int b = 0; b < nB; b++) {
            uint64_t need = (rt - hB[b]) & mm; uint64_t hb = ht - hB[b];
            for (long p = bD[need]; p < bD[need + 1]; p++) {
                uint64_t key = hb - sD[p];
                { uint64_t bi = (key * 0xC2B2AE3D27D4EB4FULL) >> 40 & bmm; if (!((bm[bi >> 6] >> (bi & 63)) & 1)) continue; }
                uint64_t i = ((key * 0x9E3779B97F4A7C15ULL) >> 24) & cm;
                while (tab[i].a != 0xffffffffu) {
                    if (tab[i].key == key) {
                        uint32_t a = tab[i].a, cj = tab[i].c, dj = iD[p];
                        gen_cfg(&cx, nC, seedC, cj, sC_, eC_); gen_cfg(&cx, nD, seedD, dj, sD_, eD_);
                        feats(&cx, SA + (long)a * n, EA + (long)a * n, fa); feats(&cx, SB + (long)b * n, EB + (long)b * n, fb); feats(&cx, sC_, eC_, fc); feats(&cx, sD_, eD_, fd);
                        int ok = 1; for (int k = 0; k < nU; k++) if ((int)fa[k] + fb[k] + fc[k] + fd[k] != target[k]) { ok = 0; break; }
                        if (ok) {
                            if (found < maxout) { memcpy(oSC + found * n, sC_, n); memcpy(oEC + found * n, eC_, n); memcpy(oSD + found * n, sD_, n); memcpy(oED + found * n, eD_, n); oAB[2 * found] = a; oAB[2 * found + 1] = b; }
                            found++;
                        }
                    }
                    i = (i + 1) & cm;
                }
            }
        }
    }
    free(bm); free(tab); free(hA); free(hB); free(bC); free(bD); free(sC); free(sD); free(iC); free(iD);
    return found;
}

// regenerate one config (for tests)
void l1_gen(int n, int q, int nu, const int *w, const int *forced, const int *poff, const int *plist, const double *dpw, int nfree, const int *freeidx, int fw,
            int nw, uint64_t seed, uint64_t j, int8_t *s, int8_t *e)
{
    Ctx cx = { n, q, 0, NULL, nu, w, forced, poff, plist, dpw, nfree, freeidx, fw, {0} };
    gen_cfg(&cx, nw, seed, j, s, e);
}
