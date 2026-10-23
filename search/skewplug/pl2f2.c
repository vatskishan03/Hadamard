// pl2f.c : FAST exact level-2 (t-sign) completion for PLAIN Paley plug-in shapes (drop-in replacement of r31/sl1/pl2.c: same input, same output format, same set of completions).
// Differences to pl2.c (all exact / sound; the solution set is mathematically defined by the key equations, see "EXACTNESS" below):
//  * the pair filter is a transposed (structure-of-arrays) SIMD kernel: for a fixed outer element x, 32 inner elements are tested at once, channel by channel, with the exact same byte quantisation
//    and limits as pl2.c (p_x + p_y <= L per channel, L from the minima of the other pair of fibres);  elements x with p_x[c] > L[c] for some channel are dropped (sound);
//  * the join of the two pair lists uses a per-class table of 32-bit tags of the 64-bit linear hash; hash-class survivors are not stored with their indices: a tag hit triggers an exact
//    resolution by rescanning the class (rare: ~ the number of solutions) followed by the exact integer verification of all nk key components from the t-vectors;
//  * (fast fibre builder for blocks with singleton orbits and contiguous lags: see build_fibre_fast).
// EXACTNESS: a completion = a quadruple (one fibre element per block) whose nI + nQ key components (pafS + q PAF_t, and Q) sum to ZERO.  Every filter used before the final exact check is a NECESSARY condition
//  of such a quadruple (tau-set filters, PSD <= T + 0.01 per channel, pair PSD limit with floor quantisation: floor(a)+floor(b) <= floor(a+b)).  The final check recomputes the keys from the t-vectors in exact
//  integer arithmetic (elem_key).  Hence the set of completions printed/counted equals the one of pl2.c (if pl2.c is correct), independent of hashing, ordering and tags (false tag hits are filtered by the exact check).
// input/output format: see pl2.c.   usage: pl2f file [maxsol_print] [tlimit_per_instance] [seed]
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include <immintrin.h>
#include <malloc.h>

#define MAXN 64
#define MAXK 64
#define MAXV 40
#ifndef PS
#define PS 64
#endif
#define LANES 32
static int n, q, eps, T, nI, nQ, kI[MAXK], kQ[MAXK], nk, nch;
static int nu[2];
static int omm[2][MAXN], opos[2][MAXN][MAXN], osg[2][MAXN][MAXN];
static double cgr[MAXN][MAXN], cgi[MAXN][MAXN];     // cg[c][z] = gamma_c * omega^(psi_c z)
static double sre[MAXN][MAXN], sim_[MAXN][MAXN];     // omega^(psi z) for psi = 1..n-1 (index psi)
static int chpsi[MAXN];
static double sqq;
static uint64_t rc[2 * MAXK];
static double now(void) { struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return t.tv_sec + 1e-9 * t.tv_nsec; }
static uint64_t rs = 0x9E3779B97F4A7C15ULL;
static uint64_t rnd(void) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; }

static int S[4][MAXN], Dfl[4][MAXN];
static int nvar[4], vorb[4][MAXV];        // variable -> orbit index
static int allowT[4][2 * MAXN + 1];       // tau allowed (index tau + MAXN)
static int wQ[4][MAXN][MAXK];             // wQ[b][z][j] = s(z+k_j) + eps s(z-k_j)
static int DV;      // psd byte divisor

typedef struct { long n; long cap; uint64_t *hs; uint8_t *psd; uint32_t *mask; float fmin[PS]; } Fib;
static void fib_free(Fib *F) { free(F->hs); free(F->psd); free(F->mask); memset(F, 0, sizeof(*F)); }
static void fib_reserve(Fib *F, long need) {
    if (need <= F->cap) return;
    long nc = F->cap ? F->cap : 1024; while (nc < need) nc *= 2;
    F->hs = realloc(F->hs, sizeof(uint64_t) * nc);
    uint8_t *np = aligned_alloc(32, (size_t)PS * nc); if (F->psd) { memcpy(np, F->psd, (size_t)PS * F->n); free(F->psd); } F->psd = np;
    F->mask = realloc(F->mask, sizeof(uint32_t) * nc);
    if (!F->hs || !F->psd || !F->mask) { fprintf(stderr, "oom fibre %ld\n", nc); exit(3); }
    F->cap = nc;
}

// exact key of a fibre element of block b with variable-sign mask msk (bit j <-> variable j; the fixed variable (last of the fixed block) is +1 = bit not in the mask)
static void elem_t(int b, uint32_t msk, int *t) {
    int kd = b == 0 ? 0 : 1; memset(t, 0, sizeof(int) * MAXN);
    for (int j = 0; j < nvar[b]; j++) { int o = vorb[b][j]; int e = (j < 32 && ((msk >> j) & 1)) ? -1 : 1; for (int m = 0; m < omm[kd][o]; m++) t[opos[kd][o][m]] = e * osg[kd][o][m]; }
}
static void elem_key(int b, uint32_t msk, long *key) {
    int t[MAXN]; elem_t(b, msk, t); const int *s = S[b];
    for (int i = 0; i < nI; i++) {
        int k = kI[i], pe = 0, ps_ = 0;
        for (int z = 0; z < n; z++) { int z2 = (z + k) % n; pe += t[z] * t[z2]; ps_ += s[z] * s[z2]; }
        key[i] = (long)ps_ + (long)q * pe;
    }
    for (int j = 0; j < nQ; j++) { long a = 0; for (int z = 0; z < n; z++) a += (long)t[z] * wQ[b][z][j]; key[nI + j] = a; }
}

// ---------------------------------------------------------------- fibre builders
// generic builder (orbits of any size, any lags): Gray code over the variables with incremental updates (as pl2.c)
static void build_fibre(int b, Fib *F, int fixfirst) {
    int kd = b == 0 ? 0 : 1; int nvb = nvar[b];
    int nbits = (fixfirst && nvb > 0) ? nvb - 1 : nvb;
    long N = 1L << nbits;
    F->n = 0; for (int t = 0; t < PS; t++) F->fmin[t] = 1e30f;
    fib_reserve(F, 1024);
    int t[MAXN]; memset(t, 0, sizeof(t));
    for (int j = 0; j < nvb; j++) { int o = vorb[b][j]; for (int m = 0; m < omm[kd][o]; m++) t[opos[kd][o][m]] = osg[kd][o][m]; }
    int *s = S[b];
    int P[MAXK], Qd[MAXK], pafS[MAXK]; int tau = 0;
    for (int i = 0; i < nI; i++) {
        int k = kI[i], pe = 0, ps_ = 0;
        for (int z = 0; z < n; z++) { int z2 = (z + k) % n; pe += t[z] * t[z2]; ps_ += s[z] * s[z2]; }
        P[i] = pe; pafS[i] = ps_;
    }
    for (int j = 0; j < nQ; j++) { int a = 0; for (int z = 0; z < n; z++) a += t[z] * wQ[b][z][j]; Qd[j] = a; }
    double ar[MAXN], ai[MAXN];
    for (int c = 0; c < nch; c++) {
        int psi = chpsi[c]; double re = 0, im = 0;
        for (int z = 0; z < n; z++) { re += s[z] * sre[psi][z] + t[z] * cgr[c][z]; im += s[z] * sim_[psi][z] + t[z] * cgi[c][z]; }
        ar[c] = re; ai[c] = im;
    }
    for (int z = 0; z < n; z++) tau += t[z];
    uint32_t mask = 0;
    for (long i = 0; i < N; i++) {
        if (i > 0) {
            int bit = __builtin_ctzl(i);
            mask ^= (1u << bit);
            int o = vorb[b][bit];
            for (int m = 0; m < omm[kd][o]; m++) {
                int z = opos[kd][o][m]; int tv = t[z];
                for (int a = 0; a < nI; a++) { int k = kI[a]; int zp = z + k; if (zp >= n) zp -= n; int zm = z - k; if (zm < 0) zm += n; P[a] -= 2 * tv * (t[zp] + t[zm]); }
                for (int j = 0; j < nQ; j++) Qd[j] -= 2 * tv * wQ[b][z][j];
                for (int c = 0; c < nch; c++) { ar[c] -= 2.0 * tv * cgr[c][z]; ai[c] -= 2.0 * tv * cgi[c][z]; }
                tau -= 2 * tv; t[z] = -tv;
            }
        }
        if (!allowT[b][tau + MAXN]) continue;
        int ok = 1; float vv[MAXN];
        for (int c = 0; c < nch; c++) { double v = ar[c] * ar[c] + ai[c] * ai[c]; if (v > T + 0.01) { ok = 0; break; } vv[c] = (float)v; }
        if (!ok) continue;
        if (F->n >= F->cap) fib_reserve(F, F->n + 1);
        uint8_t *ps = F->psd + (size_t)F->n * PS;
        for (int c = 0; c < nch; c++) { ps[c] = (uint8_t)(vv[c] / DV); if (vv[c] < F->fmin[c]) F->fmin[c] = vv[c]; }
        for (int c = nch; c < PS; c++) ps[c] = 0;
        uint64_t hs = 0;
        for (int a = 0; a < nI; a++) hs += rc[a] * (uint64_t)(int64_t)(int16_t)(pafS[a] + q * P[a]);
        for (int j = 0; j < nQ; j++) hs += rc[nI + j] * (uint64_t)(int64_t)(int16_t)Qd[j];
        F->hs[F->n] = hs; F->mask[F->n] = mask; F->n++;
    }
}

#include "pl2f2_fast.h"

// ---------------------------------------------------------------- join (v3)
// A "side" = one pair of fibres (outer list X = smaller fibre, inner list Y = larger fibre in natural order, stored transposed: SoA blocks of LANES elements).  All pairs (x,y) passing the SIMD PSD
// filter are recorded as 32-bit tags of h = hs_x + hs_y in per-class buffers (class = h & mm); the two sides are then joined class by class (double-bit bitmap + exact scan).
// Tag hits are resolved exactly (class rescan + integer verification of all key components).
typedef struct {
    int nblk; uint8_t *soa;            // [blk][c][LANES], padding lanes = 255
    uint64_t *hs;                      // [pos] (padded)
    uint32_t *ybs, *ysort;             // lazily built for resolution: class bucket starts and fibre indices sorted by class
} SoaList;
static void soa_free(SoaList *Y) { free(Y->soa); free(Y->hs); free(Y->ybs); free(Y->ysort); memset(Y, 0, sizeof(*Y)); }
static void build_soa(const Fib *F, SoaList *Y) {
    long nn = F->n; int nblk = (int)((nn + LANES - 1) / LANES); if (nblk == 0) nblk = 1; Y->nblk = nblk;
    size_t nl = (size_t)nblk * LANES;
    Y->soa = aligned_alloc(32, (size_t)nblk * nch * LANES + 32); memset(Y->soa, 255, (size_t)nblk * nch * LANES);
    Y->hs = calloc(nl + 2 * LANES, sizeof(uint64_t)); memcpy(Y->hs, F->hs, sizeof(uint64_t) * nn);
    for (long i = 0; i < nn; i++) {
        size_t blk = i / LANES, lane = i % LANES;
        const uint8_t *ps = F->psd + (size_t)i * PS; uint8_t *dst = Y->soa + blk * (size_t)nch * LANES + lane;
        for (int c = 0; c < nch; c++) dst[(size_t)c * LANES] = ps[c];
    }
}
static void soa_sort(const Fib *F, SoaList *Y, uint64_t mm, uint64_t M) {
    if (Y->ybs) return;
    long nn = F->n; uint32_t *cnt = calloc(M + 1, sizeof(uint32_t));
    for (long i = 0; i < nn; i++) cnt[F->hs[i] & mm]++;
    Y->ybs = malloc(sizeof(uint32_t) * (M + 1)); uint32_t tot = 0; for (uint64_t b = 0; b < M; b++) { Y->ybs[b] = tot; tot += cnt[b]; } Y->ybs[M] = tot;
    Y->ysort = malloc(sizeof(uint32_t) * (nn + 1)); uint32_t *fill = malloc(sizeof(uint32_t) * M); memcpy(fill, Y->ybs, sizeof(uint32_t) * M);
    for (long i = 0; i < nn; i++) Y->ysort[fill[F->hs[i] & mm]++] = (uint32_t)i;
    free(cnt); free(fill);
}

typedef struct {
    int bX, bY; const Fib *FX, *FY;
    uint32_t nx; uint32_t *xo; uint64_t *xhs; uint8_t *xlim;   // alive outer elements
    SoaList Y;
} Side;
static void side_free(Side *sd) { free(sd->xo); free(sd->xhs); free(sd->xlim); soa_free(&sd->Y); memset(sd, 0, sizeof(*sd)); }
static void side_build(Side *sd, int b1, int b2, const Fib *F, const uint8_t *L) {
    int bx = b1, by = b2; if (F[b1].n > F[b2].n) { bx = b2; by = b1; }
    sd->bX = bx; sd->bY = by; sd->FX = &F[bx]; sd->FY = &F[by];
    const Fib *FX = sd->FX; long nn = FX->n;
    sd->xo = malloc(sizeof(uint32_t) * (nn + 1)); sd->xhs = malloc(sizeof(uint64_t) * (nn + 1));
    sd->xlim = aligned_alloc(32, (size_t)PS * (nn + 1)); sd->nx = 0;
    for (long i = 0; i < nn; i++) {
        const uint8_t *ps = FX->psd + (size_t)i * PS; uint8_t *lim = sd->xlim + (size_t)sd->nx * PS; int dead = 0;
        for (int c = 0; c < nch; c++) { int l = (int)L[c] - (int)ps[c]; if (l < 0) { dead = 1; break; } lim[c] = (uint8_t)(l > 254 ? 254 : l); }
        if (dead) continue;
        sd->xo[sd->nx] = (uint32_t)i; sd->xhs[sd->nx] = FX->hs[i]; sd->nx++;
    }
    build_soa(sd->FY, &sd->Y);
}

static inline uint32_t tag_of(uint64_t h) { return (uint32_t)(h >> 32) | 1u; }
typedef struct { uint32_t *p; uint32_t n, cap; } CB;
static void cb_grow(CB *c) { uint32_t nc = c->cap ? c->cap * 2 : 256; c->p = realloc(c->p, sizeof(uint32_t) * nc); if (!c->p) { fprintf(stderr, "oom cb\n"); exit(3); } c->cap = nc; }
static inline void cb_push(CB *c, uint32_t v) { if (__builtin_expect(c->n == c->cap, 0)) cb_grow(c); c->p[c->n++] = v; }
static void cbs_free(CB *cb, uint64_t M) { for (uint64_t r = 0; r < M; r++) free(cb[r].p); free(cb); }

static long long total_compl; static int nsol_print, maxsol_print;
static void emit_solution(int inst, const Fib *Fs[4], const long idx[4]) {
    printf("SOL %d", inst);
    for (int b = 0; b < 4; b++) {
        int t[MAXN]; elem_t(b, Fs[b]->mask[idx[b]], t);
        printf(" |"); for (int z = 0; z < n; z++) printf(" %d", t[z]);
    }
    printf("\n");
}

// SIMD pair filter: 'NB' consecutive blocks of the inner list against one outer element; masks m[i] bit j = lane j of block i passes
#define FILT_BODY(NB) \
    static inline void filt##NB(const uint8_t *p, const uint8_t *lim, uint32_t *m) { \
        const size_t bs = (size_t)nch * LANES; const __m256i zero = _mm256_setzero_si256(); \
        __m256i a0 = zero, a1 = zero, a2 = zero, a3 = zero; (void)a1; (void)a2; (void)a3; \
        for (int c = 0; c < nch; c++) { \
            __m256i l = _mm256_set1_epi8((char)lim[c]); const uint8_t *pc = p + (size_t)c * LANES; \
            a0 = _mm256_or_si256(a0, _mm256_subs_epu8(_mm256_load_si256((const __m256i *)pc), l)); \
            if (NB > 1) a1 = _mm256_or_si256(a1, _mm256_subs_epu8(_mm256_load_si256((const __m256i *)(pc + bs)), l)); \
            if (NB > 2) a2 = _mm256_or_si256(a2, _mm256_subs_epu8(_mm256_load_si256((const __m256i *)(pc + 2 * bs)), l)); \
            if (NB > 3) a3 = _mm256_or_si256(a3, _mm256_subs_epu8(_mm256_load_si256((const __m256i *)(pc + 3 * bs)), l)); \
        } \
        m[0] = (uint32_t)_mm256_movemask_epi8(_mm256_cmpeq_epi8(a0, zero)); \
        if (NB > 1) m[1] = (uint32_t)_mm256_movemask_epi8(_mm256_cmpeq_epi8(a1, zero)); \
        if (NB > 2) m[2] = (uint32_t)_mm256_movemask_epi8(_mm256_cmpeq_epi8(a2, zero)); \
        if (NB > 3) m[3] = (uint32_t)_mm256_movemask_epi8(_mm256_cmpeq_epi8(a3, zero)); \
    }
FILT_BODY(1) FILT_BODY(2) FILT_BODY(3) FILT_BODY(4)

#define TILE 8
#define STG 512
static inline void flush_stg(const uint64_t *stg, long n_, uint64_t mm, int negate, CB *cb) {
    if (negate) { for (long k = 0; k < n_; k++) { uint64_t h = 0 - stg[k]; cb_push(&cb[h & mm], tag_of(h)); } }
    else { for (long k = 0; k < n_; k++) { uint64_t h = stg[k]; cb_push(&cb[h & mm], tag_of(h)); } }
}
// scan all (x, y) pairs of a side; survivors are pushed (as 32-bit tags of h or -h) into the per-class buffers cb (class = h & mm resp. (-h) & mm)
static long long scan_side(const Side *sd, uint64_t mm, int negate, CB *cb, double t0, double tlimit, int *timed_out) {
    long long ns = 0; const SoaList *Y = &sd->Y; const size_t bs = (size_t)nch * LANES;
    uint64_t stg[STG + TILE * LANES + 16]; long nst = 0;
    for (int blk0 = 0; blk0 < Y->nblk; blk0 += TILE) {
        int nbt = Y->nblk - blk0; if (nbt > TILE) nbt = TILE;
        for (uint32_t xi = 0; xi < sd->nx; xi++) {
            const uint8_t *lim = sd->xlim + (size_t)xi * PS; const uint64_t hx = sd->xhs[xi];
            for (int g = 0; g < nbt;) {
                uint32_t m[4]; int nb = nbt - g; if (nb > 4) nb = 4;
                const uint8_t *p = Y->soa + (size_t)(blk0 + g) * bs;
                if (nb == 4) filt4(p, lim, m); else if (nb == 3) filt3(p, lim, m); else if (nb == 2) filt2(p, lim, m); else filt1(p, lim, m);
                for (int i = 0; i < nb; i++) {
                    uint32_t mk = m[i]; const uint64_t *hy = Y->hs + (size_t)(blk0 + g + i) * LANES;
                    int cn = __builtin_popcount(mk);
                    uint32_t m1 = mk;
                    unsigned j0 = _tzcnt_u32(m1); m1 = _blsr_u32(m1); unsigned j1 = _tzcnt_u32(m1); m1 = _blsr_u32(m1);
                    unsigned j2 = _tzcnt_u32(m1); m1 = _blsr_u32(m1); unsigned j3 = _tzcnt_u32(m1); m1 = _blsr_u32(m1);
                    stg[nst] = hx + hy[j0]; stg[nst + 1] = hx + hy[j1]; stg[nst + 2] = hx + hy[j2]; stg[nst + 3] = hx + hy[j3];
                    if (__builtin_expect(cn <= 4, 1)) nst += cn;
                    else { nst += 4; while (m1) { unsigned j = _tzcnt_u32(m1); m1 = _blsr_u32(m1); stg[nst++] = hx + hy[j]; } }
                }
                g += nb;
            }
            if (nst >= STG) { flush_stg(stg, nst, mm, negate, cb); ns += nst; nst = 0; }
        }
        if ((blk0 & 63) == 0 && now() - t0 > tlimit) { *timed_out = 1; break; }
    }
    flush_stg(stg, nst, mm, negate, cb); ns += nst;
    return ns;
}

static long long compl_g, matches_g; static int cur_inst;
// verify the quadruple idx[4] exactly (all nk key components sum to zero)
static void check_quad(const Fib *F, const long idx[4]) {
    long key[4][2 * MAXK];
    for (int b = 0; b < 4; b++) elem_key(b, F[b].mask[idx[b]], key[b]);
    for (int j = 0; j < nk; j++) { long sum = key[0][j] + key[1][j] + key[2][j] + key[3][j]; if (sum != 0) return; }
    compl_g++;
    if (nsol_print < maxsol_print) { nsol_print++; const Fib *Fs[4] = {&F[0], &F[1], &F[2], &F[3]}; emit_solution(cur_inst, Fs, idx); }
}
typedef struct { uint32_t tag; long a, b; } PairRec;      // block-X index, block-Y index (fibre indices)
// all pairs (x,y) of the side with (h & mm) == cls (h = hs_x + hs_y, or -h if negate) and tag in hits[0..nh): appended to out
static void enum_class(Side *sd, uint64_t mm, uint64_t M, uint64_t cls, int negate, const uint32_t *hits, int nh, PairRec **out, long *nout, long *ocap) {
    soa_sort(sd->FY, &sd->Y, mm, M);
    const SoaList *Y = &sd->Y; const Fib *FX = sd->FX; const Fib *FY = sd->FY;
    uint64_t tgt = negate ? ((0 - cls) & mm) : cls;      // class of hs_x + hs_y
    for (long xi = 0; xi < FX->n; xi++) {
        uint64_t hx = FX->hs[xi]; uint32_t b = (uint32_t)((tgt - hx) & mm);
        for (uint32_t pos = Y->ybs[b]; pos < Y->ybs[b + 1]; pos++) {
            uint32_t yi = Y->ysort[pos];
            uint64_t h = hx + FY->hs[yi]; if (negate) h = 0 - h;
            uint32_t tg = tag_of(h); int lo = 0, hi = nh - 1, ok = 0;
            while (lo <= hi) { int mid = (lo + hi) >> 1; if (hits[mid] == tg) { ok = 1; break; } if (hits[mid] < tg) lo = mid + 1; else hi = mid - 1; }
            if (!ok) continue;
            if (*nout >= *ocap) { *ocap = *ocap ? *ocap * 2 : 64; *out = realloc(*out, sizeof(PairRec) * *ocap); }
            (*out)[*nout].tag = tg; (*out)[*nout].a = xi; (*out)[*nout].b = yi; (*nout)++;
        }
    }
}
static void resolve_class(Side *SL, Side *SR, const Fib *F, uint64_t mm, uint64_t M, uint64_t r, const uint32_t *hits, int nh) {
    PairRec *Lp = NULL, *Rp = NULL; long nl = 0, nr = 0, cl = 0, cr = 0;
    enum_class(SL, mm, M, r, 0, hits, nh, &Lp, &nl, &cl);
    enum_class(SR, mm, M, r, 1, hits, nh, &Rp, &nr, &cr);
    for (long i = 0; i < nl; i++) for (long j = 0; j < nr; j++) {
        if (Lp[i].tag != Rp[j].tag) continue;
        long idx[4]; idx[SL->bX] = Lp[i].a; idx[SL->bY] = Lp[i].b; idx[SR->bX] = Rp[j].a; idx[SR->bY] = Rp[j].b;
        matches_g++; check_quad(F, idx);
    }
    free(Lp); free(Rp);
}
static int cmp_u32(const void *a, const void *b) { uint32_t x = *(const uint32_t *)a, y = *(const uint32_t *)b; return x < y ? -1 : x > y; }

// join all classes: bitmap (2 bits per left tag) + exact vectorised scan on a bitmap hit; returns nothing, counts via compl_g/matches_g
static void join_classes(Side *SL, Side *SR, const Fib *F, uint64_t mm, uint64_t M, CB *cbL, CB *cbR, double t0, double tlimit, int *timed_out) {
    uint64_t *bm = NULL; size_t bmcap = 0; uint32_t *hitsb = NULL; int hcap = 0;
    for (uint64_t r = 0; r < M; r++) {
        uint32_t nL = cbL[r].n, nR = cbR[r].n; if (!nL || !nR) continue;
        int lgb = 12; while (lgb < 22 && ((size_t)1 << lgb) < 48 * (size_t)nL) lgb++;
        size_t words = ((size_t)1 << lgb) >> 6; if (words > bmcap) { free(bm); bm = malloc(sizeof(uint64_t) * words); bmcap = words; }
        memset(bm, 0, sizeof(uint64_t) * words); const int sh_ = 32 - lgb;
        const uint32_t *pl = cbL[r].p;
        for (uint32_t i = 0; i < nL; i++) {
            uint32_t tg = pl[i]; uint32_t i1 = (tg * 0x9E3779B1u) >> sh_, i2 = (tg * 0x85EBCA6Bu) >> sh_;
            bm[i1 >> 6] |= 1ull << (i1 & 63); bm[i2 >> 6] |= 1ull << (i2 & 63);
        }
        int nh = 0; const uint32_t *pr = cbR[r].p;
        for (uint32_t i = 0; i < nR; i++) {
            uint32_t tg = pr[i]; uint32_t i1 = (tg * 0x9E3779B1u) >> sh_;
            if (!((bm[i1 >> 6] >> (i1 & 63)) & 1)) continue;
            uint32_t i2 = (tg * 0x85EBCA6Bu) >> sh_;
            if (!((bm[i2 >> 6] >> (i2 & 63)) & 1)) continue;
            // exact scan of the left tags
            __m256i tv = _mm256_set1_epi32((int)tg); int found = 0; uint32_t k = 0;
            for (; k + 8 <= nL; k += 8) { if (_mm256_movemask_epi8(_mm256_cmpeq_epi32(_mm256_loadu_si256((const __m256i *)(pl + k)), tv))) { found = 1; break; } }
            if (!found) for (; k < nL; k++) if (pl[k] == tg) { found = 1; break; }
            if (!found) continue;
            if (nh >= hcap) { hcap = hcap ? hcap * 2 : 64; hitsb = realloc(hitsb, sizeof(uint32_t) * hcap); }
            hitsb[nh++] = tg;
        }
        if (nh > 1) { qsort(hitsb, nh, sizeof(uint32_t), cmp_u32); int w = 1; for (int i = 1; i < nh; i++) if (hitsb[i] != hitsb[w - 1]) hitsb[w++] = hitsb[i]; nh = w; }
        if (nh) resolve_class(SL, SR, F, mm, M, r, hitsb, nh);
        if ((r & 255) == 255 && now() - t0 > tlimit) { *timed_out = 1; break; }
    }
    free(bm); free(hitsb);
}

#include "pl2f2_main.h"
