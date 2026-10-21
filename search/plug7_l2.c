// plug7_l2.c (search/NOTES.md section 82): fast exact level-2 completion for the plug-in class Z_p x F_7 (H = squares), unbordered, q = 3 mod 4.
// Same semantics as the l2.c of search/NOTES.md section 78 (program not included; instance = four X0 sequences u_X = x0 + 6 e, completions = number of t-quadruples with
//   sum_X PAF_{s}(d) + 7 PAF_t(d) = 0,  sum_X Q_X(d) = 0,  d = 1..h),  but
//   * fibres by Gray code with incremental PAF / Q / DFT updates (no O(p^2) per element),
//   * per-block structure-of-arrays, elements bucketed by hash class so the pair loops stream memory, AVX2 PSD pair filter.
// usage: plug7_l2 file [maxsol_print] [tlimit_per_instance]      output as l2.c ("INST ..." per instance, "SOL ..." lines)
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include <immintrin.h>

#define MAXP 31
#define PS 32   // psd stride (floats), >= p, multiple of 8
static int p, h, T;
static double cs[MAXP][MAXP], sn[MAXP][MAXP];
static uint64_t rc[2 * MAXP];
static double now(void) { struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return t.tv_sec + 1e-9 * t.tv_nsec; }
static uint64_t rs = 0x9E3779B97F4A7C15ULL;
static uint64_t rnd(void) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; }

static int X0[4][MAXP], R[4][MAXP], Uu[4][MAXP];
static int vpos[4][MAXP], nv[4];
static int allowE[4][2 * MAXP + 1];
static double Rr[4][MAXP], Ri[4][MAXP];   // DFT of R at t=0..p-1
static int pafR[4][MAXP];                 // PAF_R(d), d = 1..h
static int wR[4][MAXP][MAXP];             // wR[b][j][d] = R[j+d] - R[j-d]

typedef struct { long n; uint64_t *hs; int16_t *key; uint8_t *psd; uint32_t *mask; float fmin[PS]; } Fib;

static void fib_free(Fib *F) { free(F->hs); free(F->key); free(F->psd); free(F->mask); memset(F, 0, sizeof(*F)); }

// build fibre of block b by Gray code. fixfirst: top bit fixed.
static void build_fibre(int b, Fib *F, int fixfirst) {
    int sym = (b < 2), n = nv[b];
    int nb = (fixfirst && n > 0) ? n - 1 : n;           // number of enumerated bits
    long N = 1L << nb;
    F->n = 0; for (int t = 0; t < PS; t++) F->fmin[t] = 1e30f;
    F->hs = malloc(sizeof(uint64_t) * N); F->key = malloc(sizeof(int16_t) * 2 * h * N);
    F->psd = aligned_alloc(32, PS * N); F->mask = malloc(sizeof(uint32_t) * N);
    if (!F->hs || !F->key || !F->psd || !F->mask) { fprintf(stderr, "oom fibre %ld\n", N); exit(3); }
    int E[MAXP]; memset(E, 0, sizeof(E));
    for (int j = 0; j < n; j++) { int z = vpos[b][j]; E[z] = 1; if (sym) E[p - z] = -1; }
    // initial features
    int P[MAXP], Qd[MAXP]; double ea[MAXP], eb[MAXP]; int tau = 0;
    for (int d = 1; d <= h; d++) {
        int pe = 0, im = 0;
        for (int z = 0; z < p; z++) { int z2 = z + d; if (z2 >= p) z2 -= p; pe += E[z] * E[z2]; im += E[z] * R[b][z2] - R[b][z] * E[z2]; }
        P[d] = pe; Qd[d] = im;
    }
    for (int t = 0; t <= h; t++) { double a = 0, bb = 0; for (int z = 0; z < p; z++) { a += E[z] * cs[t][z]; bb += E[z] * sn[t][z]; } ea[t] = a; eb[t] = bb; }
    for (int z = 0; z < p; z++) tau += E[z];
    const double s7 = sqrt(7.0);
    uint32_t mask = 0;
    for (long i = 0; i < N; i++) {
        if (i > 0) {
            int bit = __builtin_ctzl(i);
            mask ^= (1u << bit);
            int z0 = vpos[b][bit];
            int nf = sym ? 2 : 1; int zs[2] = {z0, p - z0};
            for (int f = 0; f < nf; f++) {
                int z = zs[f]; int ev = E[z];
                for (int d = 1; d <= h; d++) {
                    int zp = z + d; if (zp >= p) zp -= p; int zm = z - d; if (zm < 0) zm += p;
                    P[d] -= 2 * ev * (E[zp] + E[zm]);
                    Qd[d] -= 2 * ev * wR[b][z][d];
                }
                for (int t = 0; t <= h; t++) { ea[t] -= 2.0 * ev * cs[t][z]; eb[t] -= 2.0 * ev * sn[t][z]; }
                tau -= 2 * ev; E[z] = -ev;
            }
        }
        if (!allowE[b][tau + MAXP]) continue;
        // psd check
        uint8_t *ps = F->psd + (size_t)F->n * PS;
        int ok = 1; float vv[MAXP];
        for (int t = 0; t < p; t++) {
            double a, bb; if (t <= h) { a = ea[t]; bb = eb[t]; } else { a = ea[p - t]; bb = -eb[p - t]; }
            double re = Rr[b][t] - s7 * bb, im = Ri[b][t] + s7 * a;
            double v = re * re + im * im;
            if (v > T + 0.01) { ok = 0; break; }
            vv[t] = (float)v;
        }
        if (!ok) continue;
        for (int t = 0; t < p; t++) { ps[t] = (uint8_t)(vv[t] / 3.0f); if (vv[t] < F->fmin[t]) F->fmin[t] = vv[t]; }
        for (int t = p; t < PS; t++) ps[t] = 0;
        int16_t *k = F->key + (size_t)F->n * 2 * h; uint64_t hs = 0;
        for (int d = 1; d <= h; d++) { k[d - 1] = (int16_t)(pafR[b][d] + 7 * P[d]); k[h + d - 1] = (int16_t)Qd[d]; }
        for (int j = 0; j < 2 * h; j++) hs += rc[j] * (uint64_t)(int64_t)k[j];
        F->hs[F->n] = hs; F->mask[F->n] = mask; F->n++;
    }
}

typedef struct { uint64_t hs; uint32_t i, j; } Ent;

typedef struct { long n; uint64_t *hs; uint8_t *psd; uint32_t *orig; long *bs; } Sorted;   // elements sorted by class

static void sort_by_class(const Fib *F, uint64_t mm, uint64_t M, Sorted *S) {
    S->n = F->n; S->hs = malloc(sizeof(uint64_t) * (F->n + 1)); S->psd = aligned_alloc(32, PS * (F->n + 1)); S->orig = malloc(sizeof(uint32_t) * (F->n + 1));
    S->bs = calloc(M + 1, sizeof(long));
    for (long i = 0; i < F->n; i++) S->bs[(F->hs[i] & mm) + 1]++;
    for (uint64_t r = 0; r < M; r++) S->bs[r + 1] += S->bs[r];
    long *fill = malloc(sizeof(long) * M); memcpy(fill, S->bs, sizeof(long) * M);
    for (long i = 0; i < F->n; i++) {
        long q = fill[F->hs[i] & mm]++;
        S->hs[q] = F->hs[i]; S->orig[q] = (uint32_t)i; memcpy(S->psd + (size_t)q * PS, F->psd + (size_t)i * PS, PS);
    }
    free(fill);
}
static void sorted_free(Sorted *S) { free(S->hs); free(S->psd); free(S->orig); free(S->bs); }

static long long total_compl; static int nsol_print, maxsol_print;
static void emit_solution(int inst, const Fib *Fa, long ia, const Fib *Fb, long ib, const Fib *Fc, long ic, const Fib *Fd, long id) {
    const Fib *Fs[4] = {Fa, Fb, Fc, Fd}; long idx[4] = {ia, ib, ic, id};
    printf("SOL %d", inst);
    for (int k = 0; k < 4; k++) {
        int E[MAXP]; memset(E, 0, sizeof(E)); uint32_t msk = Fs[k]->mask[idx[k]];
        for (int j = 0; j < nv[k]; j++) { int z = vpos[k][j]; int e = ((msk >> j) & 1) ? -1 : 1; E[z] = e; if (k < 2) E[p - z] = -e; }
        printf(" |"); for (int z = 0; z < p; z++) printf(" %d:%d", X0[k][z], E[z]);
    }
    printf("\n");
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: plug7_l2 file [maxsol_print] [tlimit]\n"); return 1; }
    FILE *f = fopen(argv[1], "r"); if (!f) return 2;
    maxsol_print = argc > 2 ? atoi(argv[2]) : 5;
    double tlimit = argc > 3 ? atof(argv[3]) : 1e30;
    int ninst;
    if (fscanf(f, "%d %d", &p, &ninst) != 2) return 2;
    h = (p - 1) / 2; T = 28 * p - 2;
    for (int t = 0; t < p; t++) for (int z = 0; z < p; z++) { double a = 2 * M_PI * ((t * z) % p) / p; cs[t][z] = cos(a); sn[t][z] = sin(a); }
    for (int j = 0; j < 2 * h; j++) rc[j] = rnd() | 1;
    for (int inst = 0; inst < ninst; inst++) {
        for (int k = 0; k < 4; k++) for (int z = 0; z < p; z++) if (fscanf(f, "%d", &X0[k][z]) != 1) return 2;
        double t0 = now(); int bad = 0;
        for (int k = 0; k < 4; k++) {
            nv[k] = 0;
            for (int z = 0; z < p; z++) {
                int x = X0[k][z], c0 = 0, u = 0;
                switch (x) {
                    case 1: c0 = 1; u = 0; break; case -1: c0 = -1; u = 0; break;
                    case 5: c0 = -1; u = 1; break; case -5: c0 = 1; u = -1; break;
                    case 7: c0 = 1; u = 1; break; case -7: c0 = -1; u = -1; break;
                    case 6: c0 = 0; u = 1; break; case -6: c0 = 0; u = -1; break;
                    default: bad = 1;
                }
                R[k][z] = c0 - u; Uu[k][z] = u;
            }
            for (int z = (k < 2 ? 1 : 0); z < (k < 2 ? h + 1 : p); z++) if (Uu[k][z] == 0) vpos[k][nv[k]++] = z;
        }
        if (bad) { printf("INST %d bad input\n", inst); continue; }
        int fixk = nv[2] > 0 ? 2 : (nv[3] > 0 ? 3 : (nv[0] > 0 ? 0 : 1));
        // row-sum (level 0)
        {
            int sr2 = 0; for (int k = 0; k < 4; k++) { int sr = 0; for (int z = 0; z < p; z++) sr += R[k][z]; sr2 += sr * sr; }
            memset(allowE, 0, sizeof(allowE)); allowE[0][MAXP] = allowE[1][MAXP] = 1;
            int Q7 = T - sr2;
            if (Q7 >= 0 && Q7 % 7 == 0) {
                int Q = Q7 / 7;
                for (int x = -nv[2]; x <= nv[2]; x += 2) for (int y = -nv[3]; y <= nv[3]; y += 2) if (x * x + y * y == Q) { allowE[2][x + MAXP] = 1; allowE[3][y + MAXP] = 1; }
            }
        }
        for (int k = 0; k < 4; k++) {
            for (int t = 0; t < p; t++) { double a = 0, b = 0; for (int z = 0; z < p; z++) { a += R[k][z] * cs[t][z]; b += R[k][z] * sn[t][z]; } Rr[k][t] = a; Ri[k][t] = b; }
            for (int d = 1; d <= h; d++) { int a = 0; for (int z = 0; z < p; z++) a += R[k][z] * R[k][(z + d) % p]; pafR[k][d] = a; }
            for (int j = 0; j < p; j++) for (int d = 1; d <= h; d++) wR[k][j][d] = R[k][(j + d) % p] - R[k][(j - d + p) % p];
        }
        Fib F[4]; memset(F, 0, sizeof(F));
        for (int k = 0; k < 4; k++) build_fibre(k, &F[k], k == fixk);
        double tfib = now() - t0;
        // sides: (A,D) x (B,C)
        int sx = 0, sy = 3, sz = 1, sw = 2;
        double npairs1 = (double)F[sx].n * F[sy].n, npairs2 = (double)F[sz].n * F[sw].n;
        long long compl = 0, matches = 0, surv1 = 0, surv2 = 0; int timed_out = 0;
        int emptyf = 0; for (int k = 0; k < 4; k++) if (F[k].n == 0) emptyf = 1;
        if (!emptyf) {
            uint8_t lim1[PS] __attribute__((aligned(32))), lim2[PS] __attribute__((aligned(32)));
            for (int t = 0; t < PS; t++) {
                if (t < p) {
                    float l1 = T + 0.02f - F[sz].fmin[t] - F[sw].fmin[t], l2 = T + 0.02f - F[sx].fmin[t] - F[sy].fmin[t];
                    lim1[t] = (uint8_t)(l1 < 0 ? 0 : (l1 / 3.0f > 255 ? 255 : (int)(l1 / 3.0f)));
                    lim2[t] = (uint8_t)(l2 < 0 ? 0 : (l2 / 3.0f > 255 ? 255 : (int)(l2 / 3.0f)));
                } else { lim1[t] = lim2[t] = 255; }
            }
            double big = npairs1 > npairs2 ? npairs1 : npairs2;
            int lgM = 0; while ((big / (double)(1L << lgM)) > 20000.0 && lgM < 26) lgM++;
            uint64_t M = 1ull << lgM, mm = M - 1;
            Sorted SX, SY, SZ, SW;
            sort_by_class(&F[sx], mm, M, &SX); sort_by_class(&F[sy], mm, M, &SY); sort_by_class(&F[sz], mm, M, &SZ); sort_by_class(&F[sw], mm, M, &SW);
            size_t lcap = 1 << 16; Ent *list = malloc(sizeof(Ent) * lcap);
            size_t tcap = 1 << 17; Ent *tab = malloc(sizeof(Ent) * tcap);
            __m256i L1 = _mm256_load_si256((const __m256i *)lim1), L2 = _mm256_load_si256((const __m256i *)lim2); const __m256i zero = _mm256_setzero_si256();
            for (uint64_t r = 0; r < M && !timed_out; r++) {
                size_t nl = 0;
                for (long i = 0; i < SX.n; i++) {
                    // x elements are sorted by class too but we need all x for every r: loop over all x
                    uint64_t b = (r - SX.hs[i]) & mm;
                    long q0 = SY.bs[b], q1 = SY.bs[b + 1]; if (q0 == q1) continue;
                    __m256i xv = _mm256_load_si256((const __m256i *)(SX.psd + (size_t)i * PS));
                    for (long q = q0; q < q1; q++) {
                        __m256i sm = _mm256_adds_epu8(xv, _mm256_load_si256((const __m256i *)(SY.psd + (size_t)q * PS)));
                        if (_mm256_movemask_epi8(_mm256_cmpeq_epi8(_mm256_subs_epu8(sm, L1), zero)) != -1) continue;
                        if (nl >= lcap) { lcap *= 2; list = realloc(list, sizeof(Ent) * lcap); }
                        list[nl].hs = SX.hs[i] + SY.hs[q]; list[nl].i = (uint32_t)i; list[nl].j = (uint32_t)q; nl++;
                    }
                }
                surv1 += nl;
                if (!nl) continue;
                size_t need = 1; while (need < 2 * nl) need <<= 1;
                if (need > tcap) { tcap = need; tab = realloc(tab, sizeof(Ent) * tcap); }
                size_t cap = need; uint64_t cm = cap - 1;
                for (size_t q = 0; q < cap; q++) tab[q].i = 0xffffffffu;
                for (size_t q = 0; q < nl; q++) {
                    uint64_t pos = (list[q].hs * 0x9E3779B97F4A7C15ULL) >> 20 & cm;
                    while (tab[pos].i != 0xffffffffu) pos = (pos + 1) & cm;
                    tab[pos] = list[q];
                }
                uint64_t rr = (0 - r) & mm;
                for (long i = 0; i < SZ.n; i++) {
                    uint64_t b = (rr - SZ.hs[i]) & mm;
                    long q0 = SW.bs[b], q1 = SW.bs[b + 1]; if (q0 == q1) continue;
                    __m256i xv = _mm256_load_si256((const __m256i *)(SZ.psd + (size_t)i * PS));
                    for (long q = q0; q < q1; q++) {
                        __m256i sm = _mm256_adds_epu8(xv, _mm256_load_si256((const __m256i *)(SW.psd + (size_t)q * PS)));
                        if (_mm256_movemask_epi8(_mm256_cmpeq_epi8(_mm256_subs_epu8(sm, L2), zero)) != -1) continue;
                        surv2++;
                        uint64_t key = 0 - (SZ.hs[i] + SW.hs[q]);
                        uint64_t pos = (key * 0x9E3779B97F4A7C15ULL) >> 20 & cm;
                        while (tab[pos].i != 0xffffffffu) {
                            if (tab[pos].hs == key) {
                                matches++;
                                long xo = SX.orig[tab[pos].i], yo = SY.orig[tab[pos].j], zo = SZ.orig[i], wo = SW.orig[q];
                                const int16_t *kx = F[sx].key + xo * 2 * h, *ky = F[sy].key + yo * 2 * h, *kz = F[sz].key + zo * 2 * h, *kw = F[sw].key + wo * 2 * h;
                                int good = 1;
                                for (int j = 0; j < 2 * h; j++) if (kx[j] + ky[j] + kz[j] + kw[j] != 0) { good = 0; break; }
                                if (good) { compl++; if (nsol_print < maxsol_print) { nsol_print++; emit_solution(inst, &F[0], xo, &F[1], zo, &F[2], wo, &F[3], yo); } }
                            }
                            pos = (pos + 1) & cm;
                        }
                    }
                }
                if (now() - t0 > tlimit) timed_out = 1;
            }
            free(list); free(tab); sorted_free(&SX); sorted_free(&SY); sorted_free(&SZ); sorted_free(&SW);
        }
        long long mult = (nv[fixk] > 0) ? 2 : 1;
        total_compl += compl * mult;
        printf("INST %d nv=%d,%d,%d,%d fib=%ld,%ld,%ld,%ld pairs=%.3g,%.3g surv=%lld,%lld matches=%lld completions=%lld%s t=%.3f tfib=%.3f\n",
               inst, nv[0], nv[1], nv[2], nv[3], F[0].n, F[1].n, F[2].n, F[3].n, npairs1, npairs2, surv1, surv2, matches,
               compl * mult, timed_out ? " TIMEOUT" : "", now() - t0, tfib);
        fflush(stdout);
        for (int k = 0; k < 4; k++) fib_free(&F[k]);
    }
    printf("TOTAL completions %lld\n", total_compl);
    return 0;
}
