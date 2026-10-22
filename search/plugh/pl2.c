// pl2.c : fast exact level-2 (t-sign) completion for PLAIN Paley plug-in shapes, cyclic K = Z_n, any odd q (prime power), any K-symmetric orbit structure.
// Generalisation of search/plug7_l2.c (NOTES sec 82) to the Hadamard (plain, no zero position) array: given the four (s, defect-flag) of a level-1
// instance, find all t-sign assignments (t = tau * sigma on non-defect orbits, 0 on defects) with
//    (I)  sum_X [PAF_s(k) + q PAF_t(k)] = 0   for k in repsI (k != 0)
//    (Q)  sum_X Q_X(k) = 0,  Q_X(k) = sum_j t(j) (s(j+k) + eps s(j-k))   for k in repsQ
// plus the exact necessary conditions used as filters:  row sums  sum_X (sigma_X^2 + q tau_X^2) = T [and sum sigma_X tau_X = 0 for eps = +1],
// PSD  |s^(psi) + gamma t^(psi)|^2 <= T  for every channel (gamma = i sqrt(q) if eps = -1, +-sqrt(q) if eps = +1).
// The two block pairs (A,D) x (B,C) are joined by hash-class bucketing with an AVX2 PSD pair filter (as plug7_l2).  Counts are of t-quadruples
// (x2 for the global t -> -t gauge).
// input file (whitespace separated):
//   n q eps T
//   nI k...            (I-equation lags, k != 0)
//   nQ k...            (Q-equation lags, may contain 0)
//   nuAB ; per orbit: m (pos sig)*m          (kind of blocks 0,1)
//   nuCD ; same                              (blocks 2,3)
//   ninst
//   per instance: s_0[0..n-1] s_1 s_2 s_3 ; then defect flags d_0[0..nuAB-1] d_1[..] d_2[0..nuCD-1] d_3[..]
// output: "INST i nv=.. fib=.. pairs=.. surv=.. matches=.. compl=.. t=.." and "SOL i | t_0 | t_1 | t_2 | t_3" (full length n each)
// usage: pl2 file [maxsol_print] [tlimit_per_instance] [seed]
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include <immintrin.h>

#define MAXN 64
#define MAXK 64
#define MAXV 40
#ifndef PS
#define PS 64
#endif
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
static int ptb[4][MAXN];                  // t array initial (tau=+1)

typedef struct { long n; uint64_t *hs; int16_t *key; uint8_t *psd; uint32_t *mask; float fmin[PS]; } Fib;
static void fib_free(Fib *F) { free(F->hs); free(F->key); free(F->psd); free(F->mask); memset(F, 0, sizeof(*F)); }
static int DV;      // psd byte divisor

static void build_fibre(int b, Fib *F, int fixfirst) {
    int kd = b < 2 ? 0 : 1; int nvb = nvar[b];
    int nbits = (fixfirst && nvb > 0) ? nvb - 1 : nvb;
    long N = 1L << nbits;
    F->n = 0; for (int t = 0; t < PS; t++) F->fmin[t] = 1e30f;
    F->hs = malloc(sizeof(uint64_t) * N); F->key = malloc(sizeof(int16_t) * nk * N);
    F->psd = aligned_alloc(32, (size_t)PS * N); F->mask = malloc(sizeof(uint32_t) * N);
    if (!F->hs || !F->key || !F->psd || !F->mask) { fprintf(stderr, "oom fibre %ld\n", N); exit(3); }
    int t[MAXN]; memset(t, 0, sizeof(t));
    for (int j = 0; j < nvb; j++) { int o = vorb[b][j]; for (int m = 0; m < omm[kd][o]; m++) t[opos[kd][o][m]] = osg[kd][o][m]; }
    int *s = S[b];
    // initial features
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
        uint8_t *ps = F->psd + (size_t)F->n * PS;
        int ok = 1; float vv[MAXN];
        for (int c = 0; c < nch; c++) { double v = ar[c] * ar[c] + ai[c] * ai[c]; if (v > T + 0.01) { ok = 0; break; } vv[c] = (float)v; }
        if (!ok) continue;
        for (int c = 0; c < nch; c++) { ps[c] = (uint8_t)(vv[c] / DV); if (vv[c] < F->fmin[c]) F->fmin[c] = vv[c]; }
        for (int c = nch; c < PS; c++) ps[c] = 0;
        int16_t *k = F->key + (size_t)F->n * nk; uint64_t hs = 0;
        for (int a = 0; a < nI; a++) k[a] = (int16_t)(pafS[a] + q * P[a]);
        for (int j = 0; j < nQ; j++) k[nI + j] = (int16_t)Qd[j];
        for (int j = 0; j < nk; j++) hs += rc[j] * (uint64_t)(int64_t)k[j];
        F->hs[F->n] = hs; F->mask[F->n] = mask; F->n++;
    }
}

typedef struct { uint64_t hs; uint32_t i, j; } Ent;
typedef struct { long n; uint64_t *hs; uint8_t *psd; uint32_t *orig; long *bs; } Sorted;
static void sort_by_class(const Fib *F, uint64_t mm, uint64_t M, Sorted *Sx) {
    Sx->n = F->n; Sx->hs = malloc(sizeof(uint64_t) * (F->n + 1)); Sx->psd = aligned_alloc(32, (size_t)PS * (F->n + 1)); Sx->orig = malloc(sizeof(uint32_t) * (F->n + 1));
    Sx->bs = calloc(M + 1, sizeof(long));
    for (long i = 0; i < F->n; i++) Sx->bs[(F->hs[i] & mm) + 1]++;
    for (uint64_t r = 0; r < M; r++) Sx->bs[r + 1] += Sx->bs[r];
    long *fill = malloc(sizeof(long) * M); memcpy(fill, Sx->bs, sizeof(long) * M);
    for (long i = 0; i < F->n; i++) {
        long qq = fill[F->hs[i] & mm]++;
        Sx->hs[qq] = F->hs[i]; Sx->orig[qq] = (uint32_t)i; memcpy(Sx->psd + (size_t)qq * PS, F->psd + (size_t)i * PS, PS);
    }
    free(fill);
}
static void sorted_free(Sorted *Sx) { free(Sx->hs); free(Sx->psd); free(Sx->orig); free(Sx->bs); }

static long long total_compl; static int nsol_print, maxsol_print;
static void emit_solution(int inst, const Fib *Fs[4], long idx[4]) {
    printf("SOL %d", inst);
    for (int b = 0; b < 4; b++) {
        int kd = b < 2 ? 0 : 1; int t[MAXN]; memset(t, 0, sizeof(t)); uint32_t msk = Fs[b]->mask[idx[b]];
        for (int j = 0; j < nvar[b]; j++) { int o = vorb[b][j]; int e = ((msk >> j) & 1) ? -1 : 1; for (int m = 0; m < omm[kd][o]; m++) t[opos[kd][o][m]] = e * osg[kd][o][m]; }
        printf(" |"); for (int z = 0; z < n; z++) printf(" %d", t[z]);
    }
    printf("\n");
}

// achievable tau sums of a block: DP over variables with weight w = sum of sg over the orbit
static void tau_set(int b, char *set /* size 2*MAXN+1 */) {
    int kd = b < 2 ? 0 : 1; static char cur[2 * MAXN * 4 + 1]; int OFF = 2 * MAXN;
    memset(cur, 0, sizeof(cur)); cur[OFF] = 1;
    for (int j = 0; j < nvar[b]; j++) {
        int o = vorb[b][j]; int w = 0; for (int m = 0; m < omm[kd][o]; m++) w += osg[kd][o][m];
        static char nxt[2 * MAXN * 4 + 1]; memset(nxt, 0, sizeof(nxt));
        for (int x = 0; x <= 2 * OFF; x++) if (cur[x]) { if (x + w >= 0 && x + w <= 2 * OFF) nxt[x + w] = 1; if (x - w >= 0 && x - w <= 2 * OFF) nxt[x - w] = 1; }
        memcpy(cur, nxt, sizeof(cur));
    }
    memset(set, 0, 2 * MAXN + 1);
    for (int x = 0; x <= 2 * OFF; x++) if (cur[x]) { int tv = x - OFF; if (tv >= -MAXN && tv <= MAXN) set[tv + MAXN] = 1; }
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: pl2 file [maxsol_print] [tlimit] [seed]\n"); return 1; }
    FILE *f = fopen(argv[1], "r"); if (!f) return 2;
    maxsol_print = argc > 2 ? atoi(argv[2]) : 5;
    double tlimit = argc > 3 ? atof(argv[3]) : 1e30;
    if (argc > 4) rs ^= (uint64_t)atol(argv[4]) * 0xD1B54A32D192ED03ULL;
    if (fscanf(f, "%d %d %d %d", &n, &q, &eps, &T) != 4) return 2;
    if (fscanf(f, "%d", &nI) != 1) return 2; for (int i = 0; i < nI; i++) if (fscanf(f, "%d", &kI[i]) != 1) return 2;
    if (fscanf(f, "%d", &nQ) != 1) return 2; for (int i = 0; i < nQ; i++) if (fscanf(f, "%d", &kQ[i]) != 1) return 2;
    for (int kd = 0; kd < 2; kd++) {
        if (fscanf(f, "%d", &nu[kd]) != 1) return 2;
        for (int o = 0; o < nu[kd]; o++) { if (fscanf(f, "%d", &omm[kd][o]) != 1) return 2; for (int m = 0; m < omm[kd][o]; m++) if (fscanf(f, "%d %d", &opos[kd][o][m], &osg[kd][o][m]) != 2) return 2; }
    }
    int ninst; if (fscanf(f, "%d", &ninst) != 1) return 2;
    nk = nI + nQ; sqq = sqrt((double)q);
    if (n - 1 > PS) { fprintf(stderr, "PS too small\n"); return 4; }
    // channels
    nch = 0;
    for (int psi = 1; psi < n; psi++) for (int z = 0; z < n; z++) { double a = 2 * M_PI * ((psi * z) % n) / n; sre[psi][z] = cos(a); sim_[psi][z] = sin(a); }
    if (eps == -1) {
        for (int psi = 1; psi < n; psi++) { chpsi[nch] = psi; for (int z = 0; z < n; z++) { // gamma = i sqrt(q):  gamma*(c + i s) = -sqrt(q) s + i sqrt(q) c
                cgr[nch][z] = -sqq * sim_[psi][z]; cgi[nch][z] = sqq * sre[psi][z]; } nch++; }
    } else {
        for (int psi = 1; psi <= (n - 1) / 2; psi++) for (int sg = 1; sg >= -1; sg -= 2) { chpsi[nch] = psi; for (int z = 0; z < n; z++) { cgr[nch][z] = sg * sqq * sre[psi][z]; cgi[nch][z] = sg * sqq * sim_[psi][z]; } nch++; }
    }
    DV = (T + 249) / 250; if (DV < 1) DV = 1;
    for (int j = 0; j < 2 * MAXK; j++) rc[j] = rnd() | 1;
    for (int inst = 0; inst < ninst; inst++) {
        for (int b = 0; b < 4; b++) for (int z = 0; z < n; z++) if (fscanf(f, "%d", &S[b][z]) != 1) return 2;
        for (int b = 0; b < 4; b++) { int kd = b < 2 ? 0 : 1; for (int o = 0; o < nu[kd]; o++) if (fscanf(f, "%d", &Dfl[b][o]) != 1) return 2; }
        double t0 = now();
        for (int b = 0; b < 4; b++) {
            int kd = b < 2 ? 0 : 1; nvar[b] = 0;
            for (int o = 0; o < nu[kd]; o++) if (!Dfl[b][o]) { if (nvar[b] >= MAXV) { fprintf(stderr, "too many vars\n"); return 5; } vorb[b][nvar[b]++] = o; }
            for (int z = 0; z < n; z++) for (int j = 0; j < nQ; j++) { int zp = (z + kQ[j]) % n, zm = ((z - kQ[j]) % n + n) % n; wQ[b][z][j] = S[b][zp] + eps * S[b][zm]; }
        }
        // order the variables so that the fixed one (last) is in the block with the most variables (C or D preferred)
        int fixk = 2; for (int b = 0; b < 4; b++) if (nvar[b] > nvar[fixk]) fixk = b;
        if (nvar[fixk] == 0) fixk = -1;
        // row-sum / tau filters
        int sg2 = 0; int sgm[4]; for (int b = 0; b < 4; b++) { int sr = 0; for (int z = 0; z < n; z++) sr += S[b][z]; sgm[b] = sr; sg2 += sr * sr; }
        static char tset[4][2 * MAXN + 1];
        for (int b = 0; b < 4; b++) tau_set(b, tset[b]);
        memset(allowT, 0, sizeof(allowT));
        {
            // left (A,B) sums, right (C,D)
            int rem0 = T - sg2;
            if (rem0 >= 0) {
                for (int ta = -MAXN; ta <= MAXN; ta++) if (tset[0][ta + MAXN])
                for (int tb = -MAXN; tb <= MAXN; tb++) if (tset[1][tb + MAXN]) {
                    int r1 = rem0 - q * (ta * ta + tb * tb); if (r1 < 0) continue;
                    int cr1 = sgm[0] * ta + sgm[1] * tb;
                    for (int tc = -MAXN; tc <= MAXN; tc++) if (tset[2][tc + MAXN]) {
                        int r2 = r1 - q * tc * tc; if (r2 < 0) continue;
                        for (int td = -MAXN; td <= MAXN; td++) if (tset[3][td + MAXN]) {
                            if (q * td * td != r2) continue;
                            if (eps == 1 && cr1 + sgm[2] * tc + sgm[3] * td != 0) continue;
                            allowT[0][ta + MAXN] = allowT[1][tb + MAXN] = allowT[2][tc + MAXN] = allowT[3][td + MAXN] = 1;
                        }
                    }
                }
            }
        }
        int anyok = 1; for (int b = 0; b < 4; b++) { int a = 0; for (int x = 0; x <= 2 * MAXN; x++) a |= allowT[b][x]; if (!a) anyok = 0; }
        Fib F[4]; memset(F, 0, sizeof(F));
        long long compl = 0, matches = 0, surv1 = 0, surv2 = 0; int timed_out = 0; double npairs1 = 0, npairs2 = 0;
        int big = 0; for (int b = 0; b < 4; b++) if (nvar[b] - (b == fixk ? 1 : 0) > 24) big = 1;
        if (!anyok) { printf("INST %d rowfilter\n", inst); fflush(stdout); continue; }
        if (big) { printf("INST %d TOOBIG nv=%d,%d,%d,%d\n", inst, nvar[0], nvar[1], nvar[2], nvar[3]); fflush(stdout); continue; }
        // move the fixed variable to the end of its block's list
        for (int b = 0; b < 4; b++) build_fibre(b, &F[b], b == fixk);
        double tfib = now() - t0;
        int sx = 0, sy = 3, sz = 1, sw = 2;
        npairs1 = (double)F[sx].n * F[sy].n; npairs2 = (double)F[sz].n * F[sw].n;
        int emptyf = 0; for (int b = 0; b < 4; b++) if (F[b].n == 0) emptyf = 1;
        if (!emptyf) {
            uint8_t lim1[PS] __attribute__((aligned(32))), lim2[PS] __attribute__((aligned(32)));
            for (int c = 0; c < PS; c++) {
                if (c < nch) {
                    float l1 = T + 0.02f - F[sz].fmin[c] - F[sw].fmin[c], l2 = T + 0.02f - F[sx].fmin[c] - F[sy].fmin[c];
                    lim1[c] = (uint8_t)(l1 < 0 ? 0 : (l1 / DV > 255 ? 255 : (int)(l1 / DV)));
                    lim2[c] = (uint8_t)(l2 < 0 ? 0 : (l2 / DV > 255 ? 255 : (int)(l2 / DV)));
                } else { lim1[c] = lim2[c] = 255; }
            }
            double bigp = npairs1 > npairs2 ? npairs1 : npairs2;
            int lgM = 0; while ((bigp / (double)(1L << lgM)) > 20000.0 && lgM < 26) lgM++;
            uint64_t M = 1ull << lgM, mm = M - 1;
            Sorted SX, SY, SZ, SW;
            sort_by_class(&F[sx], mm, M, &SX); sort_by_class(&F[sy], mm, M, &SY); sort_by_class(&F[sz], mm, M, &SZ); sort_by_class(&F[sw], mm, M, &SW);
            size_t lcap = 1 << 16; Ent *list = malloc(sizeof(Ent) * lcap);
            size_t tcap = 1 << 17; Ent *tab = malloc(sizeof(Ent) * tcap);
            const int NV = PS / 32;
            __m256i L1[NV], L2[NV]; for (int v = 0; v < NV; v++) { L1[v] = _mm256_load_si256((const __m256i *)(lim1 + 32 * v)); L2[v] = _mm256_load_si256((const __m256i *)(lim2 + 32 * v)); }
            const __m256i zero = _mm256_setzero_si256();
            for (uint64_t r = 0; r < M && !timed_out; r++) {
                size_t nl = 0;
                for (long i = 0; i < SX.n; i++) {
                    uint64_t bk = (r - SX.hs[i]) & mm;
                    long q0 = SY.bs[bk], q1 = SY.bs[bk + 1]; if (q0 == q1) continue;
                    __m256i xv[NV]; for (int v = 0; v < NV; v++) xv[v] = _mm256_load_si256((const __m256i *)(SX.psd + (size_t)i * PS + 32 * v));
                    for (long qq = q0; qq < q1; qq++) {
                        int pass = 1;
                        for (int v = 0; v < NV && pass; v++) {
                            __m256i sm = _mm256_adds_epu8(xv[v], _mm256_load_si256((const __m256i *)(SY.psd + (size_t)qq * PS + 32 * v)));
                            if (_mm256_movemask_epi8(_mm256_cmpeq_epi8(_mm256_subs_epu8(sm, L1[v]), zero)) != -1) pass = 0;
                        }
                        if (!pass) continue;
                        if (nl >= lcap) { lcap *= 2; list = realloc(list, sizeof(Ent) * lcap); }
                        list[nl].hs = SX.hs[i] + SY.hs[qq]; list[nl].i = (uint32_t)i; list[nl].j = (uint32_t)qq; nl++;
                    }
                }
                surv1 += nl;
                if (!nl) continue;
                size_t need = 1; while (need < 2 * nl) need <<= 1;
                if (need > tcap) { tcap = need; tab = realloc(tab, sizeof(Ent) * tcap); }
                size_t cap = need; uint64_t cm = cap - 1;
                for (size_t qq = 0; qq < cap; qq++) tab[qq].i = 0xffffffffu;
                for (size_t qq = 0; qq < nl; qq++) {
                    uint64_t pos = (list[qq].hs * 0x9E3779B97F4A7C15ULL) >> 20 & cm;
                    while (tab[pos].i != 0xffffffffu) pos = (pos + 1) & cm;
                    tab[pos] = list[qq];
                }
                uint64_t rr = (0 - r) & mm;
                for (long i = 0; i < SZ.n; i++) {
                    uint64_t bk = (rr - SZ.hs[i]) & mm;
                    long q0 = SW.bs[bk], q1 = SW.bs[bk + 1]; if (q0 == q1) continue;
                    __m256i xv[NV]; for (int v = 0; v < NV; v++) xv[v] = _mm256_load_si256((const __m256i *)(SZ.psd + (size_t)i * PS + 32 * v));
                    for (long qq = q0; qq < q1; qq++) {
                        int pass = 1;
                        for (int v = 0; v < NV && pass; v++) {
                            __m256i sm = _mm256_adds_epu8(xv[v], _mm256_load_si256((const __m256i *)(SW.psd + (size_t)qq * PS + 32 * v)));
                            if (_mm256_movemask_epi8(_mm256_cmpeq_epi8(_mm256_subs_epu8(sm, L2[v]), zero)) != -1) pass = 0;
                        }
                        if (!pass) continue;
                        surv2++;
                        uint64_t key = 0 - (SZ.hs[i] + SW.hs[qq]);
                        uint64_t pos = (key * 0x9E3779B97F4A7C15ULL) >> 20 & cm;
                        while (tab[pos].i != 0xffffffffu) {
                            if (tab[pos].hs == key) {
                                matches++;
                                long xo = SX.orig[tab[pos].i], yo = SY.orig[tab[pos].j], zo = SZ.orig[i], wo = SW.orig[qq];
                                const int16_t *kx = F[sx].key + xo * nk, *ky = F[sy].key + yo * nk, *kz = F[sz].key + zo * nk, *kw = F[sw].key + wo * nk;
                                int good = 1;
                                for (int j = 0; j < nk; j++) if (kx[j] + ky[j] + kz[j] + kw[j] != 0) { good = 0; break; }
                                if (good) {
                                    compl++;
                                    if (nsol_print < maxsol_print) { nsol_print++; const Fib *Fs[4] = {&F[0], &F[1], &F[2], &F[3]}; long idx[4] = {xo, zo, wo, yo}; emit_solution(inst, Fs, idx); }
                                }
                            }
                            pos = (pos + 1) & cm;
                        }
                    }
                }
                if (now() - t0 > tlimit) timed_out = 1;
            }
            free(list); free(tab); sorted_free(&SX); sorted_free(&SY); sorted_free(&SZ); sorted_free(&SW);
        }
        long long mult = (fixk >= 0) ? 2 : 1;
        total_compl += compl * mult;
        printf("INST %d nv=%d,%d,%d,%d fib=%ld,%ld,%ld,%ld pairs=%.3g,%.3g surv=%lld,%lld matches=%lld completions=%lld%s t=%.3f tfib=%.3f\n",
               inst, nvar[0], nvar[1], nvar[2], nvar[3], F[0].n, F[1].n, F[2].n, F[3].n, npairs1, npairs2, surv1, surv2, matches,
               compl * mult, timed_out ? " TIMEOUT" : "", now() - t0, tfib);
        fflush(stdout);
        for (int b = 0; b < 4; b++) fib_free(&F[b]);
    }
    printf("TOTAL completions %lld\n", total_compl);
    return 0;
}
