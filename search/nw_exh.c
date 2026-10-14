/* nw_exh.c: exhaustive search for near-Williamson quadruples (A, B, C, D) on Z_v with a multiplier group H
   (tables from nw_prep.py).

   B, C, D are symmetric and constant on the L-orbits (L = <H, -1>), A is constant on the H-orbits, and
   PAF_A + PAF_B + PAF_C + PAF_D = 0 on every class. Every quadruple satisfies the power-spectrum condition
   |A^(t)|^2 + B^(t)^2 + C^(t)^2 + D^(t)^2 = 4v at every character t, so each partial sum is at most 4v.

   1. Lists of B, C, D with the given row sums, filtered by X^(t)^2 <= 4v.
   2. All triples (B, C, D) with B^2 + C^2 + D^2 <= 4v at every character class (C taken >= B when the row sums
      of B and C agree); the D candidates of a pair are the AND of two bitsets (the D compatible with B, and with
      C). The triples are stored in a hash table keyed by their PAF vector (duplicates kept).
   3. Every A is looked up by -PAF_A. On each L-orbit O u -O, A is either symmetric (equal values on O and -O) or
      antisymmetric (opposite values). For a split E (the L-orbits where A is symmetric) the symmetric part A_s and
      the antisymmetric part A_a are independent, and since A_s^ is real and A_a^ imaginary,
      PAF_A = PAF_{A_s} + PAF_{A_a} and |A^|^2 = A_s^2 + |A_a^|^2 exactly. So per split the members of A are the
      sums of two short lists (Gray-code enumerations with incremental PAF and Fourier values), filtered by the
      power spectrum and probed in blocks with the table slots prefetched.
   So the search is complete for the structure: every quadruple with the row sums is found (A and its reversal,
   which have the same PAF, both appear). Single-precision sums of spectra are compared with 4v + 1 (rounding
   errors stay far below 1); values computed in double precision with 4v + 1e-6.

   usage: nw_exh TABLES CAP_LOG2 [shard nshards]
   Build: cc -O3 -march=native -fopenmp nw_exh.c -o nw_exh -lm
*/
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <omp.h>

#define MAXO 64
#define MAXC 64

static int v, C, O, rA, rB, rC, rD, hsize;
static int sizesL[MAXO];
static int *QL, *QA; static double *GL, *GA; static int *orbLg, *plusg;
static double LIM;     /* 4v + 1e-6, for single values computed in double precision */
static float LIMS;     /* 4v + 1, for sums of single-precision spectra (rounding errors stay below 1e-3) */
#define QLx(c, i, j) QL[((c) * O + (i)) * O + (j)]
#define QAx(c, i, j) QA[((c) * (O - 1) + (i)) * (O - 1) + (j)]
#define GLx(c, i) GL[(c) * O + (i)]
#define GAx(c, i) GA[(c) * (O - 1) + (i)]

static double now(void) { struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return t.tv_sec + 1e-9 * t.tv_nsec; }
static void readall(FILE *f, void *p, size_t n) { if (fread(p, 1, n, f) != n) { fprintf(stderr, "short read\n"); exit(1); } }

typedef struct { int16_t paf[MAXC]; float psd[MAXC]; uint32_t mask; int y0; } Sym;

static int sym_list(int r, Sym **out) {
    int n = O - 1, cap = 1 << 16, cnt = 0;
    Sym *L = malloc(sizeof(Sym) * cap);
    int y[MAXO];
    for (int y0 = -1; y0 <= 1; y0 += 2)
        for (uint32_t m = 0; m < (1u << n); m++) {
            int s = y0;
            for (int i = 0; i < n; i++) s += ((m >> i) & 1 ? -1 : 1) * sizesL[i + 1];
            if (s != r) continue;
            y[0] = y0;
            for (int i = 0; i < n; i++) y[i + 1] = (m >> i) & 1 ? -1 : 1;
            Sym e; e.mask = m; e.y0 = y0;
            int ok = 1;
            for (int c = 0; c < C && ok; c++) {
                double F = 0; for (int i = 0; i < O; i++) F += y[i] * GLx(c, i);
                e.psd[c] = (float)(F * F);
                if (F * F > LIM) ok = 0;
            }
            if (!ok) continue;
            for (int c = 0; c < C; c++) {
                long p = 0;
                for (int i = 0; i < O; i++) { long t = 0; for (int j = 0; j < O; j++) t += (long)QLx(c, i, j) * y[j]; p += t * y[i]; }
                e.paf[c] = (int16_t)p;
            }
            if (cnt == cap) { cap *= 2; L = realloc(L, sizeof(Sym) * cap); }
            L[cnt++] = e;
        }
    *out = L;
    return cnt;
}

static uint64_t K64[MAXC];
static inline uint64_t fp(const int16_t *p) {
    uint64_t h = 0x9E3779B97F4A7C15ull;
    for (int c = 0; c < C; c++) h += (uint64_t)(uint16_t)(p[c] + 1024) * K64[c];
    h ^= h >> 31; h *= 0xD6E8FEB86659FD93ull; h ^= h >> 32;
    return h ? h : 1;
}
static uint64_t *TK, *TV; static uint64_t TMASK;
static void tins(uint64_t k, uint64_t val) {            /* duplicates allowed */
    uint64_t i = k & TMASK;
    for (;;) {
        uint64_t z = 0;
        if (__atomic_compare_exchange_n(&TK[i], &z, k, 0, __ATOMIC_RELAXED, __ATOMIC_RELAXED)) { TV[i] = val; return; }
        i = (i + 1) & TMASK;
    }
}

typedef struct { int16_t paf[MAXC]; float psd[MAXC]; uint32_t bits; } Item;

static long enum_form(int n, const int *Mf, const int *lin0, const int *base, const double *G, const double *F0,
                      const int *wsz, int wtarget, Item **out, long *cap) {
    int x[MAXO], lin[MAXC * MAXO], paf[MAXC]; double F[MAXC];
    for (int i = 0; i < n; i++) x[i] = 1;
    int ws = 0; if (wsz) for (int i = 0; i < n; i++) ws += wsz[i];
    for (int c = 0; c < C; c++) {
        int p = base[c];
        for (int i = 0; i < n; i++) { p += lin0[c * n + i]; for (int j = i + 1; j < n; j++) p += Mf[(c * n + i) * n + j]; }
        paf[c] = p;
        double f = F0[c]; for (int i = 0; i < n; i++) f += G[c * n + i]; F[c] = f;
        for (int i = 0; i < n; i++) { int t = lin0[c * n + i]; for (int j = 0; j < n; j++) if (j != i) t += Mf[(c * n + i) * n + j]; lin[c * MAXO + i] = t; }
    }
    long cnt = 0; uint32_t gray = 0;
    for (uint64_t step = 0;; step++) {
        if (!wsz || ws == wtarget) {
            int ok = 1; Item it;
            for (int c = 0; c < C; c++) { double p = F[c] * F[c]; if (p > LIM) { ok = 0; break; } it.psd[c] = (float)p; }
            if (ok) {
                for (int c = 0; c < C; c++) it.paf[c] = (int16_t)paf[c];
                it.bits = gray;
                if (cnt == *cap) { *cap = *cap * 2 + 1024; *out = realloc(*out, sizeof(Item) * (*cap)); }
                (*out)[cnt++] = it;
            }
        }
        if (n == 0 || step + 1 == (1ull << n)) break;
        int j = __builtin_ctzll(step + 1);
        int xj = x[j];
        for (int c = 0; c < C; c++) {
            paf[c] -= 2 * xj * lin[c * MAXO + j];
            F[c] -= 2 * xj * G[c * n + j];
            const int *row = &Mf[(c * n + j) * n];
            int *lc = &lin[c * MAXO];
            for (int i = 0; i < n; i++) lc[i] -= 2 * xj * row[i];
        }
        x[j] = -xj;
        if (wsz) ws -= 2 * xj * wsz[j];
        gray ^= 1u << j;
    }
    return cnt;
}

static void paf_full(const int *x, long *p) { for (int s = 0; s < v; s++) { long t = 0; for (int g = 0; g < v; g++) t += x[g] * x[(g + s) % v]; p[s] = t; } }
static void sym_seq(const Sym *e, int *y) { y[0] = e->y0; for (int g = 1; g < v; g++) { int k = orbLg[g]; y[g] = ((e->mask >> (k - 1)) & 1) ? -1 : 1; } }

static Sym *LB, *LC, *LD; static int nB, nC, nD, sameBC;

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: nw_exh TABLES CAP_LOG2 [shard nshards]\n"); return 1; }
    FILE *f = fopen(argv[1], "rb"); if (!f) { perror("tables"); return 1; }
    int hdr[8]; readall(f, hdr, sizeof hdr);
    v = hdr[0]; C = hdr[1]; O = hdr[2]; rA = hdr[3]; rB = hdr[4]; rC = hdr[5]; rD = hdr[6]; hsize = hdr[7];
    readall(f, sizesL, sizeof(int) * O);
    QL = malloc(sizeof(int) * C * O * O); readall(f, QL, sizeof(int) * C * O * O);
    QA = malloc(sizeof(int) * C * (O - 1) * (O - 1)); readall(f, QA, sizeof(int) * C * (O - 1) * (O - 1));
    GL = malloc(sizeof(double) * C * O); readall(f, GL, sizeof(double) * C * O);
    GA = malloc(sizeof(double) * C * (O - 1)); readall(f, GA, sizeof(double) * C * (O - 1));
    orbLg = malloc(sizeof(int) * v); plusg = malloc(sizeof(int) * v);
    readall(f, orbLg, sizeof(int) * v); readall(f, plusg, sizeof(int) * v);
    fclose(f);
    int caplog = atoi(argv[2]);
    int shard = argc > 4 ? atoi(argv[3]) : 0, nshards = argc > 4 ? atoi(argv[4]) : 1;
    LIM = 4.0 * v + 1e-6;
    LIMS = (float)(4.0 * v + 1.0);
    srand(12345);
    for (int c = 0; c < C; c++) K64[c] = (((uint64_t)rand() << 33) ^ ((uint64_t)rand() << 11) ^ (uint64_t)rand()) | 1;
    fprintf(stderr, "v=%d C=%d O=%d rows %d|%d,%d,%d |H|=%d threads=%d shard %d/%d\n", v, C, O, rA, rB, rC, rD, hsize, omp_get_max_threads(), shard, nshards);
    double t0 = now();
    nB = sym_list(rB, &LB);
    sameBC = (rB == rC);
    if (sameBC) { LC = LB; nC = nB; } else nC = sym_list(rC, &LC);
    nD = sym_list(rD, &LD);
    fprintf(stderr, "lists after PSD filter: B %d, C %d, D %d (%.1f s)\n", nB, nC, nD, now() - t0);
    uint64_t cap = 1ull << caplog; TMASK = cap - 1;
    TK = calloc(cap, sizeof(uint64_t)); TV = malloc(sizeof(uint64_t) * cap);
    if (!TK || !TV) { fprintf(stderr, "no memory\n"); return 1; }
    /* ---- triples ----
       For each B, the set of D with psd_B + psd_D <= 4v at every class is a bitset; likewise for C. A triple needs
       D in both sets (and the full condition), so the candidates for a pair are the AND of two bitsets. */
    long ntrip = 0, npair = 0, nchk = 0;
    t0 = now();
    int nw64 = (nD + 63) / 64;
    uint64_t *bitsB = calloc((size_t)nB * nw64, 8), *bitsC = sameBC ? bitsB : calloc((size_t)nC * nw64, 8);
    if (!bitsB || !bitsC) { fprintf(stderr, "no memory for bitsets\n"); return 1; }
    #pragma omp parallel for schedule(dynamic, 64)
    for (int i = 0; i < nB; i++)
        for (int k = 0; k < nD; k++) {
            int ok = 1; for (int c = 0; c < C; c++) if (LB[i].psd[c] + LD[k].psd[c] > LIMS) { ok = 0; break; }
            if (ok) bitsB[(size_t)i * nw64 + k / 64] |= 1ull << (k % 64);
        }
    if (!sameBC) {
        #pragma omp parallel for schedule(dynamic, 64)
        for (int j = 0; j < nC; j++)
            for (int k = 0; k < nD; k++) {
                int ok = 1; for (int c = 0; c < C; c++) if (LC[j].psd[c] + LD[k].psd[c] > LIMS) { ok = 0; break; }
                if (ok) bitsC[(size_t)j * nw64 + k / 64] |= 1ull << (k % 64);
            }
    }
    fprintf(stderr, "pair bitsets built (%.1f s)\n", now() - t0); fflush(stderr);
    #pragma omp parallel for schedule(dynamic, 8) reduction(+:ntrip, npair, nchk)
    for (int i = 0; i < nB; i++) {
        const Sym *b = &LB[i];
        const uint64_t *bb = &bitsB[(size_t)i * nw64];
        for (int j = sameBC ? i : 0; j < nC; j++) {
            const Sym *c2 = &LC[j];
            float ps[MAXC]; int ok = 1;
            for (int c = 0; c < C; c++) { ps[c] = b->psd[c] + c2->psd[c]; if (ps[c] > LIMS) { ok = 0; break; } }
            if (!ok) continue;
            npair++;
            const uint64_t *cc = &bitsC[(size_t)j * nw64];
            for (int wd = 0; wd < nw64; wd++) {
                uint64_t m = bb[wd] & cc[wd];
                while (m) {
                    int k = wd * 64 + __builtin_ctzll(m); m &= m - 1;
                    const Sym *d = &LD[k];
                    nchk++;
                    int ok2 = 1;
                    for (int c = 0; c < C; c++) if (ps[c] + d->psd[c] > LIMS) { ok2 = 0; break; }
                    if (!ok2) continue;
                    int16_t s[MAXC];
                    for (int c = 0; c < C; c++) s[c] = b->paf[c] + c2->paf[c] + d->paf[c];
                    tins(fp(s), ((uint64_t)i << 40) | ((uint64_t)j << 20) | (uint64_t)k);
                    ntrip++;
                }
            }
        }
    }
    fprintf(stderr, "triples: %ld (pairs %ld, D checks %ld) in %.1f s; table load %.3f\n", ntrip, npair, nchk, now() - t0, (double)ntrip / cap);
    if ((double)ntrip > 0.8 * cap) { fprintf(stderr, "table too small\n"); return 1; }
    fflush(stderr);
    /* ---- stream A ---- */
    int n = O - 1;
    long nE = 1L << n, hits = 0, probes = 0, ncomb = 0;
    t0 = now();
    double tlast = t0;
    long Edone = 0;
    #pragma omp parallel reduction(+:probes, ncomb)
    {
        Item *X = NULL, *Y = NULL; long capX = 0, capY = 0;
        int *MfS = malloc(sizeof(int) * C * n * n), *MfA = malloc(sizeof(int) * C * n * n);
        int *linS = malloc(sizeof(int) * C * n), *linA = malloc(sizeof(int) * C * n);
        int baseS[MAXC], baseA[MAXC]; double GS_[MAXC * MAXO], GA_[MAXC * MAXO], F0S[MAXC], F0A[MAXC];
        int wsz[MAXO];
        #pragma omp for schedule(dynamic, 32)
        for (long Em = 1; Em < nE; Em++) {
            if ((Em % nshards) != shard) continue;
            int Es[MAXO], Ea[MAXO], ns = 0, na = 0;
            for (int i = 0; i < n; i++) { if ((Em >> i) & 1) Es[ns++] = i + 1; else Ea[na++] = i + 1; }
            for (int a0 = -1; a0 <= 1; a0 += 2) {
                int need = rA - a0;
                int sz = 0; for (int i = 0; i < ns; i++) sz += sizesL[Es[i]];
                if (need > sz || need < -sz || ((sz - need) % 2)) continue;
                for (int c = 0; c < C; c++) {
                    baseS[c] = QLx(c, 0, 0);
                    for (int i = 0; i < ns; i++) baseS[c] += QLx(c, Es[i], Es[i]);
                    F0S[c] = a0 * GLx(c, 0);
                    for (int i = 0; i < ns; i++) {
                        linS[c * ns + i] = a0 * (QLx(c, 0, Es[i]) + QLx(c, Es[i], 0));
                        GS_[c * ns + i] = GLx(c, Es[i]);
                        for (int j = 0; j < ns; j++) MfS[(c * ns + i) * ns + j] = (i == j) ? 0 : QLx(c, Es[i], Es[j]) + QLx(c, Es[j], Es[i]);
                    }
                }
                for (int i = 0; i < ns; i++) wsz[i] = sizesL[Es[i]];
                long nx = enum_form(ns, MfS, linS, baseS, GS_, F0S, wsz, need, &X, &capX);
                if (!nx) continue;
                for (int c = 0; c < C; c++) {
                    baseA[c] = 0;
                    for (int i = 0; i < na; i++) baseA[c] += QAx(c, Ea[i] - 1, Ea[i] - 1);
                    F0A[c] = 0;
                    for (int i = 0; i < na; i++) {
                        linA[c * na + i] = 0;
                        GA_[c * na + i] = GAx(c, Ea[i] - 1);
                        for (int j = 0; j < na; j++) MfA[(c * na + i) * na + j] = (i == j) ? 0 : QAx(c, Ea[i] - 1, Ea[j] - 1) + QAx(c, Ea[j] - 1, Ea[i] - 1);
                    }
                }
                long ny = enum_form(na, MfA, linA, baseA, GA_, F0A, NULL, 0, &Y, &capY);
                if (!ny) continue;
                /* probes in blocks of PB with the table slots prefetched */
                enum { PB = 16 };
                uint64_t pk[PB]; long px[PB], py[PB]; int np_ = 0;
                for (long xi = 0; xi < nx; xi++) {
                    const Item *a = &X[xi];
                    for (long yi = 0; yi <= ny; yi++) {
                        if (yi < ny) {
                            const Item *b = &Y[yi];
                            ncomb++;
                            int ok = 1;
                            for (int c = 0; c < C; c++) if (a->psd[c] + b->psd[c] > LIMS) { ok = 0; break; }
                            if (ok) {
                                int16_t tg[MAXC];
                                for (int c = 0; c < C; c++) tg[c] = -(a->paf[c] + b->paf[c]);
                                uint64_t k = fp(tg);
                                __builtin_prefetch(&TK[k & TMASK]);
                                pk[np_] = k; px[np_] = xi; py[np_] = yi; np_++;
                            }
                        }
                        if (np_ == PB || (yi == ny && np_ > 0 && xi == nx - 1) || (yi == ny && np_ > 0)) {
                            for (int t = 0; t < np_; t++) {
                                uint64_t k = pk[t];
                                probes++;
                                uint64_t h = k & TMASK;
                                while (TK[h]) {
                                    if (TK[h] == k) {
                                        const Item *a2 = &X[px[t]], *b2 = &Y[py[t]];
                                        uint64_t val = TV[h];
                                        int ib = (int)(val >> 40), ic = (int)((val >> 20) & 0xFFFFF), id = (int)(val & 0xFFFFF);
                                        int A[2048], B_[2048], C_[2048], D_[2048]; long pA[2048], pB[2048], pC[2048], pD[2048];
                                        A[0] = a0;
                                        for (int g = 1; g < v; g++) {
                                            int kk = orbLg[g]; int pos;
                                            int isS = 0; for (pos = 0; pos < ns; pos++) if (Es[pos] == kk) { isS = 1; break; }
                                            if (isS) A[g] = ((a2->bits >> pos) & 1) ? -1 : 1;
                                            else { for (pos = 0; pos < na; pos++) if (Ea[pos] == kk) break; int xv = ((b2->bits >> pos) & 1) ? -1 : 1; A[g] = plusg[g] > 0 ? xv : -xv; }
                                        }
                                        sym_seq(&LB[ib], B_); sym_seq(&LC[ic], C_); sym_seq(&LD[id], D_);
                                        paf_full(A, pA); paf_full(B_, pB); paf_full(C_, pC); paf_full(D_, pD);
                                        int good = pA[0] + pB[0] + pC[0] + pD[0] == 4 * v;
                                        for (int s2 = 1; s2 < v && good; s2++) if (pA[s2] + pB[s2] + pC[s2] + pD[s2] != 0) good = 0;
                                        #pragma omp critical
                                        {
                                            if (good) {
                                                hits++;
                                                printf("HIT v=%d A=", v); for (int g = 0; g < v; g++) putchar(A[g] > 0 ? '+' : '-');
                                                printf(" B="); for (int g = 0; g < v; g++) putchar(B_[g] > 0 ? '+' : '-');
                                                printf(" C="); for (int g = 0; g < v; g++) putchar(C_[g] > 0 ? '+' : '-');
                                                printf(" D="); for (int g = 0; g < v; g++) putchar(D_[g] > 0 ? '+' : '-');
                                                printf("\n"); fflush(stdout);
                                            }
                                        }
                                    }
                                    h = (h + 1) & TMASK;
                                }
                            }
                            np_ = 0;
                        }
                    }
                }
            }
            #pragma omp atomic
            Edone++;
            if (omp_get_thread_num() == 0 && now() - tlast > 60) {
                tlast = now();
                fprintf(stderr, "progress: splits %ld of %ld (%.1f%%), hits %ld, %.0f s\n", Edone, (nE - 1) / nshards, 100.0 * Edone * nshards / (nE - 1), hits, tlast - t0);
                fflush(stderr);
            }
        }
        free(X); free(Y); free(MfS); free(MfA); free(linS); free(linA);
    }
    fprintf(stderr, "DONE combinations %ld, probes %ld, hits %ld, %.1f s\n", ncomb, probes, hits, now() - t0);
    return 0;
}
