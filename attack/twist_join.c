/*
 * Exact two-list join for the twisted Williamson-type structures of attack/wt_twist.py: every pair (a, y) with
 *
 *     VA[a][c] + v_c(y) = T[c]        for every class c,      v_c(y) = sum_{i,j} S_c[i][j] y_i y_j,
 *
 * where VA are the class values of the first block (a file of int16 rows, e.g. written by bnw_enum) and y runs over
 * every +-1 vector on the O orbits of the second block with row sum r (y_0 = +1 when fix0 = 1).
 *
 * The first list goes into an open-addressing table keyed by a random linear 64-bit hash h(v) = sum_c R_c v_c of its
 * rows. The second block is split into high orbits (0 .. O-L-1) and L low orbits. For each high part (in parallel)
 * the low part runs in Gray order; v(y) = base(high) + low(low) + cross, and the cross term is linear in the low signs,
 * so h(T - v(y)) changes by a precomputed constant per flip. Table hits are re-checked exactly. With P passes
 * (argument 4) only the table entries whose hash is = p (mod P) are held in pass p, for memory.
 *
 * input: int32 O, C, rowsum, fix0; int32 sizes[O]; int32 S[C][O][O]; int16 T[C]
 *   cc -O3 -march=native -fopenmp twist_join.c -o twist_join ; ./twist_join in.bin A.vals [K] [P]
 * Prints "SOL a mask" (bit i of mask set = y_i = -1) for up to K solutions (default 10^6) and "TOTAL n".
 */
#include <omp.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static int O, C, R0, FIX0;
static int32_t *SZ, *S;
static int16_t *T, *VA;
static long nA;
static uint64_t RH[1024];

static double now(void) { struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return t.tv_sec + 1e-9 * t.tv_nsec; }

static uint64_t sm(uint64_t *s) {
    uint64_t z = (*s += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL; z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

static uint64_t hrow(const int16_t *v) {
    uint64_t h = 0;
    for (int c = 0; c < C; c++) h += RH[c] * (uint64_t)(int64_t)v[c];
    return h;
}

static void full_values(const int *y, long *v) {
    for (int c = 0; c < C; c++) {
        long s = 0;
        const int32_t *Sc = S + (long)c * O * O;
        for (int i = 0; i < O; i++) for (int j = 0; j < O; j++) s += (long)Sc[i * O + j] * y[i] * y[j];
        v[c] = s;
    }
}

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: twist_join in.bin A.vals [K] [P]\n"); return 1; }
    FILE *f = fopen(argv[1], "rb");
    int32_t hdr[4];
    if (!f || fread(hdr, 4, 4, f) != 4) return 2;
    O = hdr[0]; C = hdr[1]; R0 = hdr[2]; FIX0 = hdr[3];
    if (O > 63 || C > 1024) { fprintf(stderr, "too large\n"); return 2; }
    SZ = malloc(4 * O); S = malloc(4L * C * O * O); T = malloc(2 * C);
    if (fread(SZ, 4, O, f) != (size_t)O || fread(S, 4, (size_t)C * O * O, f) != (size_t)C * O * O ||
        fread(T, 2, C, f) != (size_t)C) return 2;
    fclose(f);
    long K = argc > 3 ? atol(argv[3]) : 1000000;
    int P = argc > 4 ? atoi(argv[4]) : 1;

    FILE *fa = fopen(argv[2], "rb");
    if (!fa) return 3;
    fseek(fa, 0, SEEK_END); nA = ftell(fa) / (2L * C); fseek(fa, 0, SEEK_SET);
    VA = malloc(2L * C * (nA + 1));
    if (fread(VA, 2, (size_t)C * nA, fa) != (size_t)C * nA) return 3;
    fclose(fa);

    uint64_t seed = 0x5DEECE66DULL;
    for (int c = 0; c < C; c++) RH[c] = sm(&seed) | 1;
    uint64_t hT = hrow(T);

    int L = O - 1 < 22 ? O - 1 : 22;                       /* low orbits */
    int Hn = O - L;
    if (FIX0 && Hn < 1) { L = O - 1; Hn = 1; }
    long nlow = 1L << L;
    /* low tables in Gray order: hash of the low-low part of v, row sum of the low part, flipped bit */
    uint64_t *hL = malloc(8 * nlow);
    int32_t *rL = malloc(4 * nlow);
    {
        int y[64]; long v[1024]; long g[1024 * 24];
        for (int i = 0; i < L; i++) y[i] = 1;
        long rs = 0;
        for (int i = 0; i < L; i++) rs += SZ[Hn + i];
        for (int c = 0; c < C; c++) {
            const int32_t *Sc = S + (long)c * O * O;
            long s = 0;
            for (int i = 0; i < L; i++) for (int j = 0; j < L; j++) s += Sc[(Hn + i) * O + Hn + j];
            v[c] = s;
            for (int j = 0; j < L; j++) {
                long t = 0;
                for (int i = 0; i < L; i++) if (i != j) t += Sc[(Hn + i) * O + Hn + j] + Sc[(Hn + j) * O + Hn + i];
                g[c * L + j] = t;
            }
        }
        for (long k = 0; k < nlow; k++) {
            if (k) {
                int b = __builtin_ctzl(k), yo = y[b];
                for (int c = 0; c < C; c++) {
                    const int32_t *Sc = S + (long)c * O * O;
                    v[c] += -2L * yo * g[c * L + b];
                    for (int j = 0; j < L; j++)
                        if (j != b) g[c * L + j] += -2L * yo * (Sc[(Hn + b) * O + Hn + j] + Sc[(Hn + j) * O + Hn + b]);
                }
                y[b] = -yo; rs += -2L * yo * SZ[Hn + b];
            }
            uint64_t h = 0;
            for (int c = 0; c < C; c++) h += RH[c] * (uint64_t)(int64_t)v[c];
            hL[k] = h; rL[k] = (int32_t)rs;
        }
    }
    long lowsum = 0;
    for (int i = 0; i < L; i++) lowsum += SZ[Hn + i];

    long nhigh = 1L << (Hn - FIX0);
    long total = 0, printed = 0;
    fprintf(stderr, "twist_join: nA = %ld, O = %d (high %d, low %d), C = %d, passes %d\n", nA, O, Hn, L, C, P);
    double t0 = now();
    for (int pass = 0; pass < P; pass++) {
        /* table of the first list (entries of this pass) */
        long cnt = 0;
        for (long a = 0; a < nA; a++) if (hrow(VA + a * C) % P == (uint64_t)pass) cnt++;
        int kb = 1;
        while ((1L << kb) < 2 * cnt + 2) kb++;
        long tsz = 1L << kb, tmask = tsz - 1;
        uint64_t *keys = malloc(8 * tsz);
        uint32_t *idx = calloc(tsz, 4);
        for (long a = 0; a < nA; a++) {
            uint64_t h = hrow(VA + a * C);
            if (h % P != (uint64_t)pass) continue;
            long pos = (long)(h >> (64 - kb));
            while (idx[pos]) pos = (pos + 1) & tmask;
            keys[pos] = h; idx[pos] = (uint32_t)(a + 1);
        }
        fprintf(stderr, "pass %d: %ld entries in a table of 2^%d (%.1fs)\n", pass, cnt, kb, now() - t0);
        #pragma omp parallel for schedule(dynamic, 16) reduction(+:total)
        for (long hv = 0; hv < nhigh; hv++) {
            int y[64];
            long rH = 0;
            for (int i = 0; i < Hn; i++) {
                int bit = FIX0 ? (i == 0 ? 0 : (int)((hv >> (i - 1)) & 1)) : (int)((hv >> i) & 1);
                y[i] = bit ? -1 : 1; rH += (long)SZ[i] * y[i];
            }
            long need = R0 - rH;
            if (need > lowsum || need < -lowsum) continue;
            uint64_t hb = 0;
            uint64_t zeta[64];
            for (int b = 0; b < L; b++) zeta[b] = 0;
            for (int c = 0; c < C; c++) {
                const int32_t *Sc = S + (long)c * O * O;
                long s = 0;
                for (int i = 0; i < Hn; i++) for (int j = 0; j < Hn; j++) s += (long)Sc[i * O + j] * y[i] * y[j];
                hb += RH[c] * (uint64_t)(int64_t)s;
                for (int b = 0; b < L; b++) {
                    long z = 0;
                    for (int i = 0; i < Hn; i++) z += (long)(Sc[i * O + Hn + b] + Sc[(Hn + b) * O + i]) * y[i];
                    zeta[b] += RH[c] * (uint64_t)(int64_t)z;
                }
            }
            uint64_t hdot = 0;
            int ys[64];
            for (int b = 0; b < L; b++) { ys[b] = 1; hdot += zeta[b]; }
            for (long k = 0; k < nlow; k++) {
                if (k) {
                    int b = __builtin_ctzl(k);
                    hdot += (uint64_t)(-2LL * ys[b]) * zeta[b];
                    ys[b] = -ys[b];
                }
                if (rL[k] != need) continue;
                uint64_t key = hT - (hb + hL[k] + hdot);
                if (key % P != (uint64_t)pass) continue;
                long pos = (long)(key >> (64 - kb));
                while (idx[pos]) {
                    if (keys[pos] == key) {
                        long a = idx[pos] - 1;
                        int yy[64]; long v[1024];
                        for (int i = 0; i < Hn; i++) yy[i] = y[i];
                        for (int b = 0; b < L; b++) yy[Hn + b] = ys[b];
                        full_values(yy, v);
                        int ok = 1;
                        for (int c = 0; c < C && ok; c++) ok = (VA[a * C + c] + v[c] == T[c]);
                        if (ok) {
                            total++;
                            uint64_t mask = 0;
                            for (int i = 0; i < O; i++) if (yy[i] < 0) mask |= 1ULL << i;
                            #pragma omp critical
                            {
                                if (printed < K) { printf("SOL %ld %llu\n", a, (unsigned long long)mask); fflush(stdout); printed++; }
                            }
                        }
                    }
                    pos = (pos + 1) & tmask;
                }
            }
        }
        free(keys); free(idx);
        fprintf(stderr, "pass %d done: %ld solutions so far (%.1fs)\n", pass, total, now() - t0);
    }
    printf("TOTAL %ld\n", total);
    return 0;
}
