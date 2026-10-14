/*
 * Exact 4-list finder: all (a,b,c,d) in A x B x C x D with VA[a]+VB[b]+VC[c]+VD[d] = T (componentwise).
 * Left pairs (a,b), right pairs (c,d), partitioned exactly on the values of the first two classes,
 * joined on a 64-bit linear hash and re-checked exactly. Prints "SOL a b c d" for each solution
 * (up to K) and the total count.
 *
 * input: int32 nA nB nC nD C; int16 T[C]; int16 VA[nA][C] VB[nB][C] VC[nC][C] VD[nD][C]
 *   cc -O3 -march=native -fopenmp find4.c -o find4 ;  ./find4 input.bin [K]
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

static int C; static int16_t *T, *V[4]; static long nn[4];
static uint64_t R[512], *h[4], hT;
typedef struct { uint64_t h; uint32_t i, j; } P;
static int cmp(const void *x, const void *y) { uint64_t a = ((const P*)x)->h, b = ((const P*)y)->h; return a < b ? -1 : a > b; }

typedef struct { int lo0, lo1, W0, W1; long *start; uint32_t *idx; } Bk;
static Bk bucket(int L, int lo0, int lo1, int W0, int W1) {
    Bk b = {lo0, lo1, W0, W1, calloc((size_t)W0 * W1 + 1, 8), malloc(4 * (nn[L] + 1))};
    for (long i = 0; i < nn[L]; i++) b.start[(V[L][i*C] - lo0) * W1 + (V[L][i*C+1] - lo1) + 1]++;
    for (long k = 0; k < (long)W0 * W1; k++) b.start[k+1] += b.start[k];
    long *pos = malloc(8 * ((size_t)W0 * W1)); memcpy(pos, b.start, 8 * ((size_t)W0 * W1));
    for (long i = 0; i < nn[L]; i++) b.idx[pos[(V[L][i*C] - lo0) * W1 + (V[L][i*C+1] - lo1)]++] = i;
    free(pos); return b;
}

int main(int argc, char **argv) {
    FILE *f = fopen(argv[1], "rb"); int hdr[5];
    if (fread(hdr, 4, 5, f) != 5) return 1;
    for (int k = 0; k < 4; k++) nn[k] = hdr[k];
    C = hdr[4]; long K = argc > 2 ? atol(argv[2]) : 20;
    T = malloc(2 * C); if (fread(T, 2, C, f) != (size_t)C) return 1;
    for (int k = 0; k < 4; k++) { V[k] = malloc(2 * nn[k] * C + 2); if (fread(V[k], 2, nn[k] * C, f) != (size_t)(nn[k] * C)) return 1; }
    uint64_t s = 0x9E3779B97F4A7C15ULL;
    for (int c = 0; c < C; c++) { s ^= s << 13; s ^= s >> 7; s ^= s << 17; R[c] = s | 1; }
    hT = 0; for (int c = 0; c < C; c++) hT += R[c] * (uint64_t)(int64_t)T[c];
    int lo0 = 1 << 30, hi0 = -(1 << 30), lo1 = 1 << 30, hi1 = -(1 << 30);
    for (int k = 0; k < 4; k++) {
        h[k] = malloc(8 * (nn[k] + 1));
        for (long i = 0; i < nn[k]; i++) {
            uint64_t x = 0; for (int c = 0; c < C; c++) x += R[c] * (uint64_t)(int64_t)V[k][i*C+c]; h[k][i] = x;
            int a = V[k][i*C], b = V[k][i*C+1];
            if (a < lo0) lo0 = a; if (a > hi0) hi0 = a; if (b < lo1) lo1 = b; if (b > hi1) hi1 = b;
        }
    }
    int W0 = hi0 - lo0 + 1, W1 = hi1 - lo1 + 1;
    Bk B[4]; for (int k = 0; k < 4; k++) B[k] = bucket(k, lo0, lo1, W0, W1);
    long total = 0, printed = 0;
    int K0 = 2 * W0 - 1, K1 = 2 * W1 - 1;
    #pragma omp parallel for schedule(dynamic, 1) reduction(+:total)
    for (int key = 0; key < K0 * K1; key++) {
        int k0 = key / K1 + 2 * lo0, k1 = key % K1 + 2 * lo1;     /* VA+VB on classes 0,1 */
        int t0 = T[0] - k0, t1 = T[1] - k1;                       /* VC+VD on classes 0,1 */
        long nl = 0, nr = 0;
        for (int a0 = lo0; a0 <= hi0; a0++) for (int a1 = lo1; a1 <= hi1; a1++) {
            int ka = (a0 - lo0) * W1 + (a1 - lo1);
            int b0 = k0 - a0, b1 = k1 - a1;
            if (b0 >= lo0 && b0 <= hi0 && b1 >= lo1 && b1 <= hi1) {
                int kb = (b0 - lo0) * W1 + (b1 - lo1);
                nl += (B[0].start[ka+1] - B[0].start[ka]) * (B[1].start[kb+1] - B[1].start[kb]);
            }
            int c0 = t0 - a0, c1 = t1 - a1;
            if (c0 >= lo0 && c0 <= hi0 && c1 >= lo1 && c1 <= hi1) {
                int kc = (c0 - lo0) * W1 + (c1 - lo1);
                nr += (B[2].start[ka+1] - B[2].start[ka]) * (B[3].start[kc+1] - B[3].start[kc]);
            }
        }
        if (!nl || !nr) continue;
        P *Lp = malloc(sizeof(P) * nl), *Rp = malloc(sizeof(P) * nr); long u = 0, w = 0;
        for (int a0 = lo0; a0 <= hi0; a0++) for (int a1 = lo1; a1 <= hi1; a1++) {
            int ka = (a0 - lo0) * W1 + (a1 - lo1);
            int b0 = k0 - a0, b1 = k1 - a1;
            if (b0 >= lo0 && b0 <= hi0 && b1 >= lo1 && b1 <= hi1) {
                int kb = (b0 - lo0) * W1 + (b1 - lo1);
                for (long x = B[0].start[ka]; x < B[0].start[ka+1]; x++)
                    for (long y = B[1].start[kb]; y < B[1].start[kb+1]; y++)
                        Lp[u++] = (P){h[0][B[0].idx[x]] + h[1][B[1].idx[y]], B[0].idx[x], B[1].idx[y]};
            }
            int c0 = t0 - a0, c1 = t1 - a1;
            if (c0 >= lo0 && c0 <= hi0 && c1 >= lo1 && c1 <= hi1) {
                int kc = (c0 - lo0) * W1 + (c1 - lo1);
                for (long x = B[2].start[ka]; x < B[2].start[ka+1]; x++)
                    for (long y = B[3].start[kc]; y < B[3].start[kc+1]; y++)
                        Rp[w++] = (P){hT - h[2][B[2].idx[x]] - h[3][B[3].idx[y]], B[2].idx[x], B[3].idx[y]};
            }
        }
        qsort(Lp, nl, sizeof(P), cmp);
        long cnt = 0;
        for (long r = 0; r < nr; r++) {
            long a = 0, b = nl;
            while (a < b) { long m = (a + b) / 2; if (Lp[m].h < Rp[r].h) a = m + 1; else b = m; }
            for (; a < nl && Lp[a].h == Rp[r].h; a++) {
                long ia = Lp[a].i, ib = Lp[a].j, ic = Rp[r].i, id = Rp[r].j; int ok = 1;
                for (int c = 0; c < C && ok; c++)
                    ok = V[0][ia*C+c] + V[1][ib*C+c] + V[2][ic*C+c] + V[3][id*C+c] == T[c];
                if (ok) {
                    cnt++;
                    #pragma omp critical
                    { if (printed < K) { printf("SOL %ld %ld %ld %ld\n", ia, ib, ic, id); fflush(stdout); printed++; } }
                }
            }
        }
        total += cnt; free(Lp); free(Rp);
    }
    printf("TOTAL %ld\n", total);
    return 0;
}
