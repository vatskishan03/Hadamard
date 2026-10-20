/*
 * Exact count of quadruples a in A, b, c, d in B with VA[a] + VB[b] + VB[c] + VB[d] = T
 * (componentwise over C classes). Pairs are partitioned exactly on the values of the
 * first two classes, then joined on a 64-bit linear hash and every hit is re-checked
 * exactly. Optional: print up to K solutions as index quadruples.
 *
 * input file: int32 nA, nB, C; int16 T[C]; int16 VA[nA][C]; int16 VB[nB][C]
 *   cc -O3 -march=native -fopenmp count4.c -o count4
 *   ./count4 input.bin [K]
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

static int nA, nB, C;
static int16_t *T, *VA, *VB;
static uint64_t *hA, *hB, hT, R[256];

typedef struct { uint64_t h; uint32_t i, j; } P;
static int cmp(const void *x, const void *y) {
    uint64_t a = ((const P *)x)->h, b = ((const P *)y)->h; return a < b ? -1 : a > b;
}

int main(int argc, char **argv) {
    FILE *f = fopen(argv[1], "rb");
    int hdr[3]; if (fread(hdr, 4, 3, f) != 3) return 1;
    nA = hdr[0]; nB = hdr[1]; C = hdr[2];
    (void)argc;
    T = malloc(2 * C); VA = malloc(2L * nA * C); VB = malloc(2L * nB * C);
    if (fread(T, 2, C, f) != (size_t)C || fread(VA, 2, (size_t)nA * C, f) != (size_t)nA * C ||
        fread(VB, 2, (size_t)nB * C, f) != (size_t)nB * C) return 1;
    uint64_t s = 0x9E3779B97F4A7C15ULL;
    for (int c = 0; c < C; c++) { s ^= s << 13; s ^= s >> 7; s ^= s << 17; R[c] = s | 1; }
    hA = malloc(8L * nA); hB = malloc(8L * nB); hT = 0;
    for (int c = 0; c < C; c++) hT += R[c] * (uint64_t)(int64_t)T[c];
    for (long i = 0; i < nA; i++) { uint64_t h = 0; for (int c = 0; c < C; c++) h += R[c] * (uint64_t)(int64_t)VA[i * C + c]; hA[i] = h; }
    for (long i = 0; i < nB; i++) { uint64_t h = 0; for (int c = 0; c < C; c++) h += R[c] * (uint64_t)(int64_t)VB[i * C + c]; hB[i] = h; }
    /* buckets on (v0, v1) */
    int lo0 = 1 << 30, hi0 = -(1 << 30), lo1 = 1 << 30, hi1 = -(1 << 30);
    for (long i = 0; i < nA; i++) { int a = VA[i*C], b = VA[i*C+1]; if (a<lo0) lo0=a; if (a>hi0) hi0=a; if (b<lo1) lo1=b; if (b>hi1) hi1=b; }
    for (long i = 0; i < nB; i++) { int a = VB[i*C], b = VB[i*C+1]; if (a<lo0) lo0=a; if (a>hi0) hi0=a; if (b<lo1) lo1=b; if (b>hi1) hi1=b; }
    int W0 = hi0 - lo0 + 1, W1 = hi1 - lo1 + 1, NBk = W0 * W1;
    long *startA = calloc(NBk + 1, 8), *startB = calloc(NBk + 1, 8);
    uint32_t *idxA = malloc(4L * nA), *idxB = malloc(4L * nB);
    for (long i = 0; i < nA; i++) startA[(VA[i*C]-lo0) * W1 + (VA[i*C+1]-lo1) + 1]++;
    for (long i = 0; i < nB; i++) startB[(VB[i*C]-lo0) * W1 + (VB[i*C+1]-lo1) + 1]++;
    for (int k = 0; k < NBk; k++) { startA[k+1] += startA[k]; startB[k+1] += startB[k]; }
    long *posA = malloc(8L * NBk), *posB = malloc(8L * NBk);
    memcpy(posA, startA, 8L * NBk); memcpy(posB, startB, 8L * NBk);
    for (long i = 0; i < nA; i++) idxA[posA[(VA[i*C]-lo0) * W1 + (VA[i*C+1]-lo1)]++] = i;
    for (long i = 0; i < nB; i++) idxB[posB[(VB[i*C]-lo0) * W1 + (VB[i*C+1]-lo1)]++] = i;
    long total = 0;
    /* key = (k0, k1) = value of (VA[a]+VB[b]) on classes 0,1, range [2lo, 2hi] */
    #pragma omp parallel for schedule(dynamic, 1) reduction(+:total)
    for (int key = 0; key < (2*W0-1) * (2*W1-1); key++) {
        int k0 = key / (2*W1-1) + 2*lo0, k1 = key % (2*W1-1) + 2*lo1;
        int t0 = T[0] - k0, t1 = T[1] - k1;           /* VB[c]+VB[d] must equal (t0, t1) */
        long nl = 0, nr = 0;
        for (int a0 = lo0; a0 <= hi0; a0++) for (int a1 = lo1; a1 <= hi1; a1++) {
            int ka = (a0-lo0)*W1 + (a1-lo1);
            int b0 = k0 - a0, b1 = k1 - a1;
            if (!(b0 < lo0 || b0 > hi0 || b1 < lo1 || b1 > hi1)) {
                int kb = (b0-lo0)*W1 + (b1-lo1);
                nl += (startA[ka+1]-startA[ka]) * (startB[kb+1]-startB[kb]);
            }
            int c0 = t0 - a0, c1 = t1 - a1;            /* right pairs: c in B-bucket (a0,a1), d in (c0,c1) */
            if (!(c0 < lo0 || c0 > hi0 || c1 < lo1 || c1 > hi1)) {
                int kc = (c0-lo0)*W1 + (c1-lo1);
                nr += (startB[ka+1]-startB[ka]) * (startB[kc+1]-startB[kc]);
            }
        }
        if (!nl || !nr) continue;   /* (counts and generation now agree) */
        P *Lp = malloc(sizeof(P) * nl); uint64_t *Rp = malloc(8 * nr);
        long u = 0, w = 0;
        for (int a0 = lo0; a0 <= hi0; a0++) for (int a1 = lo1; a1 <= hi1; a1++) {
            int ka = (a0-lo0)*W1 + (a1-lo1);
            int b0 = k0 - a0, b1 = k1 - a1;
            if (!(b0 < lo0 || b0 > hi0 || b1 < lo1 || b1 > hi1)) {
                int kb = (b0-lo0)*W1 + (b1-lo1);
                for (long x = startA[ka]; x < startA[ka+1]; x++) for (long y = startB[kb]; y < startB[kb+1]; y++)
                    Lp[u++] = (P){hA[idxA[x]] + hB[idxB[y]], idxA[x], idxB[y]};
            }
            int c0 = t0 - a0, c1 = t1 - a1;
            if (!(c0 < lo0 || c0 > hi0 || c1 < lo1 || c1 > hi1)) {
                int kc = (c0-lo0)*W1 + (c1-lo1);
                for (long x = startB[ka]; x < startB[ka+1]; x++) for (long y = startB[kc]; y < startB[kc+1]; y++)
                    Rp[w++] = hT - hB[idxB[x]] - hB[idxB[y]];
            }
        }
        qsort(Lp, nl, sizeof(P), cmp);
        long cnt = 0;
        for (long r = 0; r < nr; r++) {
            long a = 0, b = nl;
            while (a < b) { long m = (a + b) / 2; if (Lp[m].h < Rp[r]) a = m + 1; else b = m; }
            while (a < nl && Lp[a].h == Rp[r]) { cnt++; a++; }
        }
        total += cnt;
        free(Lp); free(Rp);
    }
    printf("%ld\n", total);
    return 0;
}
