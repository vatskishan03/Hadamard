/*
 * Four-list k-tree (Wagner) for class-value vectors: find (a, b, c, d) with VA[a] + VB[b] + VC[c] + VD[d] = T on all
 * D classes, among sampled lists.
 *
 * For each split value sigma of the merge classes S (|S| <= 5, values within +-2047), the pairs (a, b) with (VA[a] + VB[b])_S = sigma and
 * the pairs (c, d) with (VC[c] + VD[d])_S = (T - sigma)_S are formed through a sorted index of VB and VD on S, and
 * joined on a 64-bit linear hash of the other classes; hash matches are re-checked on all classes. The values of
 * sigma are taken in order of increasing distance from the centre mu_A + mu_B + (T - mu)/2 (on S), stepping by the
 * lattice of each merge class (the gcd of the differences of its values within each list, combined).
 *
 * input: int32 nA nB nC nD D; int16 T[D]; int16 VA[nA][D] VB[nB][D] VC[nC][D] VD[nD][D]
 *   cc -O3 -march=native -fopenmp kt4.c -o kt4 ;  ./kt4 input.bin S0,S1,.. NSIG [K]
 * Prints "SOL a b c d" (at most K, default 1e9), progress on stderr, and "TOTAL found pairsAB pairsCD".
 */
#include <omp.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

static int D, NS, Sx[5];
static long MAXPAIRS = 300000000L;
static int16_t *T, *V[4];
static long n[4];
static uint64_t Rh[256];

static double now(void) { struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return t.tv_sec + 1e-9 * t.tv_nsec; }
static void rd(void *p, size_t sz, size_t k, FILE *f) { if (fread(p, sz, k, f) != k) { fprintf(stderr, "kt4: short read\n"); exit(1); } }

static inline uint64_t pk(int x) { return (uint64_t)((x + 2048) & 0xFFF); }
static inline uint64_t skey(const int16_t *v, const int *off) {          /* pack the S-values (12 bits each) */
    uint64_t k = 0;
    for (int s = 0; s < NS; s++) k = (k << 12) | pk(v[Sx[s]] + off[s]);
    return k;
}
static inline int isS(int c) { for (int s = 0; s < NS; s++) if (Sx[s] == c) return 1; return 0; }

typedef struct { long m; uint64_t *key; long *start; long *idx; } Index;      /* list sorted by S-key, CSR */

static int cmpk_base;
static uint64_t *cmp_keys;
static int cmpk(const void *a, const void *b) {
    uint64_t x = cmp_keys[*(const long *)a], y = cmp_keys[*(const long *)b];
    return x < y ? -1 : x > y;
}

static Index build_index(int L, const int *off) {
    Index I; long N = n[L];
    uint64_t *k = malloc(8 * N); long *o = malloc(sizeof(long) * N);
    for (long i = 0; i < N; i++) { k[i] = skey(V[L] + (size_t)D * i, off); o[i] = i; }
    cmp_keys = k; qsort(o, N, sizeof(long), cmpk);
    I.idx = o; I.key = malloc(8 * (N + 1)); I.start = malloc(sizeof(long) * (N + 1)); I.m = 0;
    for (long i = 0; i < N; i++) {
        if (i == 0 || k[o[i]] != k[o[i - 1]]) { I.key[I.m] = k[o[i]]; I.start[I.m] = i; I.m++; }
    }
    I.start[I.m] = N; free(k);
    return I;
}

static long find(const Index *I, uint64_t key) {                         /* group index or -1 */
    long lo = 0, hi = I->m - 1;
    while (lo <= hi) { long mid = (lo + hi) >> 1; if (I->key[mid] < key) lo = mid + 1; else if (I->key[mid] > key) hi = mid - 1; else return mid; }
    return -1;
}

int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage: kt4 input.bin S0,S1,.. NSIG [K]\n"); return 2; }
    FILE *f = fopen(argv[1], "rb"); if (!f) { perror(argv[1]); return 1; }
    int hdr[5]; rd(hdr, 4, 5, f); for (int i = 0; i < 4; i++) n[i] = hdr[i]; D = hdr[4];
    T = malloc(2 * D); rd(T, 2, D, f);
    for (int i = 0; i < 4; i++) { V[i] = malloc(2 * (size_t)D * n[i]); rd(V[i], 2, (size_t)D * n[i], f); }
    fclose(f);
    NS = 0; for (char *p = strtok(argv[2], ","); p && NS < 5; p = strtok(NULL, ",")) Sx[NS++] = atoi(p);
    long NSIG = atol(argv[3]); long K = argc > 4 ? atol(argv[4]) : 1000000000L;
    if (argc > 5) MAXPAIRS = atol(argv[5]);
    uint64_t seed = 0x243F6A8885A308D3ULL;
    for (int c = 0; c < D; c++) { seed ^= seed << 13; seed ^= seed >> 7; seed ^= seed << 17; Rh[c] = isS(c) ? 0 : (seed | 1); }
    /* per-list hash over the non-S classes */
    uint64_t *h[4];
    for (int L = 0; L < 4; L++) {
        h[L] = malloc(8 * n[L]);
        #pragma omp parallel for
        for (long i = 0; i < n[L]; i++) { uint64_t x = 0; const int16_t *v = V[L] + (size_t)D * i; for (int c = 0; c < D; c++) x += Rh[c] * (uint64_t)(int64_t)v[c]; h[L][i] = x; }
    }
    uint64_t hT = 0; for (int c = 0; c < D; c++) hT += Rh[c] * (uint64_t)(int64_t)T[c];
    /* means, lattices and offsets for the S classes */
    double mu[4][5]; int g[5], vmin[5], vmax[5];
    for (int s = 0; s < NS; s++) {
        g[s] = 0; vmin[s] = 1 << 30; vmax[s] = -(1 << 30);
        for (int L = 0; L < 4; L++) {
            double acc = 0; int first = V[L][Sx[s]];
            for (long i = 0; i < n[L]; i++) {
                int x = V[L][(size_t)D * i + Sx[s]]; acc += x;
                int d = abs(x - first); while (d) { int t = g[s] % d; g[s] = d; d = t; }
                if (x < vmin[s]) vmin[s] = x; if (x > vmax[s]) vmax[s] = x;
            }
            mu[L][s] = acc / n[L];
        }
        if (g[s] == 0) g[s] = 1;
    }
    int off[5]; for (int s = 0; s < NS; s++) {
        off[s] = 0;
        if (vmin[s] < -1000 || vmax[s] > 1000) { fprintf(stderr, "kt4: class %d out of range\n", Sx[s]); return 1; }
    }
    Index IB = build_index(1, off), ID = build_index(3, off);
    fprintf(stderr, "kt4: lists %ld %ld %ld %ld, D=%d, S=", n[0], n[1], n[2], n[3], D);
    for (int s = 0; s < NS; s++) fprintf(stderr, "%d(g=%d) ", Sx[s], g[s]);
    fprintf(stderr, ", groups B %ld D %ld\n", IB.m, ID.m);
    /* candidate sigmas: lattice points near the centre, by increasing normalised distance */
    double cen[5], sd[5];
    for (int s = 0; s < NS; s++) {
        double tot = 0; for (int L = 0; L < 4; L++) tot += mu[L][s];
        cen[s] = mu[0][s] + mu[1][s] + 0.5 * (T[Sx[s]] - tot);
        sd[s] = 1.0;
    }
    int R = 1; long ncand = 1; while (ncand < 4 * NSIG) { R++; ncand = 1; for (int s = 0; s < NS; s++) ncand *= (2 * R + 1); }
    long nc = 1; for (int s = 0; s < NS; s++) nc *= (2 * R + 1);
    int *cand = malloc(sizeof(int) * nc * NS); double *dist = malloc(sizeof(double) * nc); long *ord = malloc(sizeof(long) * nc);
    /* residue of sigma_s: sum of a value of A and of B is congruent to (vA + vB) mod g */
    int base[5]; for (int s = 0; s < NS; s++) base[s] = V[0][Sx[s]] + V[1][Sx[s]];
    for (long t = 0; t < nc; t++) {
        long q = t; double dd = 0;
        for (int s = 0; s < NS; s++) {
            int k = (int)(q % (2 * R + 1)) - R; q /= (2 * R + 1);
            long c0 = lround((cen[s] - base[s]) / g[s]); int sig = base[s] + (int)(c0 + k) * g[s];
            cand[t * NS + s] = sig; dd += (sig - cen[s]) * (sig - cen[s]);
        }
        dist[t] = dd; ord[t] = t;
    }
    /* sort candidates by distance (simple insertion into an index array by qsort on dist) */
    { double *dk = dist; for (long i = 1; i < nc; i++) { long x = ord[i]; long j = i - 1; while (j >= 0 && dk[ord[j]] > dk[x]) { ord[j + 1] = ord[j]; j--; } ord[j + 1] = x; } }
    long found = 0; long long totAB = 0, totCD = 0; double t0 = now();
    for (long ci = 0; ci < NSIG && ci < nc && found < K; ci++) {
        const int *sig = cand + ord[ci] * NS;
        /* pairs (a, b) with (VA[a] + VB[b])_S = sigma: hash table keyed by hA + hB */
        long cap = 1 << 20, cnt = 0;
        uint64_t *hk = NULL; uint32_t *ha = NULL, *hb = NULL;
        /* count pass */
        long *cntA = calloc(n[0] + 1, sizeof(long));
        #pragma omp parallel for schedule(static)
        for (long a = 0; a < n[0]; a++) {
            const int16_t *va = V[0] + (size_t)D * a;
            uint64_t k = 0; for (int s = 0; s < NS; s++) k = (k << 12) | pk(sig[s] - va[Sx[s]]);
            long gi = find(&IB, k); cntA[a + 1] = gi < 0 ? 0 : IB.start[gi + 1] - IB.start[gi];
        }
        for (long a = 0; a < n[0]; a++) cntA[a + 1] += cntA[a];
        cnt = cntA[n[0]];
        if (cnt > MAXPAIRS) { fprintf(stderr, "kt4: sigma %ld: %ld pairs exceed the cap, skipped\n", ci + 1, cnt); free(cntA); continue; }
        cap = 1; while (cap < 2 * cnt + 2) cap <<= 1;
        hk = malloc(8 * cap); ha = malloc(4 * cap); hb = malloc(4 * cap);
        memset(hk, 0xff, 8 * cap);
        uint64_t mask = cap - 1;
        #pragma omp parallel for schedule(dynamic, 1024)
        for (long a = 0; a < n[0]; a++) {
            if (cntA[a + 1] == cntA[a]) continue;
            const int16_t *va = V[0] + (size_t)D * a; uint64_t k = 0;
            for (int s = 0; s < NS; s++) k = (k << 12) | pk(sig[s] - va[Sx[s]]);
            long gi = find(&IB, k);
            for (long p = IB.start[gi]; p < IB.start[gi + 1]; p++) {
                long b = IB.idx[p]; uint64_t x = h[0][a] + h[1][b]; if (x == ~0ULL) x--;
                uint64_t slot = (x * 0x9E3779B97F4A7C15ULL) & mask;
                while (1) {
                    uint64_t cur = __atomic_load_n(&hk[slot], __ATOMIC_RELAXED);
                    if (cur == ~0ULL && __atomic_compare_exchange_n(&hk[slot], &cur, x, 0, __ATOMIC_RELAXED, __ATOMIC_RELAXED)) {
                        ha[slot] = (uint32_t)a; hb[slot] = (uint32_t)b; break;
                    }
                    slot = (slot + 1) & mask;
                }
            }
        }
        free(cntA);
        /* probe with pairs (c, d) with (VC[c] + VD[d])_S = (T - sigma)_S */
        long long ncd = 0; long fnd = 0;
        #pragma omp parallel for schedule(dynamic, 1024) reduction(+:ncd, fnd)
        for (long c = 0; c < n[2]; c++) {
            const int16_t *vc = V[2] + (size_t)D * c; uint64_t k = 0;
            for (int s = 0; s < NS; s++) k = (k << 12) | pk(T[Sx[s]] - sig[s] - vc[Sx[s]]);
            long gi = find(&ID, k); if (gi < 0) continue;
            for (long p = ID.start[gi]; p < ID.start[gi + 1]; p++) {
                long d = ID.idx[p]; ncd++;
                uint64_t x = hT - h[2][c] - h[3][d]; if (x == ~0ULL) x--;
                uint64_t slot = (x * 0x9E3779B97F4A7C15ULL) & mask;
                while (hk[slot] != ~0ULL) {
                    if (hk[slot] == x) {
                        long a = ha[slot], b = hb[slot]; const int16_t *va = V[0] + (size_t)D * a, *vb = V[1] + (size_t)D * b, *vd = V[3] + (size_t)D * d;
                        int ok = 1; for (int q = 0; q < D; q++) if (va[q] + vb[q] + vc[q] + vd[q] != T[q]) { ok = 0; break; }
                        if (ok) {
                            fnd++;
                            #pragma omp critical
                            { if (found < K) { printf("SOL %ld %ld %ld %ld\n", a, b, c, d); fflush(stdout); } found++; }
                        }
                    }
                    slot = (slot + 1) & mask;
                }
            }
        }
        totAB += cnt; totCD += ncd;
        fprintf(stderr, "kt4: sigma %ld/%ld (", ci + 1, NSIG);
        for (int s = 0; s < NS; s++) fprintf(stderr, "%d%s", sig[s], s + 1 < NS ? "," : "");
        fprintf(stderr, "): pairs AB %ld CD %lld, found %ld (total %ld) [%.0fs]\n", cnt, ncd, fnd, found, now() - t0);
        free(hk); free(ha); free(hb);
    }
    printf("TOTAL %ld %lld %lld\n", found, totAB, totCD);
    return 0;
}
