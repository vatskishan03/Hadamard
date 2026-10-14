// Three-level exact join for cyclotomic-product Goethals-Seidel structures (search/cyclo_product.py).
//
// Each of the four blocks is a list of sign vectors grouped by a coarse key (the "level-B" vector, DB int16
// components); every member also carries a fine vector (the "level-C" vector, DC int16 components). A solution
// is one member of each list whose level-B vectors and level-C vectors both sum to zero.
//
// Group pairs of lists 0 and 1 are stored in a hash table keyed by their level-B sum. The program walks every
// group pair of lists 2 and 3, looks up the complementary level-B sums and, for each matching quadruple of
// groups, joins the four groups exactly on the level-C vectors (members of lists 0 x 1 in a small hash table,
// members of lists 2 x 3 probed against it). Every hit is checked component by component.
//
//   join3 input.bin [K] [maxquads]   ->  "SOL i0 i1 i2 i3" (member ids as stored), then "TOTAL found quads work"
//
// Input (native byte order): int32 DB, int32 DC; then per list: int32 nG, int32 0, int64 N, int16 B[nG][DB],
// int64 offs[nG + 1], int16 C[N][DC], int64 id[N] (members of group g are offs[g] .. offs[g+1]-1).
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <omp.h>

typedef struct { int nG; long long N; int16_t *B; long long *offs; int16_t *C; long long *id; } List;
static List L[4];
static int DB, DC;

static void rd(void *p, size_t sz, size_t n, FILE *f) { if (fread(p, sz, n, f) != n) { fprintf(stderr, "join3: short read\n"); exit(1); } }
static inline uint64_t mix(uint64_t x) { x ^= x >> 33; x *= 0xff51afd7ed558ccdULL; x ^= x >> 33; x *= 0xc4ceb9fe1a85ec53ULL; x ^= x >> 33; return x; }
static inline uint64_t hvec(const int16_t *a, const int16_t *b, int d, int neg) {
    uint64_t x = 0x9e3779b97f4a7c15ULL;
    for (int i = 0; i < d; i++) { int v = a[i] + b[i]; if (neg) v = -v; x = mix(x + (uint64_t)(uint16_t)v); }
    return x == ~0ULL ? x - 1 : x;
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: join3 input.bin [K] [maxquads]\n"); return 2; }
    FILE *f = fopen(argv[1], "rb"); if (!f) { perror(argv[1]); return 1; }
    long K = argc > 2 ? atol(argv[2]) : 1; long long maxquads = argc > 3 ? atoll(argv[3]) : -1;
    rd(&DB, 4, 1, f); rd(&DC, 4, 1, f);
    for (int r = 0; r < 4; r++) {
        int nG, pad; long long N; rd(&nG, 4, 1, f); rd(&pad, 4, 1, f); rd(&N, 8, 1, f);
        L[r].nG = nG; L[r].N = N;
        L[r].B = malloc(sizeof(int16_t) * (size_t)DB * nG); rd(L[r].B, 2, (size_t)DB * nG, f);
        L[r].offs = malloc(8 * ((size_t)nG + 1)); rd(L[r].offs, 8, (size_t)nG + 1, f);
        L[r].C = malloc(sizeof(int16_t) * (size_t)DC * N); rd(L[r].C, 2, (size_t)DC * N, f);
        L[r].id = malloc(8 * (size_t)N); rd(L[r].id, 8, (size_t)N, f);
    }
    fclose(f);
    // level-B hash of the group pairs of lists 0 and 1
    long long npair = (long long)L[0].nG * L[1].nG;
    uint64_t size = 1; while (size < 2 * (uint64_t)npair) size <<= 1;
    uint64_t mask = size - 1;
    uint64_t *hk = malloc(8 * size); uint32_t *hv = malloc(4 * size);
    if (!hk || !hv) { fprintf(stderr, "join3: out of memory (%lld group pairs)\n", npair); return 1; }
    memset(hk, 0xff, 8 * size);
    for (int a = 0; a < L[0].nG; a++)
        for (int b = 0; b < L[1].nG; b++) {
            uint64_t x = hvec(L[0].B + (size_t)DB * a, L[1].B + (size_t)DB * b, DB, 0), h = x & mask;
            while (hk[h] != ~0ULL) h = (h + 1) & mask;
            hk[h] = x; hv[h] = (uint32_t)((uint64_t)a * L[1].nG + b);
        }
    fprintf(stderr, "join3: DB=%d DC=%d, lists %lld %lld %lld %lld, level-B group pairs %lld\n", DB, DC,
            L[0].N, L[1].N, L[2].N, L[3].N, npair);
    long found = 0; long long quads = 0, work = 0;
    int n2 = L[2].nG, n3 = L[3].nG;
    #pragma omp parallel reduction(+:quads, work)
    {
        size_t cap = 1 << 16; uint64_t *tk = malloc(8 * cap); long long *t0 = malloc(8 * cap), *t1 = malloc(8 * cap);
        int16_t need[64];
        #pragma omp for schedule(dynamic, 8)
        for (int g2 = 0; g2 < n2; g2++) {
            if (found >= K) continue;
            for (int g3 = 0; g3 < n3; g3++) {
                if (found >= K || (maxquads >= 0 && quads >= maxquads)) break;
                const int16_t *b2 = L[2].B + (size_t)DB * g2, *b3 = L[3].B + (size_t)DB * g3;
                for (int i = 0; i < DB; i++) need[i] = -(b2[i] + b3[i]);
                uint64_t x = 0x9e3779b97f4a7c15ULL;
                for (int i = 0; i < DB; i++) x = mix(x + (uint64_t)(uint16_t)need[i]);
                if (x == ~0ULL) x--;
                for (uint64_t h = x & mask; hk[h] != ~0ULL; h = (h + 1) & mask) {
                    if (hk[h] != x) continue;
                    int g0 = (int)(hv[h] / L[1].nG), g1 = (int)(hv[h] % L[1].nG);
                    const int16_t *b0 = L[0].B + (size_t)DB * g0, *b1 = L[1].B + (size_t)DB * g1;
                    int ok = 1;
                    for (int i = 0; i < DB; i++) if (b0[i] + b1[i] != need[i]) { ok = 0; break; }
                    if (!ok) continue;
                    quads++;
                    long long a0 = L[0].offs[g0], a1 = L[0].offs[g0 + 1], c0 = L[1].offs[g1], c1 = L[1].offs[g1 + 1];
                    long long e0 = L[2].offs[g2], e1 = L[2].offs[g2 + 1], f0 = L[3].offs[g3], f1 = L[3].offs[g3 + 1];
                    long long n01 = (a1 - a0) * (c1 - c0), n23 = (e1 - e0) * (f1 - f0);
                    work += n01 + n23;
                    size_t want = 1; while (want < 2 * (size_t)n01) want <<= 1;
                    if (want > cap) { cap = want; tk = realloc(tk, 8 * cap); t0 = realloc(t0, 8 * cap); t1 = realloc(t1, 8 * cap); }
                    uint64_t tm = want - 1; memset(tk, 0xff, 8 * want);
                    for (long long i = a0; i < a1; i++)
                        for (long long j = c0; j < c1; j++) {
                            uint64_t y = hvec(L[0].C + (size_t)DC * i, L[1].C + (size_t)DC * j, DC, 0), s = y & tm;
                            while (tk[s] != ~0ULL) s = (s + 1) & tm;
                            tk[s] = y; t0[s] = i; t1[s] = j;
                        }
                    for (long long i = e0; i < e1; i++)
                        for (long long j = f0; j < f1; j++) {
                            const int16_t *ci = L[2].C + (size_t)DC * i, *cj = L[3].C + (size_t)DC * j;
                            uint64_t y = hvec(ci, cj, DC, 1);
                            for (uint64_t s = y & tm; tk[s] != ~0ULL; s = (s + 1) & tm) {
                                if (tk[s] != y) continue;
                                const int16_t *ca = L[0].C + (size_t)DC * t0[s], *cb = L[1].C + (size_t)DC * t1[s];
                                int good = 1;
                                for (int d = 0; d < DC; d++) if (ca[d] + cb[d] + ci[d] + cj[d] != 0) { good = 0; break; }
                                if (!good) continue;
                                #pragma omp critical
                                {
                                    if (found < K) { printf("SOL %lld %lld %lld %lld\n", L[0].id[t0[s]], L[1].id[t1[s]], L[2].id[i], L[3].id[j]); fflush(stdout); }
                                    found++;
                                }
                            }
                        }
                }
            }
            if (omp_get_thread_num() == 0 && g2 % 256 == 0)
                fprintf(stderr, "join3: group %d/%d, quadruples %lld, work %.3g, found %ld\n", g2, n2, quads, (double)work, found);
        }
        free(tk); free(t0); free(t1);
    }
    printf("TOTAL %ld quads %lld work %lld\n", found, quads, work);
    return 0;
}
