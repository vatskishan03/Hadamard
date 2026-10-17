/* Exhaustive search for skew-propus quadruples of odd length n (search/NOTES.md section 31):
 *   A3 symmetric, A1 (used twice), A0 of skew type (A0(0) = 1, A0(-x) = -A0(x)), all +-1 on Z_n, with
 *   PAF_A3(s) + 2 PAF_A1(s) + PAF_A0(s) = 0 for s = 1 .. (n - 1)/2.
 * Row sums a3, a1 with a3^2 + 2 a1^2 = 4n - 1 (the row sum of A0 is 1); negating A3 or A1 keeps a solution, so both
 * are taken positive.
 *
 * Meet in the middle. Side 1: A3 up to the multipliers (Z_n^* acts on the solutions and fixes the symmetry types), and
 * every A0; a pair survives the power-spectrum test |A3^(k)|^2 + |A0^(k)|^2 <= 4n and is stored under the
 * autocorrelation vector -(PAF_A3 + PAF_A0)/2 that A1 must have (exact 7-bit fields, no hashing collisions are trusted).
 * Side 2: every A1 with the row sum a1 and A1(0) = -1 (every class under shifts has such a member), its
 * autocorrelations from popcounts, looked up in the table. Each hit is checked from scratch and printed as
 *     SOL <A3> <A1> <A0>     (+/- strings)
 * The search is exhaustive: no solution exists for (n, a3, a1) if none is printed.
 *
 *   cc -O3 -march=native -fopenmp -o sp_mitm search/sp_mitm.c -lm
 *   ./sp_mitm n a3 a1 [threads]
 */
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

static int n, h;
static uint64_t full;

static inline uint64_t rot(uint64_t m, int s) { return ((m << s) | (m >> (n - s))) & full; }

/* autocorrelation of the +-1 sequence with -1 exactly on the bits of m */
static inline void paf_mask(uint64_t m, int *p) {
    for (int s = 1; s <= h; s++) p[s] = n - 2 * __builtin_popcountll(m ^ rot(m, s));
}

typedef struct { uint64_t k0, k1; uint32_t a3, a0; } entry;
static entry *tab;
static uint64_t cap, used;
static uint64_t *bloom, bmask;         /* two-probe Bloom filter in front of the table */

static inline uint64_t mix(uint64_t a, uint64_t b) {
    uint64_t x = a * 0x9E3779B97F4A7C15ULL ^ (b + 0x632BE59BD9B4E019ULL + (a << 6) + (a >> 2));
    x ^= x >> 31; x *= 0xBF58476D1CE4E5B9ULL; x ^= x >> 29;
    return x;
}

/* pack PAF values p[1..h] (each in [-n, n], n < 64) into 7-bit fields */
static inline void pack(const int *p, uint64_t *k0, uint64_t *k1) {
    uint64_t a = 0, b = 0;
    for (int s = 1; s <= h; s++) {
        uint64_t f = (uint64_t)(p[s] + 64) & 127;
        if (s <= 9) a |= f << (7 * (s - 1)); else b |= f << (7 * (s - 10));
    }
    *k0 = a; *k1 = b;
}

static void insert(uint64_t k0, uint64_t k1, uint32_t a3, uint32_t a0) {
    uint64_t i = mix(k0, k1) & (cap - 1);
    while (tab[i].k0 | tab[i].k1) i = (i + 1) & (cap - 1);
    tab[i].k0 = k0; tab[i].k1 = k1; tab[i].a3 = a3; tab[i].a0 = a0;
    used++;
    uint64_t h = mix(k1, k0), b1 = h & bmask, b2 = (h >> 32) & bmask;
    bloom[b1 >> 6] |= 1ULL << (b1 & 63); bloom[b2 >> 6] |= 1ULL << (b2 & 63);
}

static inline int maybe(uint64_t k0, uint64_t k1) {
    uint64_t h = mix(k1, k0), b1 = h & bmask, b2 = (h >> 32) & bmask;
    return (bloom[b1 >> 6] >> (b1 & 63) & 1) && (bloom[b2 >> 6] >> (b2 & 63) & 1);
}

static int seq_of_mask(uint64_t m, int *x) { int r = 0; for (int i = 0; i < n; i++) { x[i] = (m >> i & 1) ? -1 : 1; r += x[i]; } return r; }

static int full_check(uint64_t m3, uint64_t m1, uint64_t m0) {
    int x3[64], x1[64], x0[64];
    seq_of_mask(m3, x3); seq_of_mask(m1, x1); seq_of_mask(m0, x0);
    for (int i = 1; i < n; i++) if (x3[i] != x3[n - i] || x0[i] != -x0[n - i]) return 0;
    if (x0[0] != 1) return 0;
    for (int s = 1; s < n; s++) {
        long t = 0;
        for (int i = 0; i < n; i++) t += x3[i] * x3[(i + s) % n] + 2 * x1[i] * x1[(i + s) % n] + x0[i] * x0[(i + s) % n];
        if (t) return 0;
    }
    return 1;
}

static void print_seq(uint64_t m) { for (int i = 0; i < n; i++) putchar((m >> i & 1) ? '-' : '+'); }

int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage: %s n a3 a1 [threads]\n", argv[0]); return 1; }
    n = atoi(argv[1]); int a3 = atoi(argv[2]), a1 = atoi(argv[3]);
    if (argc > 4) omp_set_num_threads(atoi(argv[4]));
    if (n % 2 == 0 || n > 63 || a3 * a3 + 2 * a1 * a1 != 4 * n - 1) { fprintf(stderr, "bad parameters\n"); return 1; }
    h = (n - 1) / 2;
    full = (n == 64) ? ~0ULL : ((1ULL << n) - 1);
    int units[64], nu = 0;
    for (int u = 1; u < n; u++) { int a = u, b = n; while (b) { int t = a % b; a = b; b = t; } if (a == 1) units[nu++] = u; }

    /* side 1a: symmetric A3 with row sum a3, canonical under multipliers (minimal mask) */
    int nh = h;                          /* pairs {i, n - i}, i = 1..h */
    uint64_t *r3 = malloc(sizeof(uint64_t) << (nh + 1)); int n3 = 0;
    for (uint64_t pat = 0; pat < (1ULL << (nh + 1)); pat++) {
        uint64_t m = (pat & 1) ? 1 : 0;   /* bit 0: A3(0) = -1 */
        for (int i = 1; i <= nh; i++) if (pat >> i & 1) m |= (1ULL << i) | (1ULL << (n - i));
        if (n - 2 * __builtin_popcountll(m) != a3) continue;
        int canon = 1;
        for (int k = 1; k < nu && canon; k++) {
            uint64_t mm = 0;
            for (int i = 0; i < n; i++) if (m >> i & 1) mm |= 1ULL << ((long)i * units[k] % n);
            if (mm < m) canon = 0;
        }
        if (canon) r3[n3++] = m;
    }
    /* side 1b: skew-type A0 */
    int n0 = 1 << nh;
    uint64_t *r0 = malloc(sizeof(uint64_t) * n0);
    for (int pat = 0; pat < n0; pat++) {
        uint64_t m = 0;
        for (int i = 1; i <= nh; i++) m |= (pat >> (i - 1) & 1) ? (1ULL << i) : (1ULL << (n - i));
        r0[pat] = m;
    }
    /* autocorrelations and power spectra */
    int *p3 = malloc(sizeof(int) * n3 * (h + 1)), *p0 = malloc(sizeof(int) * n0 * (h + 1));
    float *f3 = malloc(sizeof(float) * n3 * (h + 1)), *f0 = malloc(sizeof(float) * n0 * (h + 1));
    double cs[64][64];
    for (int k = 0; k < n; k++) for (int s = 0; s < n; s++) cs[k][s] = cos(2 * M_PI * k * s / n);
    for (int which = 0; which < 2; which++) {
        int cnt = which ? n0 : n3; uint64_t *r = which ? r0 : r3; int *p = which ? p0 : p3; float *f = which ? f0 : f3;
        for (int j = 0; j < cnt; j++) {
            paf_mask(r[j], p + j * (h + 1));
            for (int k = 1; k <= h; k++) {        /* |X^(k)|^2 = n + 2 sum_{s=1}^{h} PAF(s) cos(2 pi k s / n) */
                double v = n;
                for (int s = 1; s <= h; s++) v += 2 * p[j * (h + 1) + s] * cs[k][s];
                f[j * (h + 1) + k] = (float)v;
            }
        }
    }
    fprintf(stderr, "n = %d, rows (%d, %d, 1): %d classes of A3, %d A0\n", n, a3, a1, n3, n0);

    /* count survivors, then build the table */
    float lim = 4.0f * n + 0.01f;
    uint64_t surv = 0;
    #pragma omp parallel for reduction(+:surv) schedule(dynamic)
    for (int i = 0; i < n3; i++)
        for (int j = 0; j < n0; j++) {
            int ok = 1;
            for (int k = 1; k <= h && ok; k++) if (f3[i * (h + 1) + k] + f0[j * (h + 1) + k] > lim) ok = 0;
            if (!ok) continue;
            for (int s = 1; s <= h && ok; s++) { int t = p3[i * (h + 1) + s] + p0[j * (h + 1) + s]; if (t & 1) ok = 0; }
            if (ok) surv++;
        }
    cap = 1; while (cap < 2 * surv + 16) cap <<= 1;
    fprintf(stderr, "%llu pairs pass the spectrum test; table of %llu slots (%.2f GB)\n",
            (unsigned long long)surv, (unsigned long long)cap, cap * sizeof(entry) / 1e9);
    tab = calloc(cap, sizeof(entry));
    uint64_t bbits = 1ULL << 20; while (bbits < 8 * surv) bbits <<= 1;     /* load <= 1/8 per probe bit */
    bmask = bbits - 1; bloom = calloc(bbits / 64, 8);
    if (!tab || !bloom) { fprintf(stderr, "out of memory\n"); return 1; }
    for (int i = 0; i < n3; i++)
        for (int j = 0; j < n0; j++) {
            int ok = 1, need[64];
            for (int k = 1; k <= h && ok; k++) if (f3[i * (h + 1) + k] + f0[j * (h + 1) + k] > lim) ok = 0;
            if (!ok) continue;
            for (int s = 1; s <= h && ok; s++) {
                int t = p3[i * (h + 1) + s] + p0[j * (h + 1) + s];
                if (t & 1) ok = 0; else need[s] = -t / 2;
            }
            if (!ok) continue;
            uint64_t k0, k1; pack(need, &k0, &k1);
            if (!(k0 | k1)) k0 = 1ULL << 63;   /* never: fields are offset by 64 */
            insert(k0, k1, i, j);
        }
    fprintf(stderr, "table built (%llu entries)\n", (unsigned long long)used);

    /* side 2: A1 with -1 at 0 and (n - a1)/2 - 1 further -1s among positions 1..n-1 */
    int w = (n - a1) / 2 - 1;
    long long hits = 0, sols = 0, done = 0;
    if (w == 0) {                            /* a single -1, at 0 */
        int p[64]; paf_mask(1ULL, p);
        uint64_t k0, k1; pack(p, &k0, &k1);
        for (uint64_t idx = mix(k0, k1) & (cap - 1); tab[idx].k0 | tab[idx].k1; idx = (idx + 1) & (cap - 1))
            if (tab[idx].k0 == k0 && tab[idx].k1 == k1) {
                hits++;
                uint64_t m3 = r3[tab[idx].a3], m0 = r0[tab[idx].a0];
                if (full_check(m3, 1ULL, m0)) {
                    sols++;
                    printf("SOL "); print_seq(m3); putchar(' '); print_seq(1ULL); putchar(' '); print_seq(m0); putchar('\n');
                }
            }
    }
    /* split on the two smallest further positions f < g: A1 has -1 at 0, f and g, the rest above g */
    int npair = 0, (*pr)[2] = malloc(sizeof(int[2]) * 64 * 64);
    if (w >= 2)
        for (int f = 1; f < n; f++) for (int g = f + 1; g < n; g++) if (n - 1 - g >= w - 2) { pr[npair][0] = f; pr[npair][1] = g; npair++; }
    if (w == 1)
        for (int f = 1; f < n; f++) { pr[npair][0] = f; pr[npair][1] = 0; npair++; }
    #pragma omp parallel for schedule(dynamic, 1) reduction(+:hits, sols)
    for (int pi = 0; pi < npair; pi++) {
        int f = pr[pi][0], g = pr[pi][1];
        int k = (w >= 2) ? w - 2 : 0;          /* positions to choose in (g, n) */
        int lo = (w >= 2) ? g + 1 : n;         /* first candidate */
        int m = n - lo;
        uint64_t base = 1ULL | (1ULL << f) | (w >= 2 ? (1ULL << g) : 0);
        int c[64];
        for (int i = 0; i < k; i++) c[i] = i;
        for (;;) {
            uint64_t mask = base;
            for (int i = 0; i < k; i++) mask |= 1ULL << (lo + c[i]);
            int p[64]; paf_mask(mask, p);
            uint64_t k0, k1; pack(p, &k0, &k1);
            if (maybe(k0, k1)) {
                uint64_t idx = mix(k0, k1) & (cap - 1);
                while (tab[idx].k0 | tab[idx].k1) {
                    if (tab[idx].k0 == k0 && tab[idx].k1 == k1) {
                        hits++;
                        uint64_t m3 = r3[tab[idx].a3], m0 = r0[tab[idx].a0];
                        if (full_check(m3, mask, m0)) {
                            sols++;
                            #pragma omp critical
                            { printf("SOL "); print_seq(m3); putchar(' '); print_seq(mask); putchar(' '); print_seq(m0); putchar('\n'); fflush(stdout); }
                        }
                    }
                    idx = (idx + 1) & (cap - 1);
                }
            }
            if (k == 0) break;
            int i = k - 1;
            while (i >= 0 && c[i] == m - k + i) i--;
            if (i < 0) break;
            c[i]++;
            for (int j = i + 1; j < k; j++) c[j] = c[j - 1] + 1;
        }
        #pragma omp atomic
        done++;
        if (done % 50 == 0) {
            #pragma omp critical
            { fprintf(stderr, "%lld of %d starts done\n", done, npair); }
        }
    }
    fprintf(stderr, "STAT n=%d rows=(%d,%d,1) A3classes=%d pairs=%llu table_hits=%lld solutions=%lld\n", n, a3, a1, n3,
            (unsigned long long)used, hits, sols);
    printf("DONE n=%d rows=(%d,%d,1) solutions=%lld\n", n, a3, a1, sols);
    return 0;
}
