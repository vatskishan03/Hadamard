/* Tabu search for cyclic sequence families without multiplier symmetry:
 *   sum_k w_k PAF_{X_k}(s) = T   for s = 1 .. (v-1)/2,
 * each block X_k a +-1 sequence of odd length v of type
 *   f  free,
 *   s  symmetric (X(-x) = X(x)),
 *   k  skew type (X(0) = 1, X(-x) = -X(x) for x != 0),
 * with weight w_k in {1, 2}. Families:
 *   skew Goethals-Seidel  "k1,f1,f1,f1"  T = 0  (skew-Hadamard matrix of order 4v)
 *   propus                "s1,f2,f1"     T = 0  (symmetric Hadamard matrix of order 4v)
 *   Williamson            "s1,s1,s1,s1"  T = 0
 *   Goethals-Seidel       "f1,f1,f1,f1"  T = 0  (Hadamard matrix of order 4v)
 *   Legendre pair         "f1,f1"        T = -2; D-optimal design "f1,f1" T = 2.
 *
 * Energy E = sum_s (P(s) - T)^2. Moves: one entry of a free block, one pair {x, -x} of a symmetric or skew block
 * (or X(0) of a symmetric block). Each step evaluates every move exactly and makes the best one that is not tabu
 * (tabu for TT..2TT steps; aspiration on a new best); restart after R steps without a new best.
 *
 *   gs_ls v types T seed maxsec [R] [TT] [nsol]
 * prints "SOL <X_1> ... <X_k>" (+/- strings) per solution and "STAT ..." at the end.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>

static uint64_t rs;
static inline uint64_t rnd(void) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; }

#define MAXB 8
static int v, h, nb, T;
static char typ[MAXB];
static int w[MAXB];
static int *X[MAXB];
static int *P, *dP;
static int M, *mb, *mx;            /* move -> block, position */

static void full_paf(void) {
    for (int s = 1; s <= h; s++) {
        long t = 0;
        for (int k = 0; k < nb; k++) {
            long u = 0;
            for (int y = 0; y < v; y++) u += X[k][y] * X[k][(y + s) % v];
            t += w[k] * u;
        }
        P[s] = (int)t;
    }
}

static long energy(void) { long e = 0; for (int s = 1; s <= h; s++) e += (long)(P[s] - T) * (P[s] - T); return e; }

static void delta(int m, int *out) {
    int k = mb[m], x = mx[m], *a = X[k], wk = w[k];
    if (typ[k] == 'f') {
        int ax = a[x];
        for (int s = 1; s <= h; s++) {
            int p = x + s; if (p >= v) p -= v;
            int q = x - s; if (q < 0) q += v;
            out[s] = -2 * wk * ax * (a[p] + a[q]);
        }
    } else if (x == 0) {                          /* symmetric block, centre entry */
        int ax = a[0];
        for (int s = 1; s <= h; s++) out[s] = -4 * wk * ax * a[s];
    } else {                                      /* pair {x, -x}: d(x) = -2a(x), d(-x) = -2a(-x) */
        int nx = v - x, ax = a[x], anx = a[nx];
        int t2 = (2 * x) % v; if (t2 > h) t2 = v - t2;
        int second = 4 * ax * anx;                /* d(x) d(-x) */
        for (int s = 1; s <= h; s++) {
            int xp = x + s; if (xp >= v) xp -= v;
            int xm = x - s; if (xm < 0) xm += v;
            int np = nx + s; if (np >= v) np -= v;
            int nm = nx - s; if (nm < 0) nm += v;
            /* first order: sum_y a(y) d(y+s) + d(y) a(y+s) = d(x)(a(x-s) + a(x+s)) + d(-x)(a(-x-s) + a(-x+s)),
               where a() is the old sequence; second order d(x)d(-x) at s = +-2x (both orientations fold to t2) */
            int f = -2 * ax * (a[xm] + a[xp]) - 2 * anx * (a[nm] + a[np]);
            out[s] = wk * (f + (s == t2 ? second : 0));
        }
    }
}

static void apply(int m) {
    int k = mb[m], x = mx[m];
    X[k][x] = -X[k][x];
    if (typ[k] != 'f' && x) X[k][v - x] = -X[k][v - x];
}

static void randomize(void) {
    for (int k = 0; k < nb; k++) {
        int *a = X[k];
        if (typ[k] == 'f') { for (int x = 0; x < v; x++) a[x] = (rnd() & 1) ? 1 : -1; }
        else if (typ[k] == 's') { for (int x = 0; x <= h; x++) { a[x] = (rnd() & 1) ? 1 : -1; if (x) a[v - x] = a[x]; } }
        else { a[0] = 1; for (int x = 1; x <= h; x++) { a[x] = (rnd() & 1) ? 1 : -1; a[v - x] = -a[x]; } }
    }
    full_paf();
}

static void put(const int *x) { for (int i = 0; i < v; i++) putchar(x[i] > 0 ? '+' : '-'); }

int main(int argc, char **argv) {
    if (argc < 6) { fprintf(stderr, "usage: gs_ls v types T seed maxsec [R] [TT] [nsol]\n"); return 1; }
    v = atoi(argv[1]);
    if (v < 3 || v % 2 == 0) { fprintf(stderr, "v must be odd\n"); return 1; }
    h = (v - 1) / 2;
    nb = 0;
    for (char *p = argv[2]; *p && nb < MAXB; ) {
        typ[nb] = p[0]; w[nb] = p[1] - '0';
        if (!strchr("fsk", typ[nb]) || w[nb] < 1 || w[nb] > 2) { fprintf(stderr, "bad block type\n"); return 1; }
        nb++; p += 2; if (*p == ',') p++;
    }
    T = atoi(argv[3]);
    rs = 0x9E3779B97F4A7C15ULL ^ (uint64_t)atoll(argv[4]) * 0x2545F4914F6CDD1DULL; if (!rs) rs = 1;
    double maxsec = atof(argv[5]);
    long R = argc > 6 ? atol(argv[6]) : 20000;
    int TT = argc > 7 ? atoi(argv[7]) : 0;
    int nsol = argc > 8 ? atoi(argv[8]) : 1;
    M = 0;
    mb = calloc(nb * v, sizeof(int)); mx = calloc(nb * v, sizeof(int));
    for (int k = 0; k < nb; k++) {
        X[k] = calloc(v, sizeof(int));
        if (typ[k] == 'f') for (int x = 0; x < v; x++) { mb[M] = k; mx[M] = x; M++; }
        else if (typ[k] == 's') for (int x = 0; x <= h; x++) { mb[M] = k; mx[M] = x; M++; }
        else for (int x = 1; x <= h; x++) { mb[M] = k; mx[M] = x; M++; }
    }
    if (TT <= 0) TT = M / 8 + 1;
    P = calloc(h + 1, sizeof(int)); dP = calloc(h + 1, sizeof(int));
    long *tabu = calloc(M, sizeof(long));
    int *cand = calloc(M, sizeof(int));
    clock_t t0 = clock();
    long step = 0, restarts = 0, found = 0, sincebest = 0, E, best;
    randomize(); E = energy(); best = E;
    while (found < nsol) {
        if ((step & 1023) == 0 && (double)(clock() - t0) / CLOCKS_PER_SEC > maxsec) break;
        long bdE = 1L << 62; int nc = 0;
        for (int m = 0; m < M; m++) {
            delta(m, dP);
            long dE = 0;
            for (int s = 1; s <= h; s++) dE += (long)dP[s] * (2L * (P[s] - T) + dP[s]);
            if (tabu[m] > step && E + dE >= best) continue;
            if (dE < bdE) { bdE = dE; nc = 0; cand[nc++] = m; }
            else if (dE == bdE) cand[nc++] = m;
        }
        if (!nc) { step++; continue; }
        int m = cand[rnd() % nc];
        delta(m, dP);
        for (int s = 1; s <= h; s++) P[s] += dP[s];
        apply(m);
        E += bdE;
        tabu[m] = step + TT + (long)(rnd() % (TT + 1));
        step++;
        if (E < best) { best = E; sincebest = 0; } else sincebest++;
        if (E == 0) {
            full_paf();
            if (energy() != 0) { fprintf(stderr, "internal error: energy mismatch\n"); return 2; }
            printf("SOL");
            for (int k = 0; k < nb; k++) { putchar(' '); put(X[k]); }
            printf(" step %ld restart %ld\n", step, restarts);
            fflush(stdout);
            found++;
            randomize(); E = energy(); best = E; sincebest = 0; restarts++;
            memset(tabu, 0, M * sizeof(long));
            continue;
        }
        if (sincebest > R) {
            randomize(); E = energy(); best = E; sincebest = 0; restarts++;
            memset(tabu, 0, M * sizeof(long));
        }
    }
    printf("STAT v %d types %s steps %ld restarts %ld found %ld sec %.1f\n", v, argv[2], step, restarts, found,
           (double)(clock() - t0) / CLOCKS_PER_SEC);
    return 0;
}
