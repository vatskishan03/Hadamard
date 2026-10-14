/*
 * Simulated annealing for four +-1 sequences on Z_m with prescribed symmetry:
 *
 *     sum_k PAF_k(g) = T   for every g != 0,
 *
 * blocks of kind 'S' symmetric (x(-t) = x(t)), kind 'A' arbitrary. T = -4 is the bordered Szöllősi / Spence
 * condition on a group of order N - 1 (search/bsz.py), T = 0 the cyclic Szöllősi condition on Z_N (search/ssxx.py).
 * No multiplier symmetry is imposed. A move flips one position of an arbitrary block or one pair {t, -t} of a
 * symmetric block; the change of every PAF is computed in O(m) and the energy sum_g (S(g) - T)^2 updated exactly.
 * Restarts with a geometric cooling schedule; prints the four sequences of the first zero-energy state.
 *
 *     sls_sz m KINDS T seed seconds [T0 T1 sweeps]       e.g.  sls_sz 106 SSAA -4 1 3600
 */
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static uint64_t rs;
static inline uint64_t rnd(void) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; }
static inline double urand(void) { return (rnd() >> 11) * (1.0 / 9007199254740992.0); }

int m, T;
char kinds[5];
int x[4][4096];
long S[4096];                 /* S[g] = sum_k PAF_k(g) */
long dl[4096];

static long paf(int *v, int g) { long s = 0; for (int t = 0; t < m; t++) s += v[t] * v[(t + g) % m]; return s; }

static void full(void) {
    for (int g = 0; g < m; g++) { S[g] = 0; for (int k = 0; k < 4; k++) S[g] += paf(x[k], g); }
}

static long energy(void) { long e = 0; for (int g = 1; g < m; g++) { long d = S[g] - T; e += d * d; } return e; }

/* delta of PAF at every g != 0 when position t of v flips (v unchanged here); accumulate into dl */
static inline void flip_delta(int *v, int t) {
    int xt = v[t];
    for (int g = 1; g < m; g++) {
        int a = t + g; if (a >= m) a -= m;
        int b = t - g; if (b < 0) b += m;
        dl[g] += -2L * xt * (v[a] + v[b]);
    }
}

int main(int argc, char **argv) {
    if (argc < 6) { fprintf(stderr, "usage: sls_sz m KINDS T seed seconds [T0 T1 sweeps]\n"); return 1; }
    m = atoi(argv[1]); strncpy(kinds, argv[2], 4); T = atoi(argv[3]);
    rs = 0x9E3779B97F4A7C15ULL ^ (uint64_t)atoll(argv[4]) * 0x2545F4914F6CDD1DULL; if (!rs) rs = 1;
    double secs = atof(argv[5]);
    double T0 = argc > 6 ? atof(argv[6]) : 40.0, T1 = argc > 7 ? atof(argv[7]) : 0.5;
    long sweeps = argc > 8 ? atol(argv[8]) : 4000;
    /* move list: (block, position) with position the orbit representative for symmetric blocks */
    int nmv = 0, mvk[4 * 4096], mvt[4 * 4096];
    for (int k = 0; k < 4; k++)
        for (int t = 0; t < m; t++) {
            if (kinds[k] == 'S' && t > (m - t) % m) continue;
            mvk[nmv] = k; mvt[nmv] = t; nmv++;
        }
    time_t start = time(NULL);
    long restarts = 0, best_all = -1;
    while (difftime(time(NULL), start) < secs) {
        restarts++;
        for (int k = 0; k < 4; k++)
            for (int t = 0; t < m; t++) {
                if (kinds[k] == 'S' && t > (m - t) % m) { x[k][t] = x[k][(m - t) % m]; continue; }
                x[k][t] = (rnd() & 1) ? 1 : -1;
            }
        full();
        long E = energy(), best = E;
        long steps = sweeps * (long)nmv;
        double lt0 = log(T0), lt1 = log(T1);
        for (long s = 0; s < steps; s++) {
            double temp = exp(lt0 + (lt1 - lt0) * (double)s / steps);
            int i = rnd() % nmv, k = mvk[i], t = mvt[i];
            int t2 = (m - t) % m;
            memset(dl, 0, sizeof(long) * m);
            flip_delta(x[k], t);
            if (kinds[k] == 'S' && t2 != t) {      /* second flip, applied after the first */
                x[k][t] = -x[k][t];
                flip_delta(x[k], t2);
                x[k][t] = -x[k][t];
            }
            long dE = 0;
            for (int g = 1; g < m; g++) { long d = S[g] - T; dE += dl[g] * (2 * d + dl[g]); }
            if (dE <= 0 || urand() < exp(-dE / temp)) {
                x[k][t] = -x[k][t];
                if (kinds[k] == 'S' && t2 != t) x[k][t2] = -x[k][t2];
                for (int g = 1; g < m; g++) S[g] += dl[g];
                E += dE;
                if (E < best) best = E;
                if (E == 0) {
                    printf("FOUND m=%d kinds=%s T=%d restart=%ld step=%ld seconds=%.0f\n", m, kinds, T, restarts, s,
                           difftime(time(NULL), start));
                    for (int kk = 0; kk < 4; kk++) {
                        for (int tt = 0; tt < m; tt++) putchar(x[kk][tt] > 0 ? '+' : '-');
                        putchar('\n');
                    }
                    fflush(stdout);
                    return 0;
                }
            }
        }
        if (best_all < 0 || best < best_all) best_all = best;
        if (restarts % 20 == 0) { fprintf(stderr, "restart %ld best %ld (overall %ld) %.0fs\n", restarts, best, best_all,
                                          difftime(time(NULL), start)); }
    }
    printf("NONE restarts=%ld best=%ld\n", restarts, best_all);
    return 1;
}
