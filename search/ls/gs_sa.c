/* Simulated annealing for Goethals-Seidel quadruples of odd length v (no symmetry).
 * Energy F = sum_{s=1}^{v-1} S[s]^2,  S[s] = sum_i PAF_i(s).  Single-entry flips, Metropolis at T.
 * usage: ./gs_sa v T maxflips seed            prints SOLVED flips time  and the 4 sequences */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
static uint64_t rs;
static inline uint64_t rnd(void){ rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; }
int main(int argc, char **argv){
    int v = atoi(argv[1]); double T = atof(argv[2]); long maxf = atol(argv[3]); rs = 0x9E3779B97F4A7C15ULL ^ (uint64_t)atol(argv[4]) * 0x2545F4914F6CDD1DULL;
    int *a = malloc(sizeof(int) * 4 * 2 * v);           /* doubled for wraparound: a[i][x] for x in [0,2v) */
    long *S = calloc(v, sizeof(long)); int *d = malloc(sizeof(int) * v);
    for (int i = 0; i < 4; i++) for (int x = 0; x < v; x++){ int b = (rnd() & 1) ? 1 : -1; a[i*2*v + x] = a[i*2*v + x + v] = b; }
    for (int i = 0; i < 4; i++) for (int s = 1; s < v; s++){ long t = 0; for (int x = 0; x < v; x++) t += a[i*2*v+x] * a[i*2*v+x+s]; S[s] += t; }
    long F = 0; for (int s = 1; s < v; s++) F += S[s]*S[s];
    double *acc = malloc(sizeof(double) * 4096);       /* acceptance threshold table for dF in units of 16 */
    for (int k = 0; k < 4096; k++) acc[k] = exp(-16.0 * k / T);
    clock_t t0 = clock();
    for (long f = 0; f < maxf; f++){
        int i = rnd() % 4, x = rnd() % v; int *ai = a + i*2*v; int ax = ai[x];
        long dF = 0;
        for (int s = 1; s < v; s++){
            int u = ai[x + s] + ai[x - s + v];             /* a(x+s) + a(x-s) */
            int ds = -2 * ax * u; d[s] = ds;
            dF += 2 * S[s] * ds + ds * ds;
        }
        int ok = dF <= 0;
        if (!ok){ long k = dF / 16; ok = k < 4096 && (rnd() >> 11) * (1.0 / 9007199254740992.0) < acc[k]; }
        if (ok){
            for (int s = 1; s < v; s++) S[s] += d[s];
            ai[x] = ai[x + v] = -ax; F += dF;
            if (F == 0){
                printf("SOLVED flips=%ld time=%.2f\n", f, (double)(clock() - t0) / CLOCKS_PER_SEC);
                for (int j = 0; j < 4; j++){ for (int y = 0; y < v; y++) putchar(a[j*2*v+y] > 0 ? '+' : '-'); putchar('\n'); }
                return 0;
            }
        }
    }
    printf("FAIL F=%ld time=%.2f\n", F, (double)(clock() - t0) / CLOCKS_PER_SEC);
    return 1;
}
