/* Replicated simulated annealing ("local entropy", Baldassi et al.) for Goethals-Seidel quadruples of length v.
 * R replicas, each with energy F_a = sum_s S_a[s]^2.  Coupling: -gamma * sum_a <x_a, c> where c = sign of the
 * replica average (recomputed lazily) -> pulls replicas into dense regions.  Metropolis single flips,
 * temperature T fixed, gamma ramped geometrically.  usage: ./gs_rsa v R T gamma0 gammamul steps seed */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <time.h>
static uint64_t rs;
static inline uint64_t rnd(void){ rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; }
static inline double unif(void){ return (rnd() >> 11) * (1.0 / 9007199254740992.0); }
int main(int argc, char **argv){
    int v = atoi(argv[1]), R = atoi(argv[2]); double T = atof(argv[3]), g = atof(argv[4]), gm = atof(argv[5]);
    long steps = atol(argv[6]); rs = 0x9E3779B97F4A7C15ULL ^ (atol(argv[7]) + 11) * 0x2545F4914F6CDD1DULL;
    int n = 4 * v;
    int *x = malloc(sizeof(int) * R * 4 * 3 * v); int *S = calloc((size_t)R * v, sizeof(int)); long *F = calloc(R, sizeof(long));
    int *sum = calloc(n, sizeof(int));                 /* sum over replicas of x_a(i,y) */
    #define X(a,i,y) x[(((long)(a) * 4 + (i)) * 3 * v) + v + (y)]
    for (int a = 0; a < R; a++){
        for (int i = 0; i < 4; i++) for (int y = 0; y < v; y++){ int b = (rnd() & 1) ? 1 : -1; for (int c = -1; c <= 1; c++) X(a,i,y + c*v) = b; sum[i*v + y] += b; }
        for (int i = 0; i < 4; i++) for (int s = 1; s < v; s++){ int t = 0; for (int y = 0; y < v; y++) t += X(a,i,y) * X(a,i,y+s); S[(long)a*v + s] += t; }
        for (int s = 1; s < v; s++) F[a] += (long)S[(long)a*v+s] * S[(long)a*v+s];
    }
    clock_t t0 = clock(); long best = 1L << 62;
    for (long it = 0; it < steps; it++){
        if (it % 100000 == 0) g *= gm;
        int a = rnd() % R, i = rnd() % 4, y = rnd() % v; int xo = X(a,i,y);
        long dF = 0; int *Sa = S + (long)a * v;
        for (int s = 1; s < v; s++){ int ds = -2 * xo * (X(a,i,y+s) + X(a,i,y-s)); dF += 2L * Sa[s] * ds + ds * ds; }
        /* coupling: energy -g * x * (sum - x) (attraction to other replicas' average) */
        double dC = 2.0 * g * xo * (sum[i*v + y] - xo);
        double dE = dF + dC;
        if (dE <= 0 || unif() < exp(-dE / T)){
            for (int s = 1; s < v; s++) Sa[s] += -2 * xo * (X(a,i,y+s) + X(a,i,y-s));
            for (int c = -1; c <= 1; c++) X(a,i,y + c*v) = -xo;
            sum[i*v + y] -= 2 * xo; F[a] += dF;
            if (F[a] < best) best = F[a];
            if (F[a] == 0){
                printf("SOLVED v=%d it=%ld replica=%d time=%.2f\n", v, it, a, (double)(clock() - t0) / CLOCKS_PER_SEC); return 0; }
        }
    }
    printf("FAIL v=%d best=%ld time=%.2f\n", v, best, (double)(clock() - t0) / CLOCKS_PER_SEC);
    return 1;
}
