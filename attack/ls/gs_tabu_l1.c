/* Tabu search for Goethals-Seidel quadruples of odd length v.  F = sum_{s=1}^{v-1} S[s]^2.
 * Each step evaluates all 4v single flips (O(v^2)), takes the best non-tabu one (aspiration if it
 * improves the best F), tabu tenure uniform in [tmin, tmax]. Restart after maxstall steps without
 * improving the walk's best.  usage: ./gs_tabu v maxsteps seed [tmin tmax maxstall] */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <string.h>
static uint64_t rs;
static inline uint64_t rnd(void){ rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; }
int main(int argc, char **argv){
    int v = atoi(argv[1]); long maxsteps = atol(argv[2]); rs = 0x9E3779B97F4A7C15ULL ^ ((uint64_t)atol(argv[3]) + 1) * 0x2545F4914F6CDD1DULL;
    int tmin = argc > 4 ? atoi(argv[4]) : v/8, tmax = argc > 5 ? atoi(argv[5]) : v/4; long maxstall = argc > 6 ? atol(argv[6]) : 20L*v;
    int n = 4 * v;
    int *a = malloc(sizeof(int) * 4 * 3 * v);   /* a[i][x] for x in [0,3v) : periodic copies */
    int *S = calloc(v, sizeof(int)); long *tabu = calloc(n, sizeof(long));
    clock_t t0 = clock(); long step = 0, restarts = 0;
    while (step < maxsteps){
        for (int i = 0; i < 4; i++) for (int x = 0; x < v; x++){ int b = (rnd() & 1) ? 1 : -1; for (int c = 0; c < 3; c++) a[i*3*v + x + c*v] = b; }
        memset(S, 0, sizeof(int) * v);
        for (int i = 0; i < 4; i++) for (int s = 1; s < v; s++){ int t = 0; for (int x = 0; x < v; x++) t += a[i*3*v+v+x] * a[i*3*v+v+x+s]; S[s] += t; }
        long F = 0; for (int s = 1; s < v; s++) F += S[s] < 0 ? -S[s] : S[s];
        long best = F, stall = 0; memset(tabu, 0, sizeof(long) * n);
        while (step < maxsteps && stall < maxstall){
            long bd = 1L << 62; int bk = -1, cnt = 0;
            for (int i = 0; i < 4; i++){
                int *ai = a + i*3*v + v;
                for (int x = 0; x < v; x++){
                    int ax = ai[x]; long dF = 0;
                    for (int s = 1; s < v; s++){ int ds = -2 * ax * (ai[x+s] + ai[x-s]); int t = S[s] + ds; dF += (t < 0 ? -t : t) - (S[s] < 0 ? -S[s] : S[s]); }
                    int k = i * v + x;
                    if (tabu[k] > step && F + dF >= best) continue;
                    if (dF < bd){ bd = dF; bk = k; cnt = 1; }
                    else if (dF == bd && (rnd() % ++cnt) == 0) bk = k;
                }
            }
            if (bk < 0){ step++; continue; }
            int i = bk / v, x = bk % v; int *ai = a + i*3*v + v; int ax = ai[x];
            for (int s = 1; s < v; s++) S[s] += -2 * ax * (ai[x+s] + ai[x-s]);
            for (int c = -1; c <= 1; c++) ai[x + c*v] = -ax;
            F += bd; step++;
            tabu[bk] = step + tmin + (long)(rnd() % (tmax - tmin + 1));
            if (F < best){ best = F; stall = 0; } else stall++;
            if (F == 0){
                printf("SOLVED steps=%ld restarts=%ld time=%.2f\n", step, restarts, (double)(clock() - t0) / CLOCKS_PER_SEC);
                for (int j = 0; j < 4; j++){ for (int y = 0; y < v; y++) putchar(a[j*3*v+v+y] > 0 ? '+' : '-'); putchar('\n'); }
                return 0;
            }
        }
        restarts++;
    }
    printf("FAIL steps=%ld restarts=%ld time=%.2f\n", step, restarts, (double)(clock() - t0) / CLOCKS_PER_SEC);
    return 1;
}
