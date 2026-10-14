/* Tabu search for base sequences BS(m+1, m): a, b of length m+1, c, d of length m (+-1),
 * N_a(s) + N_b(s) + N_c(s) + N_d(s) = 0 for s = 1..m (aperiodic autocorrelation).
 * Energy F = sum_s S[s]^2.  Steepest descent over all single flips with tabu; restarts after maxstall.
 * usage: ./bs_tabu m maxsteps seed tmin tmax maxstall          prints SOLVED + the four sequences */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
static uint64_t rs;
static inline uint64_t rnd(void){ rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; }
int main(int argc, char **argv){
    int m = atoi(argv[1]); long maxsteps = atol(argv[2]); rs = 0x9E3779B97F4A7C15ULL ^ ((uint64_t)atol(argv[3]) + 3) * 0x2545F4914F6CDD1DULL;
    int tmin = atoi(argv[4]), tmax = atoi(argv[5]); long maxstall = atol(argv[6]);
    int len[4] = {m + 1, m + 1, m, m}; int W = m + 1, PAD = m + 2;
    int *buf = calloc(4 * (W + 2 * PAD), sizeof(int));        /* zero padding on both sides */
    int *q[4]; for (int i = 0; i < 4; i++) q[i] = buf + i * (W + 2 * PAD) + PAD;
    int *S = calloc(m + 2, sizeof(int)); long *tabu = calloc(4 * W, sizeof(long));
    clock_t t0 = clock(); long step = 0, restarts = 0, bestever = 1L << 62;
    while (step < maxsteps){
        for (int i = 0; i < 4; i++) for (int x = 0; x < len[i]; x++) q[i][x] = (rnd() & 1) ? 1 : -1;
        memset(S, 0, sizeof(int) * (m + 2));
        for (int i = 0; i < 4; i++) for (int s = 1; s <= m; s++){ int t = 0; for (int x = 0; x + s < len[i]; x++) t += q[i][x] * q[i][x + s]; S[s] += t; }
        long F = 0; for (int s = 1; s <= m; s++) F += (long)S[s] * S[s];
        long best = F, stall = 0; memset(tabu, 0, sizeof(long) * 4 * W);
        while (step < maxsteps && stall < maxstall){
            long bd = 1L << 62; int bi = -1, bx = -1, cnt = 0;
            for (int i = 0; i < 4; i++){
                int *a = q[i];
                for (int x = 0; x < len[i]; x++){
                    int ax = a[x]; long dF = 0;
                    for (int s = 1; s <= m; s++){ int ds = -2 * ax * (a[x + s] + a[x - s]); dF += 2L * S[s] * ds + ds * ds; }
                    long k = (long)i * W + x;
                    if (tabu[k] > step && F + dF >= best) continue;
                    if (dF < bd){ bd = dF; bi = i; bx = x; cnt = 1; } else if (dF == bd && (rnd() % ++cnt) == 0){ bi = i; bx = x; }
                }
            }
            if (bi < 0){ step++; continue; }
            int *a = q[bi]; int ax = a[bx];
            for (int s = 1; s <= m; s++) S[s] += -2 * ax * (a[bx + s] + a[bx - s]);
            a[bx] = -ax; F += bd; step++;
            tabu[(long)bi * W + bx] = step + tmin + (long)(rnd() % (tmax - tmin + 1));
            if (F < best){ best = F; stall = 0; } else stall++;
            if (F < bestever) bestever = F;
            if (F == 0){
                printf("SOLVED m=%d steps=%ld restarts=%ld time=%.2f\n", m, step, restarts, (double)(clock() - t0) / CLOCKS_PER_SEC);
                for (int i = 0; i < 4; i++){ for (int x = 0; x < len[i]; x++) putchar(q[i][x] > 0 ? '+' : '-'); putchar('\n'); }
                return 0;
            }
        }
        restarts++;
    }
    printf("FAIL m=%d steps=%ld restarts=%ld bestever=%ld time=%.2f\n", m, step, restarts, bestever, (double)(clock() - t0) / CLOCKS_PER_SEC);
    return 1;
}
