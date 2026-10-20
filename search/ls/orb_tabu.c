/* Tabu search for Goethals-Seidel quadruples constant on the orbits of a multiplier group.
 * Input file (int32): k C, w[C], M[C][k][k] (M_c = S_c + S_c^T with zero diagonal), D[C] (diag part:
 * sum_O S_c[O][O], constant).  V_c = sum_i (D_c + y_i^T S_c y_i ... ) computed via fields.
 * Energy F = sum_c w_c V_c^2,  V_c = sum_i PAF_i(s_c) (target 0).
 * usage: ./orb_tabu file maxsteps seed tmin tmax maxstall */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
static uint64_t rs;
static inline uint64_t rnd(void){ rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; }
int main(int argc, char **argv){
    FILE *f = fopen(argv[1], "rb"); int k, C;
    if (fread(&k, 4, 1, f) != 1 || fread(&C, 4, 1, f) != 1) return 2;
    int *w = malloc(4 * C); if (fread(w, 4, C, f) != (size_t)C) return 2;
    int *M = malloc(4L * C * k * k); if (fread(M, 4, (size_t)C * k * k, f) != (size_t)C * k * k) return 2;
    int *Dg = malloc(4 * C); if (fread(Dg, 4, C, f) != (size_t)C) return 2;
    long maxsteps = atol(argv[2]); rs = 0x9E3779B97F4A7C15ULL ^ ((uint64_t)atol(argv[3]) + 7) * 0x2545F4914F6CDD1DULL;
    int tmin = atoi(argv[4]), tmax = atoi(argv[5]); long maxstall = atol(argv[6]);
    int n = 4 * k;
    int *y = malloc(4 * n);                 /* y[i*k + O] */
    int *g = malloc(4L * 4 * C * k);        /* g[(i*C + c)*k + O] = sum_O' M_c[O][O'] y_i[O'] */
    long *V = malloc(8 * C); long *tabu = malloc(8 * n);
    clock_t t0 = clock(); long step = 0, restarts = 0, bestever = 1L << 62;
    while (step < maxsteps){
        for (int j = 0; j < n; j++) y[j] = (rnd() & 1) ? 1 : -1;
        for (int i = 0; i < 4; i++) for (int c = 0; c < C; c++) for (int O = 0; O < k; O++){
            long t = 0; const int *row = M + ((long)c * k + O) * k; for (int P = 0; P < k; P++) t += row[P] * y[i*k + P];
            g[(i*C + c)*k + O] = (int)t; }
        for (int c = 0; c < C; c++){ long t = 0; for (int i = 0; i < 4; i++){ long q = 0; for (int O = 0; O < k; O++) q += (long)y[i*k+O] * g[(i*C+c)*k+O]; t += Dg[c] + q / 2; } V[c] = t; }
        long F = 0; for (int c = 0; c < C; c++) F += (long)w[c] * V[c] * V[c];
        long best = F, stall = 0; memset(tabu, 0, 8 * n);
        while (step < maxsteps && stall < maxstall){
            long bd = 1L << 62; int bj = -1, cnt = 0;
            for (int i = 0; i < 4; i++) for (int O = 0; O < k; O++){
                int j = i*k + O; int yo = y[j]; long dF = 0;
                for (int c = 0; c < C; c++){ long d = -2L * yo * g[(i*C + c)*k + O]; dF += w[c] * (2 * V[c] * d + d * d); }
                if (tabu[j] > step && F + dF >= best) continue;
                if (dF < bd){ bd = dF; bj = j; cnt = 1; } else if (dF == bd && (rnd() % ++cnt) == 0) bj = j;
            }
            if (bj < 0){ step++; continue; }
            int i = bj / k, O = bj % k; int yo = y[bj];
            for (int c = 0; c < C; c++){
                V[c] += -2L * yo * g[(i*C + c)*k + O];
                const int *col = M + (long)c * k * k;          /* M symmetric: column O = row O */
                int *gc = g + (i*C + c)*k; const int *row = col + (long)O * k;
                for (int P = 0; P < k; P++) gc[P] += -2 * yo * row[P];
            }
            y[bj] = -yo; F += bd; step++;
            tabu[bj] = step + tmin + (long)(rnd() % (tmax - tmin + 1));
            if (F < best){ best = F; stall = 0; } else stall++;
            if (F < bestever) bestever = F;
            if (F == 0){
                printf("SOLVED steps=%ld restarts=%ld time=%.2f\n", step, restarts, (double)(clock() - t0) / CLOCKS_PER_SEC);
                for (int j = 0; j < n; j++) printf("%d ", y[j]); printf("\n"); return 0;
            }
        }
        restarts++;
    }
    printf("FAIL steps=%ld restarts=%ld bestever=%ld time=%.2f\n", step, restarts, bestever, (double)(clock() - t0) / CLOCKS_PER_SEC);
    return 1;
}
