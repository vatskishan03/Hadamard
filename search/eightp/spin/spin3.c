// spin3.c -- exact search of the d = 3 block-symmetric Spence structure (A, B, muB, mu^2B) over Z_v (v odd), Gray-code enumeration of B.
// input file: C OB NA
//   size_B[OB]
//   SB[c][a][b]  (C x OB x OB, ints)      W_c(x) = sum_{a,b} SB[c][a][b] x_a x_b
//   NA lines: C ints = the target vector T_c = -4 - PAF_A(g_c) of one admissible A (row sum +1)
// COUNT = number of pairs (B, distinct A PAF vector): the NA lines are the DISTINCT target vectors, so admissible A with equal PAF vectors count once (zero / nonzero is exact).
// output: "COUNT n" and up to 20 solutions "SOL <gray index of B> <A line index>" ; B = x with x_i = 1 - 2*bit_i of the index.
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
typedef long long ll; typedef unsigned long long u64;
static int C, OB, NA;
int main(int argc, char **argv) {
    FILE *f = fopen(argv[1], "r"); if (!f) return 1;
    if (fscanf(f, "%d %d %d", &C, &OB, &NA) != 3) return 1;
    int *size = malloc(sizeof(int) * OB);
    for (int i = 0; i < OB; i++) if (fscanf(f, "%d", &size[i]) != 1) return 1;
    int *SB = malloc(sizeof(int) * C * OB * OB);
    for (long i = 0; i < (long)C * OB * OB; i++) if (fscanf(f, "%d", &SB[i]) != 1) return 1;
    ll *T = malloc(sizeof(ll) * NA * C);
    for (long i = 0; i < (long)NA * C; i++) if (fscanf(f, "%lld", &T[i]) != 1) return 1;
    fclose(f);
    u64 coef[64]; u64 st = 88172645463325252ULL;
    for (int c = 0; c < C; c++) { st ^= st << 13; st ^= st >> 7; st ^= st << 17; coef[c] = st | 1; }
    // hash table of A targets
    size_t tsz = 1; while (tsz < (size_t)NA * 4) tsz <<= 1;
    u64 *key = calloc(tsz, 8); int *val = malloc(sizeof(int) * tsz); char *used = calloc(tsz, 1);
    for (int a = 0; a < NA; a++) {
        u64 h = 0; for (int c = 0; c < C; c++) h += coef[c] * (u64)T[(long)a * C + c];
        size_t p = (h * 0x9E3779B97F4A7C15ULL) >> 20; p &= tsz - 1;
        while (used[p] && key[p] != h) p = (p + 1) & (tsz - 1);
        if (!used[p]) { used[p] = 1; key[p] = h; val[p] = a; }       // duplicate A vectors: keep the first (counted once)
    }
    // state: x (+-1), u[c][b] = sum_a (SB[c][a][b] + SB[c][b][a]) x_a ; W[c] = x^T SB x ; rho = sum size x
    int *x = malloc(sizeof(int) * OB); for (int i = 0; i < OB; i++) x[i] = 1;
    ll *W = calloc(C, sizeof(ll)); ll *u = calloc((size_t)C * OB, sizeof(ll));
    int *SS = malloc(sizeof(int) * C * OB * OB);        // symmetrised SB + SB^T
    for (int c = 0; c < C; c++) for (int a = 0; a < OB; a++) for (int b = 0; b < OB; b++) SS[((size_t)c * OB + a) * OB + b] = SB[((size_t)c * OB + a) * OB + b] + SB[((size_t)c * OB + b) * OB + a];
    for (int c = 0; c < C; c++) {
        ll w = 0; for (int a = 0; a < OB; a++) for (int b = 0; b < OB; b++) w += SB[((size_t)c * OB + a) * OB + b];
        W[c] = w;
        for (int b = 0; b < OB; b++) { ll s = 0; for (int a = 0; a < OB; a++) s += (ll)SS[((size_t)c * OB + a) * OB + b] * x[a]; u[(size_t)c * OB + b] = s; }
    }
    ll rho = 0; for (int i = 0; i < OB; i++) rho += size[i];
    u64 hW = 0; for (int c = 0; c < C; c++) hW += coef[c] * (u64)W[c];
    ll count = 0; int printed = 0;
    u64 N = 1ULL << OB;
    for (u64 g = 0; g < N; g++) {
        if (rho == 1) {
            size_t p = (hW * 0x9E3779B97F4A7C15ULL) >> 20; p &= tsz - 1;
            while (used[p]) {
                if (key[p] == hW) {
                    // verify exactly
                    int a = val[p], ok = 1;
                    for (int c = 0; c < C; c++) if (W[c] != T[(long)a * C + c]) { ok = 0; break; }
                    if (ok) { count++; if (printed < 20) { printed++; u64 gray = g ^ (g >> 1); printf("SOL %llu %d\n", gray, a); } }
                    break;
                }
                p = (p + 1) & (tsz - 1);
            }
        }
        if (g + 1 == N) break;
        int i = __builtin_ctzll(g + 1);            // bit flipped by the Gray code step g -> g+1
        // flip x[i]
        int xi = x[i];
        for (int c = 0; c < C; c++) {
            const int *Sc = SS + ((size_t)c * OB + i) * OB;
            ll s = u[(size_t)c * OB + i] - (ll)Sc[i] * xi;     // sum_{b != i} SS[i][b] x_b   (u[c][i] = sum_a SS[a][i] x_a ; SS symmetric)
            ll dW = -2LL * xi * s;
            W[c] += dW; hW += coef[c] * (u64)dW;
            ll d = -2LL * xi;                                  // change of x_i
            ll *uc = u + (size_t)c * OB;
            for (int b = 0; b < OB; b++) uc[b] += (ll)SS[((size_t)c * OB + i) * OB + b] * d;
        }
        rho += -2LL * xi * size[i]; x[i] = -xi;
    }
    printf("COUNT %lld\n", count);
    return 0;
}
