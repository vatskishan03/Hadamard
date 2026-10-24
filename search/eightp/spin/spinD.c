// spinD.c -- block-symmetric Spence structures, ONE H-invariant list X with blocks given by multiplier images.
//   mode 1 (d = 4): blocks (X, mu X, mu^2 X, mu^3 X): W(x) = sum_j PAF_X(mu^j g) = -4 for all classes g   (single list, count of x with row sum +1)
//   mode 2 (d = 2): blocks (X, mu X, Y, mu Y):       W(x) + W(y) = -4,  W(x) = PAF_X(g) + PAF_X(mu g)    (self-join of the list; ordered pairs (x,y))
// input file: mode C OB / size[OB] / S[c][a][b] (W_c(x) = sum_ab S[c][a][b] x_a x_b) / T[c] (target, -4)
// output: up to 20 lines "SOL gray [gray2]" and "COUNT n"; mode 2 counts are 64-bit-hash matches (the driver verifies the printed ones exactly)
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
typedef long long ll; typedef unsigned long long u64;
typedef struct { u64 h, g; } P;
static int cmpP(const void *a, const void *b) { const P *p = a, *q = b; return (p->h > q->h) - (p->h < q->h); }
int main(int argc, char **argv) {
    FILE *f = fopen(argv[1], "r"); if (!f) return 1;
    int mode, C, OB; if (fscanf(f, "%d %d %d", &mode, &C, &OB) != 3) return 1;
    int *size = malloc(sizeof(int) * OB); for (int i = 0; i < OB; i++) if (fscanf(f, "%d", &size[i]) != 1) return 1;
    int *SB = malloc(sizeof(int) * (size_t)C * OB * OB); for (size_t i = 0; i < (size_t)C * OB * OB; i++) if (fscanf(f, "%d", &SB[i]) != 1) return 1;
    ll *T = malloc(sizeof(ll) * C); for (int c = 0; c < C; c++) if (fscanf(f, "%lld", &T[c]) != 1) return 1;
    fclose(f);
    u64 coef[256]; u64 st = 88172645463325252ULL;
    for (int c = 0; c < C; c++) { st ^= st << 13; st ^= st >> 7; st ^= st << 17; coef[c] = st | 1; }
    u64 hT = 0; for (int c = 0; c < C; c++) hT += coef[c] * (u64)T[c];
    int *SS = malloc(sizeof(int) * (size_t)C * OB * OB);
    for (int c = 0; c < C; c++) for (int a = 0; a < OB; a++) for (int b = 0; b < OB; b++) SS[((size_t)c * OB + a) * OB + b] = SB[((size_t)c * OB + a) * OB + b] + SB[((size_t)c * OB + b) * OB + a];
    int *x = malloc(sizeof(int) * OB); for (int i = 0; i < OB; i++) x[i] = 1;
    ll *W = calloc(C, sizeof(ll)); ll *u = calloc((size_t)C * OB, sizeof(ll));
    for (int c = 0; c < C; c++) {
        ll w = 0; for (int a = 0; a < OB; a++) for (int b = 0; b < OB; b++) w += SB[((size_t)c * OB + a) * OB + b];
        W[c] = w;
        for (int b = 0; b < OB; b++) { ll s = 0; for (int a = 0; a < OB; a++) s += (ll)SS[((size_t)c * OB + a) * OB + b] * x[a]; u[(size_t)c * OB + b] = s; }
    }
    ll rho = 0; for (int i = 0; i < OB; i++) rho += size[i];
    u64 hW = 0; for (int c = 0; c < C; c++) hW += coef[c] * (u64)W[c];
    u64 N = 1ULL << OB;
    ll count = 0; int printed = 0;
    P *L = NULL; size_t nL = 0, capL = 0;
    for (u64 g = 0; g < N; g++) {
        if (rho == 1) {
            if (mode == 1) {
                if (hW == hT) {
                    int ok = 1; for (int c = 0; c < C; c++) if (W[c] != T[c]) { ok = 0; break; }
                    if (ok) { count++; if (printed < 20) { printed++; printf("SOL %llu\n", g ^ (g >> 1)); } }
                }
            } else {
                if (nL == capL) { capL = capL ? capL * 2 : (1 << 16); L = realloc(L, sizeof(P) * capL); if (!L) { printf("OOM\n"); return 2; } }
                L[nL].h = hW; L[nL].g = g ^ (g >> 1); nL++;
            }
        }
        if (g + 1 == N) break;
        int i = __builtin_ctzll(g + 1);
        int xi = x[i];
        for (int c = 0; c < C; c++) {
            const int *Sc = SS + ((size_t)c * OB + i) * OB;
            ll s = u[(size_t)c * OB + i] - (ll)Sc[i] * xi;
            ll dW = -2LL * xi * s;
            W[c] += dW; hW += coef[c] * (u64)dW;
            ll d = -2LL * xi; ll *uc = u + (size_t)c * OB;
            for (int b = 0; b < OB; b++) uc[b] += (ll)Sc[b] * d;
        }
        rho += -2LL * xi * size[i]; x[i] = -xi;
    }
    if (mode == 2) {
        qsort(L, nL, sizeof(P), cmpP);
        for (size_t i = 0; i < nL; i++) {
            u64 need = hT - L[i].h;
            size_t lo = 0, hi = nL;
            while (lo < hi) { size_t mid = (lo + hi) >> 1; if (L[mid].h < need) lo = mid + 1; else hi = mid; }
            for (size_t k = lo; k < nL && L[k].h == need; k++) { count++; if (printed < 20) { printed++; printf("SOL %llu %llu\n", L[i].g, L[k].g); } }
        }
        printf("LIST %zu\n", nL);
    }
    printf("COUNT %lld\n", count);
    return 0;
}
