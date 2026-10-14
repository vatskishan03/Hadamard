/* Exhaustive search for Legendre pairs (M, N) over an abelian group G = Z_m1 x ... x Z_mr of odd order v with
 * M of skew type (m[0] = +1, m[-g] = -m[g]) and N symmetric (n[-g] = n[g]), row sums 1 and 1:
 *     PAF_m(s) + PAF_n(s) = -2   for every s != 0   (enough: one s from each pair {s, -s}).
 * Such a pair gives a Szekeres-type skew-Hadamard matrix of order 2(v+1) with an involutive anti-automorphism,
 * hence a symmetric Hadamard matrix of order (2v+2)(2v+1) (attack/sym_mm1.py szekeres).
 *
 *   szek_group m1,m2,... [maxsol]      (e.g. szek_group 3,15)
 *
 * Output: "SOL m n" with +/- strings indexed by the mixed-radix integer g = (...(g1 m2 + g2) m3 + ...).
 * N is enumerated completely and hashed on its PAF vector; every skew-type m with m[rep_1] = +1 (m -> m^T) that
 * passes the PSD filter |m^(chi)|^2 <= 2v + 2 is looked up.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <complex.h>

static int v, h, r, mod[8];
static int *addt, *negt, rep[512];
static int16_t *NP;
static uint64_t *NM;
static int64_t nN, *head, *nxt;
static uint64_t hmask;

static void coords(int g, int *c) { for (int i = r - 1; i >= 0; i--) { c[i] = g % mod[i]; g /= mod[i]; } }
static int index_of(const int *c) { int g = 0; for (int i = 0; i < r; i++) g = g * mod[i] + c[i]; return g; }

static uint64_t hashv(const int16_t *p) {
    uint64_t x = 1469598103934665603ULL;
    for (int i = 0; i < h; i++) { x ^= (uint16_t)p[i]; x *= 1099511628211ULL; }
    return x;
}

static void paf(const int8_t *a, int16_t *out) {
    for (int i = 0; i < h; i++) {
        int s = rep[i], t = 0;
        for (int g = 0; g < v; g++) t += a[g] * a[addt[g * v + s]];
        out[i] = (int16_t)t;
    }
}

int main(int argc, char **argv) {
    char *p = argv[1];
    r = 0; v = 1;
    while (*p) { mod[r] = (int)strtol(p, &p, 10); v *= mod[r++]; if (*p == ',') p++; }
    long maxsol = argc > 2 ? atol(argv[2]) : 1000000;
    if (v % 2 == 0) { fprintf(stderr, "v must be odd\n"); return 1; }
    h = (v - 1) / 2;
    addt = malloc(sizeof(int) * v * v); negt = malloc(sizeof(int) * v);
    for (int g = 0; g < v; g++) {
        int a[8], b[8], c[8];
        coords(g, a);
        for (int i = 0; i < r; i++) c[i] = (mod[i] - a[i]) % mod[i];
        negt[g] = index_of(c);
        for (int k = 0; k < v; k++) { coords(k, b); for (int i = 0; i < r; i++) c[i] = (a[i] + b[i]) % mod[i]; addt[g * v + k] = index_of(c); }
    }
    { int seen[4096] = {0}, k = 0; for (int g = 1; g < v; g++) if (!seen[g]) { seen[g] = seen[negt[g]] = 1; rep[k++] = g; } }
    int n0 = (h % 2 == 0) ? 1 : -1, k = (h % 2 == 0) ? h / 2 : (h - 1) / 2;
    double cnt = 1; for (int i = 0; i < k; i++) cnt = cnt * (h - i) / (i + 1);
    nN = (int64_t)(cnt + 0.5);
    fprintf(stderr, "G = Z_%d", mod[0]); for (int i = 1; i < r; i++) fprintf(stderr, " x Z_%d", mod[i]);
    fprintf(stderr, ", v = %d: N has n[0] = %d and %d of %d pairs, %lld sets\n", v, n0, k, h, (long long)nN);
    NP = malloc(sizeof(int16_t) * h * nN); NM = malloc(sizeof(uint64_t) * nN);
    uint64_t hs = 1; while (hs < 2 * (uint64_t)nN) hs <<= 1; hmask = hs - 1;
    head = malloc(sizeof(int64_t) * hs); nxt = malloc(sizeof(int64_t) * nN);
    for (uint64_t i = 0; i < hs; i++) head[i] = -1;
    int8_t a[4096];
    int idx[64];
    for (int i = 0; i < k; i++) idx[i] = i;
    int64_t c = 0;
    while (1) {
        for (int g = 0; g < v; g++) a[g] = 1;
        a[0] = n0;
        uint64_t mask = 0;
        for (int i = 0; i < k; i++) { int g = rep[idx[i]]; a[g] = a[negt[g]] = -1; mask |= 1ULL << idx[i]; }
        paf(a, NP + c * h); NM[c] = mask;
        uint64_t hv = hashv(NP + c * h) & hmask;
        nxt[c] = head[hv]; head[hv] = c; c++;
        int i = k - 1;
        while (i >= 0 && idx[i] == h - k + i) i--;
        if (i < 0) break;
        idx[i]++;
        for (int j = i + 1; j < k; j++) idx[j] = idx[j - 1] + 1;
    }
    /* characters: chi_t(g) = exp(2 pi i sum_j t_j g_j / m_j); sine table S[t][i] = sin(theta_t(rep_i)) */
    double *S = malloc(sizeof(double) * v * h);
    for (int t = 0; t < v; t++) {
        int ct[8], cg[8]; coords(t, ct);
        for (int i = 0; i < h; i++) {
            coords(rep[i], cg);
            double th = 0; for (int j = 0; j < r; j++) th += (double)ct[j] * cg[j] / mod[j];
            S[t * h + i] = sin(2 * M_PI * th);
        }
    }
    long sols = 0;
    int16_t T[2048], P[2048];
    uint64_t total = 1ULL << (h - 1), passed = 0;
    for (uint64_t bits = 0; bits < total; bits++) {
        a[0] = 1;
        for (int i = 0; i < h; i++) {
            int s = i == 0 ? 1 : (((bits >> (i - 1)) & 1) ? -1 : 1);
            a[rep[i]] = s; a[negt[rep[i]]] = -s;
        }
        int bad = 0;
        for (int t = 1; t < v && !bad; t++) {
            double x = 0;
            for (int i = 0; i < h; i++) x += a[rep[i]] * S[t * h + i];
            if (1 + 4 * x * x > 2 * v + 2 + 1e-6) bad = 1;
        }
        if (bad) continue;
        passed++;
        paf(a, P);
        for (int s = 0; s < h; s++) T[s] = (int16_t)(-2 - P[s]);
        uint64_t hv = hashv(T) & hmask;
        for (int64_t j = head[hv]; j >= 0; j = nxt[j]) {
            if (memcmp(NP + j * h, T, sizeof(int16_t) * h)) continue;
            printf("SOL ");
            for (int g = 0; g < v; g++) putchar(a[g] > 0 ? '+' : '-');
            printf(" ");
            int8_t b[4096];
            for (int g = 0; g < v; g++) b[g] = 1;
            b[0] = n0;
            for (int i = 0; i < h; i++) if (NM[j] >> i & 1) b[rep[i]] = b[negt[rep[i]]] = -1;
            for (int g = 0; g < v; g++) putchar(b[g] > 0 ? '+' : '-');
            printf("\n");
            fflush(stdout);
            if (++sols >= maxsol) { fprintf(stderr, "stopped at %ld solutions\n", sols); return 0; }
        }
    }
    fprintf(stderr, "done: %llu skew sequences, %llu pass the PSD filter, %ld solutions\n",
            (unsigned long long)total, (unsigned long long)passed, sols);
    printf("DONE %ld\n", sols);
    return 0;
}
