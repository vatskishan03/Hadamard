/* Exhaustive search for Legendre pairs (M, N) of odd length v in Z_v with M of skew type and N symmetric
 * ("Szekeres pairs"): m[0] = +1, m[-g] = -m[g] for g != 0; n[-g] = n[g]; row sums 1 and 1;
 *     PAF_m(s) + PAF_n(s) = -2   for s = 1 .. (v-1)/2.
 * Such a pair gives a Szekeres-type skew-Hadamard matrix of order 2(v+1) with an involutive anti-automorphism,
 * hence a symmetric Hadamard matrix of order (2v+2)(2v+1) (search/sym_mm1.py).
 *
 *   szek_pair v [maxsol]      (prints "SOL m n" lines with +/- strings; M up to m -> m^T: s_1 = +1)
 *
 * N is enumerated completely (row sum 1 fixes the number of pairs {g, -g} in N) and hashed on its PAF vector;
 * every skew-type m with m[1] = +1 is looked up.  PSD filter on m: |m^(k)|^2 <= 2v + 2 for k != 0.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

static int v, h;
static int16_t *NP;          /* PAF vectors of N, h entries each */
static uint32_t *NM;         /* N as bitmask over pairs (bit g-1 set: n[g] = n[-g] = -1) */
static int64_t nN;
static int64_t *head, *nxt;
static uint64_t hmask;

static uint64_t hashv(const int16_t *p) {
    uint64_t x = 1469598103934665603ULL;
    for (int i = 0; i < h; i++) { x ^= (uint16_t)p[i]; x *= 1099511628211ULL; }
    return x;
}

static void paf(const int8_t *a, int16_t *out) {
    for (int s = 1; s <= h; s++) {
        int t = 0;
        for (int g = 0; g < v; g++) { int k = g + s; if (k >= v) k -= v; t += a[g] * a[k]; }
        out[s - 1] = (int16_t)t;
    }
}

int main(int argc, char **argv) {
    v = atoi(argv[1]);
    long maxsol = argc > 2 ? atol(argv[2]) : 1000000;
    h = (v - 1) / 2;
    /* n[0] and number of pairs k: 1 = n0 + 2 (h - 2k)  ->  n0 = 1 requires h = 2k, n0 = -1 requires h - 2k = 1 */
    int n0 = (h % 2 == 0) ? 1 : -1, k = (h % 2 == 0) ? h / 2 : (h - 1) / 2;
    /* count C(h, k) */
    double cnt = 1; for (int i = 0; i < k; i++) cnt = cnt * (h - i) / (i + 1);
    nN = (int64_t)(cnt + 0.5);
    fprintf(stderr, "v = %d, h = %d: N has n[0] = %d and %d of %d pairs, %lld sets\n", v, h, n0, k, h, (long long)nN);
    NP = malloc(sizeof(int16_t) * h * nN); NM = malloc(sizeof(uint32_t) * nN);
    uint64_t hs = 1; while (hs < 2 * (uint64_t)nN) hs <<= 1; hmask = hs - 1;
    head = malloc(sizeof(int64_t) * hs); nxt = malloc(sizeof(int64_t) * nN);
    for (uint64_t i = 0; i < hs; i++) head[i] = -1;
    int8_t a[512];
    int idx[64];
    for (int i = 0; i < k; i++) idx[i] = i;
    int64_t c = 0;
    while (1) {
        for (int g = 0; g < v; g++) a[g] = 1;
        a[0] = n0;
        uint32_t mask = 0;
        for (int i = 0; i < k; i++) { int g = idx[i] + 1; a[g] = a[v - g] = -1; mask |= 1u << idx[i]; }
        paf(a, NP + c * h); NM[c] = mask;
        uint64_t hv = hashv(NP + c * h) & hmask;
        nxt[c] = head[hv]; head[hv] = c; c++;
        int i = k - 1;
        while (i >= 0 && idx[i] == h - k + i) i--;
        if (i < 0) break;
        idx[i]++;
        for (int j = i + 1; j < k; j++) idx[j] = idx[j - 1] + 1;
    }
    if (c != nN) { fprintf(stderr, "count mismatch %lld %lld\n", (long long)c, (long long)nN); return 1; }
    /* cosine/sine tables for the PSD filter */
    double *sn = malloc(sizeof(double) * v * v);
    for (int kk = 0; kk < v; kk++) for (int g = 0; g < v; g++) sn[kk * v + g] = sin(2 * M_PI * (double)kk * g / v);
    long sols = 0;
    int16_t T[256], P[256];
    uint64_t total = 1ULL << (h - 1), passed = 0;
    for (uint64_t bits = 0; bits < total; bits++) {
        a[0] = 1;
        a[1] = 1; a[v - 1] = -1;
        for (int g = 2; g <= h; g++) { int s = (bits >> (g - 2)) & 1 ? -1 : 1; a[g] = s; a[v - g] = -s; }
        /* PSD: m^(k) = 1 + i * 2 sum_{g=1..h} m[g] sin(2 pi g k / v) */
        int bad = 0;
        for (int kk = 1; kk <= h && !bad; kk++) {
            double t = 0;
            for (int g = 1; g <= h; g++) t += a[g] * sn[kk * v + g];
            if (1 + 4 * t * t > 2 * v + 2 + 1e-6) bad = 1;
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
            int8_t b[512];
            for (int g = 0; g < v; g++) b[g] = 1;
            b[0] = n0;
            for (int i = 0; i < h; i++) if (NM[j] >> i & 1) b[i + 1] = b[v - i - 1] = -1;
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
