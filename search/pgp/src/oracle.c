/* oracle.c -- independent brute-force reference for pgp.c (shares no tables with it; only the labelling of the orbits:
 * orbit j of Z_w^* under H is { g^(j + m t) }, g the least primitive root).
 *
 *   oracle rows w h OUT      all rows (X, Y): H-invariant ternary X, Y on Z_w, PAF_X + PAF_Y = 0 off 0, weight_X + weight_Y = w,
 *                            sums > 0 (sign-normalised); PAF from the definition on Z_w.  Prints "X Y" as strings over +-0,
 *                            first the infinity point (the point 0 of Z_w) then the orbits 0..m-1.  Count on stdout.
 *   oracle pgp w h OUT       all pairs (a, b) of H-invariant +-1 sequences on Z_{2w} (CRT: index y <-> (y mod 2, y mod w)) with
 *                            PAF_a + PAF_b = 0 off 0, PAF from the definition on the cyclic group Z_{2w}.  Prints
 *                            "a b" with a as the layer-0 string (infinity, orbits) then the layer-1 string.  Count on stdout.
 *
 *   cc -O3 -o oracle oracle.c
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

static int is_prime(int n) { if (n < 2) return 0; for (int d = 2; d * d <= n; d++) if (n % d == 0) return 0; return 1; }
static int primroot(int p) {
    for (int a = 2; a < p; a++) {
        int ok = 1;
        for (int q = 2; q < p; q++) if ((p - 1) % q == 0 && is_prime(q)) { long x = 1; for (int i = 0; i < (p - 1) / q; i++) x = x * a % p; if (x == 1) { ok = 0; break; } }
        if (ok) return a;
    }
    return -1;
}

typedef struct { uint64_t k0, k1, k2, k3; uint32_t id; } rec_t;     /* up to 32 bytes of key */

typedef struct { uint64_t h; uint32_t id; } he;
static int cmp_he(const void *a, const void *b) { uint64_t x = ((const he *)a)->h, y = ((const he *)b)->h; return (x > y) - (x < y); }

int main(int argc, char **argv) {
    if (argc < 5) { fprintf(stderr, "usage: oracle rows|pgp w h OUT\n"); return 2; }
    int mode = !strcmp(argv[1], "rows") ? 0 : 1;
    int w = atoi(argv[2]), h = atoi(argv[3]);
    FILE *out = fopen(argv[4], "w");
    int g = primroot(w), m = (w - 1) / h;
    int *dlog = calloc(w, sizeof(int)); long x = 1;
    for (int i = 0; i < w - 1; i++) { dlog[x] = i; x = x * g % w; }
    int n1 = mode == 0 ? m + 1 : 2 * (m + 1);              /* number of free symbols per sequence */
    long long total = 1; for (int i = 0; i < n1; i++) total *= (mode == 0 ? 3 : 2);
    int len = mode == 0 ? w : 2 * w;
    int nshift = mode == 0 ? w - 1 : 2 * w - 1;
    /* list all sequences with their PAF vector over shifts 1..nshift (full vector, no symmetry used) */
    int8_t *paf = malloc((size_t)total * nshift);
    long long *wt = malloc(total * sizeof(long long)), *sm = malloc(total * sizeof(long long));
    int *seq = malloc(sizeof(int) * len);
    for (long long id = 0; id < total; id++) {
        int sym[64]; long long t = id;
        for (int i = 0; i < n1; i++) { if (mode == 0) { sym[i] = (int)(t % 3) - 1; t /= 3; } else { sym[i] = 2 * (int)(t % 2) - 1; t /= 2; } }
        long long weight = 0, sum = 0;
        if (mode == 0) {
            seq[0] = sym[0];
            for (int y = 1; y < w; y++) seq[y] = sym[1 + dlog[y] % m];
        } else {
            for (int y = 0; y < 2 * w; y++) {
                int e = y % 2, r = y % w;                    /* CRT: y <-> (e, r) */
                seq[y] = r == 0 ? sym[e * (m + 1)] : sym[e * (m + 1) + 1 + dlog[r] % m];
            }
        }
        for (int y = 0; y < len; y++) { weight += seq[y] * seq[y]; sum += seq[y]; }
        wt[id] = weight; sm[id] = sum;
        for (int s = 1; s <= nshift; s++) { int p = 0; for (int y = 0; y < len; y++) p += seq[y] * seq[(y + s) % len]; paf[(size_t)id * nshift + s - 1] = (int8_t)p; }
    }
    /* hash join on the full PAF vector */
    /* key = hash of negated vector; use sorting of (hash, id) with full comparison afterwards */
    he *H = malloc(sizeof(he) * total);
    for (long long id = 0; id < total; id++) {
        uint64_t hh = 1469598103934665603ull;
        for (int s = 0; s < nshift; s++) { hh ^= (uint8_t)paf[(size_t)id * nshift + s]; hh *= 1099511628211ull; }
        H[id].h = hh; H[id].id = (uint32_t)id;
    }
    long long count = 0;
    /* sort by hash */
    {
        qsort(H, total, sizeof(he), cmp_he);
        for (long long ia = 0; ia < total; ia++) {
            long long a = ia;                                       /* sequence a */
            if (mode == 0 && sm[a] <= 0) continue;                  /* row normalisation: both sums positive */
            uint64_t hh = 1469598103934665603ull;
            for (int s = 0; s < nshift; s++) { hh ^= (uint8_t)(-paf[(size_t)a * nshift + s]); hh *= 1099511628211ull; }
            long long lo = 0, hi = total;
            while (lo < hi) { long long mid = (lo + hi) / 2; if (H[mid].h < hh) lo = mid + 1; else hi = mid; }
            for (long long j = lo; j < total && H[j].h == hh; j++) {
                long long b = H[j].id; int ok = 1;
                for (int s = 0; s < nshift; s++) if (paf[(size_t)a * nshift + s] + paf[(size_t)b * nshift + s] != 0) { ok = 0; break; }
                if (!ok) continue;
                if (mode == 0) { if (sm[b] <= 0 || wt[a] + wt[b] != w) continue; }
                else { /* PAF(0) = 2*len automatically for +-1 sequences: a, b each have weight len */ }
                count++;
                char sa[200], sb[200]; int la = 0, lb = 0;
                for (int pass = 0; pass < 2; pass++) {
                    long long id = pass ? b : a; char *dst = pass ? sb : sa; int *l = pass ? &lb : &la; long long t = id;
                    for (int i = 0; i < n1; i++) {
                        int v; if (mode == 0) { v = (int)(t % 3) - 1; t /= 3; } else { v = 2 * (int)(t % 2) - 1; t /= 2; }
                        if (mode == 1 && i == m + 1) dst[(*l)++] = ' ';
                        dst[(*l)++] = v > 0 ? '+' : (v < 0 ? '-' : '0');
                    }
                    dst[*l] = 0;
                }
                fprintf(out, "%s %s\n", sa, sb);
            }
        }
    }
    printf("%lld\n", count);
    fclose(out);
    return 0;
}
