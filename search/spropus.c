/* Skew-propus quadruples of odd length n invariant under a multiplier group H (|H| odd, -1 not in H):
 *   A3 symmetric (constant on the orbits of <H, -1>), A1 constant on the H-orbits, A0 of skew type
 *   (a0[0] = 1, a0[-g] = -a0[g], constant on the H-orbits), with PAF_A3 + 2 PAF_A1 + PAF_A0 = 4n delta.
 * Such a quadruple gives a symmetric OD(4n; 1, 4n - 1) in the Balonin array and so a symmetric Hadamard matrix of
 * order 4n (4n - 1) (search/sym_mm1.py, search/NOTES.md section 31).
 *
 *   cc -O3 -march=native -o spropus spropus.c
 *   spropus n g1[,g2,...] [maxsol]
 *
 * PAFs are computed on the orbit level from N(o, o'; s) = #{g in o : g + s in o'} at one shift s per class of
 * <H, -1>. Every A1 with row sum a1 > 0 (a3^2 + 2 a1^2 = 4n - 1, a3 > 0) is hashed on its class vector; every pair
 * (A3, A0) looks up the needed A1. Hits are expanded and checked on Z_n before they are printed.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

static int n, O, P, *orb, *osz, *rep, *negof, cls[512], ncls;
static int *Nt;                       /* Nt[(c * O + o) * O + o2] */
static int pairs[512], npairs, zero;

static void classvec(const int8_t *x, int *out) {   /* x: sign per orbit */
    for (int c = 0; c < ncls; c++) {
        int t = 0;
        const int *T = Nt + (size_t)c * O * O;
        for (int o = 0; o < O; o++) {
            int r = 0;
            for (int o2 = 0; o2 < O; o2++) r += T[o * O + o2] * x[o2];
            t += x[o] * r;
        }
        out[c] = t;
    }
}

static uint64_t hashv(const int *p) {
    uint64_t h = 1469598103934665603ULL;
    for (int i = 0; i < ncls; i++) { h ^= (uint32_t)(p[i] + 100000); h *= 1099511628211ULL; }
    return h;
}

static void expand(const int8_t *x, int8_t *s) { for (int g = 0; g < n; g++) s[g] = x[orb[g]]; }

static int check(const int8_t *a3, const int8_t *a1, const int8_t *a0) {
    for (int s = 1; s < n; s++) {
        int t = 0;
        for (int g = 0; g < n; g++) { int k = (g + s) % n; t += a3[g] * a3[k] + 2 * a1[g] * a1[k] + a0[g] * a0[k]; }
        if (t) return 0;
    }
    return 1;
}

int main(int argc, char **argv) {
    n = atoi(argv[1]);
    int gens[16], ng = 0;
    for (char *p = argv[2]; *p; ) { gens[ng++] = (int)strtol(p, &p, 10); if (*p == ',') p++; }
    long maxsol = argc > 3 ? atol(argv[3]) : 1;
    /* H */
    char *inH = calloc(n, 1); inH[1] = 1;
    for (int changed = 1; changed; ) {
        changed = 0;
        for (int h = 0; h < n; h++) if (inH[h]) for (int i = 0; i < ng; i++) { int k = (int)((long)h * gens[i] % n); if (!inH[k]) { inH[k] = 1; changed = 1; } }
    }
    if (inH[n - 1]) { fprintf(stderr, "-1 is in H\n"); return 1; }
    int Hn = 0, Hl[4096]; for (int h = 0; h < n; h++) if (inH[h]) Hl[Hn++] = h;
    orb = malloc(sizeof(int) * n); for (int g = 0; g < n; g++) orb[g] = -1;
    rep = malloc(sizeof(int) * n); osz = calloc(n, sizeof(int));
    O = 0;
    for (int g = 0; g < n; g++) if (orb[g] < 0) { for (int i = 0; i < Hn; i++) orb[(long)g * Hl[i] % n] = O; rep[O++] = g; }
    for (int g = 0; g < n; g++) osz[orb[g]]++;
    negof = malloc(sizeof(int) * O);
    for (int o = 0; o < O; o++) negof[o] = orb[(n - rep[o]) % n];
    zero = orb[0];
    npairs = 0; for (int o = 0; o < O; o++) if (o != zero && o < negof[o]) pairs[npairs++] = o;
    /* classes of <H, -1> on Z_n \ 0: one representative shift per pair of orbits */
    ncls = 0; for (int i = 0; i < npairs; i++) cls[ncls++] = rep[pairs[i]];
    Nt = calloc((size_t)ncls * O * O, sizeof(int));
    for (int c = 0; c < ncls; c++) for (int g = 0; g < n; g++) Nt[((size_t)c * O + orb[g]) * O + orb[(g + cls[c]) % n]]++;
    int reps[64][2], nr = 0;
    for (int a3 = 1; a3 < 2 * n; a3 += 2) for (int a1 = 1; a1 < 2 * n; a1 += 2) if (a3 * a3 + 2 * a1 * a1 == 4 * n - 1) { reps[nr][0] = a3; reps[nr][1] = a1; nr++; }
    fprintf(stderr, "n = %d, |H| = %d: %d orbits, %d orbit pairs, %d classes; row sums", n, Hn, O, npairs, ncls);
    for (int i = 0; i < nr; i++) fprintf(stderr, " (%d, %d)", reps[i][0], reps[i][1]);
    fprintf(stderr, "\n");
    if (O > 34) { fprintf(stderr, "too many orbits\n"); return 1; }
    long sols = 0;
    int8_t x[64], s3[4096], s1[4096], s0[4096];
    int cv[512], cw[512];
    for (int r = 0; r < nr && sols < maxsol; r++) {
        int a3 = reps[r][0], a1 = reps[r][1];
        /* A1 table */
        uint64_t cap = 1 << 20, cnt = 0;
        int *vals = malloc(sizeof(int) * ncls * cap); uint64_t *codes = malloc(sizeof(uint64_t) * cap);
        for (uint64_t b = 0; b < (1ULL << O); b++) {
            int sum = 0;
            for (int o = 0; o < O; o++) { x[o] = (b >> o & 1) ? -1 : 1; sum += x[o] * osz[o]; }
            if (sum != a1) continue;
            if (cnt == cap) { cap *= 2; vals = realloc(vals, sizeof(int) * ncls * cap); codes = realloc(codes, sizeof(uint64_t) * cap); }
            classvec(x, vals + cnt * ncls); codes[cnt++] = b;
        }
        uint64_t hs = 1; while (hs < 2 * cnt + 2) hs <<= 1;
        int64_t *head = malloc(sizeof(int64_t) * hs), *nxt = malloc(sizeof(int64_t) * (cnt + 1));
        for (uint64_t i = 0; i < hs; i++) head[i] = -1;
        for (uint64_t i = 0; i < cnt; i++) { uint64_t h = hashv(vals + i * ncls) & (hs - 1); nxt[i] = head[h]; head[h] = (int64_t)i; }
        fprintf(stderr, "  (a3, a1) = (%d, %d): %llu A1\n", a3, a1, (unsigned long long)cnt);
        /* A3 x A0 */
        int8_t y3[64], y0[64];
        for (uint64_t b3 = 0; b3 < (1ULL << (npairs + 1)) && sols < maxsol; b3++) {
            y3[zero] = (b3 & 1) ? -1 : 1;
            int sum = y3[zero];
            for (int i = 0; i < npairs; i++) { int s = (b3 >> (i + 1) & 1) ? -1 : 1; y3[pairs[i]] = y3[negof[pairs[i]]] = s; sum += 2 * s * osz[pairs[i]]; }
            if (sum != a3) continue;
            int c3[512]; classvec(y3, c3);
            for (uint64_t b0 = 0; b0 < (1ULL << npairs) && sols < maxsol; b0++) {
                y0[zero] = 1;
                for (int i = 0; i < npairs; i++) { int s = (b0 >> i & 1) ? -1 : 1; y0[pairs[i]] = s; y0[negof[pairs[i]]] = -s; }
                classvec(y0, cv);
                int ok = 1;
                for (int c = 0; c < ncls; c++) { int t = -(c3[c] + cv[c]); if (t & 1) { ok = 0; break; } cw[c] = t / 2; }
                if (!ok) continue;
                uint64_t h = hashv(cw) & (hs - 1);
                for (int64_t j = head[h]; j >= 0; j = nxt[j]) {
                    if (memcmp(vals + j * ncls, cw, sizeof(int) * ncls)) continue;
                    for (int o = 0; o < O; o++) x[o] = (codes[j] >> o & 1) ? -1 : 1;
                    expand(y3, s3); expand(x, s1); expand(y0, s0);
                    if (!check(s3, s1, s0)) { fprintf(stderr, "WRONG\n"); return 2; }
                    printf("SOL n=%d A3 ", n); for (int g = 0; g < n; g++) putchar(s3[g] > 0 ? '+' : '-');
                    printf(" A1 "); for (int g = 0; g < n; g++) putchar(s1[g] > 0 ? '+' : '-');
                    printf(" A0 "); for (int g = 0; g < n; g++) putchar(s0[g] > 0 ? '+' : '-');
                    printf("\n"); fflush(stdout);
                    if (++sols >= maxsol) break;
                }
            }
        }
        free(vals); free(codes); free(head); free(nxt);
    }
    printf("DONE n=%d |H|=%d solutions %ld%s\n", n, Hn, sols, sols >= maxsol ? " (stopped)" : "");
    return 0;
}
