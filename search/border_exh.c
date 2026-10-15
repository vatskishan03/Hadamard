/* Exhaustive border search for the character-bordered Szollosi array (search/tbsz.py, solve_border_exh).

   Candidates: n quaternary vectors v_i of length dim = 2q (exponents e, entry i^e), the normalised quaternary
   vectors of W = V_T + V_T. The border needs k = 2t rows u_a and k columns c_b among them, rows pairwise
   orthogonal, columns pairwise orthogonal, and every u_a Mc^* c_b / 2q a fourth root of unity. The caller passes
   P_i = 2 v_i Mc^* (Gaussian integers), so the test is sum_l P_il i^(e_jl) in {M, iM, -M, -iM} with M = 4q.

   Input (little endian): int32 n, dim, k, M, nfirst; n*dim uint8 exponents; n*dim pairs int32 (re, im) of P;
   nfirst int32 indices of the rows tried first (all n when nfirst = 0).
   Or, with separate row and column candidates (n = -1): int32 -1, nr, nc, dim, k, M, nfirst; nr*dim uint8 row
   exponents; nr*dim pairs of P for the rows; nc*dim uint8 column exponents; nfirst int32 first rows. Then the search
   is exhaustive over those candidates (a restriction of the full problem, e.g. to a fixed first row).
   Output: "SOLUTION rows r_1 .. r_k cols c_1 .. c_k" (at most maxsols of them), then "DONE ...".
   The search is exhaustive: rows depth-first with the common compatible columns intersected as bitsets, then a
   k-clique of orthogonal columns among those. "solutions 0" proves that no border exists for this Mc.

       cc -O3 -march=native -fopenmp border_exh.c -o border_exh;  border_exh FILE [maxsols]                        */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static int n, nc, dim, k, M, nfirst, W, WC;
static uint64_t *B, *O, *OC;
static long long nodes = 0, sols = 0, maxsols = 1;
static int R[128], Cc[128];

static int popc(const uint64_t *a) {
    int s = 0;
    for (int w = 0; w < W; w++) s += __builtin_popcountll(a[w]);
    return s;
}

static int popcc(const uint64_t *a) {
    int s = 0;
    for (int w = 0; w < WC; w++) s += __builtin_popcountll(a[w]);
    return s;
}

/* a clique of size need - depth in the orthogonality graph among the set bits of S (indices increasing) */
static int col_clique(const uint64_t *S, int depth, int need) {
    if (depth == need) return 1;
    if (popcc(S) < need - depth) return 0;
    uint64_t *T = malloc(sizeof(uint64_t) * WC);
    for (int w = 0; w < WC; w++) {
        uint64_t x = S[w];
        while (x) {
            int b = __builtin_ctzll(x);
            x &= x - 1;
            int j = w * 64 + b;
            const uint64_t *Oj = OC + (size_t)j * WC;
            for (int v = 0; v < WC; v++) T[v] = v < w ? 0 : (S[v] & Oj[v]);
            T[w] &= b == 63 ? 0 : (~0ULL << (b + 1));
            Cc[depth] = j;
            if (col_clique(T, depth + 1, need)) { free(T); return 1; }
        }
    }
    free(T);
    return 0;
}

static void report(void) {
    sols++;
    printf("SOLUTION rows");
    for (int a = 0; a < k; a++) printf(" %d", R[a]);
    printf(" cols");
    for (int a = 0; a < k; a++) printf(" %d", Cc[a]);
    printf("\n");
    fflush(stdout);
}

/* cand: rows still allowed (orthogonal to the chosen ones, later index); com: their common compatible columns */
static int row_dfs(int depth, const uint64_t *cand, const uint64_t *com) {
    nodes++;
    if (depth == k) {
        if (col_clique(com, 0, k)) { report(); return sols >= maxsols; }
        return 0;
    }
    uint64_t *ncd = malloc(sizeof(uint64_t) * W), *ncom = malloc(sizeof(uint64_t) * WC);
    int stop = 0;
    for (int w = 0; w < W && !stop; w++) {
        uint64_t x = cand[w];
        while (x && !stop) {
            int b = __builtin_ctzll(x);
            x &= x - 1;
            int r = w * 64 + b;
            const uint64_t *Br = B + (size_t)r * WC, *Or = O + (size_t)r * W;
            int cnt = 0;
            for (int v = 0; v < WC; v++) { ncom[v] = com[v] & Br[v]; cnt += __builtin_popcountll(ncom[v]); }
            if (cnt < k) continue;
            for (int v = 0; v < W; v++) ncd[v] = v < w ? 0 : (cand[v] & Or[v]);
            ncd[w] &= b == 63 ? 0 : (~0ULL << (b + 1));
            if (popc(ncd) < k - depth - 1) continue;
            R[depth] = r;
            if (row_dfs(depth + 1, ncd, ncom)) stop = 1;
        }
    }
    free(ncd);
    free(ncom);
    return stop;
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: border_exh FILE [maxsols]\n"); return 1; }
    if (argc >= 3) maxsols = atoll(argv[2]);
    FILE *f = fopen(argv[1], "rb");
    if (!f) return 2;
    int first;
    if (fread(&first, 4, 1, f) != 1) return 3;
    int sep = first < 0;
    int hdr[6];
    if (sep) {
        if (fread(hdr, 4, 6, f) != 6) return 3;
        n = hdr[0]; nc = hdr[1]; dim = hdr[2]; k = hdr[3]; M = hdr[4]; nfirst = hdr[5];
    } else {
        if (fread(hdr, 4, 4, f) != 4) return 3;
        n = first; nc = first; dim = hdr[0]; k = hdr[1]; M = hdr[2]; nfirst = hdr[3];
    }
    W = (n + 63) / 64;
    WC = (nc + 63) / 64;
    uint8_t *E = malloc((size_t)n * dim);
    int32_t *P = malloc(sizeof(int32_t) * 2 * (size_t)n * dim);
    if (fread(E, 1, (size_t)n * dim, f) != (size_t)n * dim) return 4;
    if (fread(P, 4, 2 * (size_t)n * dim, f) != 2 * (size_t)n * dim) return 5;
    uint8_t *EC = E;
    if (sep) {
        EC = malloc((size_t)nc * dim);
        if (fread(EC, 1, (size_t)nc * dim, f) != (size_t)nc * dim) return 4;
    }
    int *firsts = malloc(sizeof(int) * (nfirst ? nfirst : n));
    if (nfirst) {
        if (fread(firsts, 4, nfirst, f) != (size_t)nfirst) return 6;
    } else {
        nfirst = n;
        for (int i = 0; i < n; i++) firsts[i] = i;
    }
    fclose(f);
    B = calloc((size_t)n * WC, 8);
    O = calloc((size_t)n * W, 8);
    OC = sep ? calloc((size_t)nc * WC, 8) : O;
    if (!B || !O || !OC) { fprintf(stderr, "out of memory\n"); return 7; }
    clock_t t0 = clock();
    /* i^e as (re, im) */
    static const int RE[4] = {1, 0, -1, 0}, IM[4] = {0, 1, 0, -1};
#pragma omp parallel for schedule(dynamic, 64)
    for (int i = 0; i < n; i++) {
        const int32_t *Pi = P + 2 * (size_t)i * dim;
        const uint8_t *Ei = E + (size_t)i * dim;
        uint64_t *Bi = B + (size_t)i * WC, *Oi = O + (size_t)i * W;
        for (int j = 0; j < nc; j++) {                /* compatibility with the column candidates */
            const uint8_t *Ej = EC + (size_t)j * dim;
            long long sr = 0, si = 0;
            for (int l = 0; l < dim; l++) {
                int e = Ej[l];
                long long pr = Pi[2 * l], pi = Pi[2 * l + 1];
                sr += pr * RE[e] - pi * IM[e];
                si += pr * IM[e] + pi * RE[e];
            }
            if ((sr == 0 && (si == M || si == -M)) || (si == 0 && (sr == M || sr == -M)))
                Bi[j >> 6] |= 1ULL << (j & 63);
        }
        for (int j = 0; j < n; j++) {                 /* orthogonality of rows: exponent differences */
            const uint8_t *Ej = E + (size_t)j * dim;
            int cr = 0, ci = 0;
            for (int l = 0; l < dim; l++) {
                int d = (Ei[l] - Ej[l]) & 3;
                cr += RE[d];
                ci += IM[d];
            }
            if (cr == 0 && ci == 0) Oi[j >> 6] |= 1ULL << (j & 63);
        }
    }
    if (sep) {
#pragma omp parallel for schedule(dynamic, 64)
        for (int i = 0; i < nc; i++) {
            const uint8_t *Ei = EC + (size_t)i * dim;
            uint64_t *Oi = OC + (size_t)i * WC;
            for (int j = 0; j < nc; j++) {
                const uint8_t *Ej = EC + (size_t)j * dim;
                int cr = 0, ci = 0;
                for (int l = 0; l < dim; l++) {
                    int d = (Ei[l] - Ej[l]) & 3;
                    cr += RE[d];
                    ci += IM[d];
                }
                if (cr == 0 && ci == 0) Oi[j >> 6] |= 1ULL << (j & 63);
            }
        }
    }
    double tprep = (double)(clock() - t0) / CLOCKS_PER_SEC;
    uint64_t *cand = calloc(W, 8), *com = calloc(WC, 8);
    int stop = 0;
    for (int fi = 0; fi < nfirst && !stop; fi++) {
        int r = firsts[fi];
        const uint64_t *Br = B + (size_t)r * WC, *Or = O + (size_t)r * W;
        int cnt = 0;
        for (int v = 0; v < WC; v++) { com[v] = Br[v]; cnt += __builtin_popcountll(com[v]); }
        if (cnt < k) continue;
        /* a first row from the list, the others over all candidates with increasing index */
        memcpy(cand, Or, sizeof(uint64_t) * W);
        R[0] = r;
        if (row_dfs(1, cand, com)) stop = 1;
    }
    printf("DONE n %d k %d firsts %d nodes %lld solutions %lld prep %.1f seconds %.1f\n", n, k, nfirst, nodes, sols,
           tprep, (double)(clock() - t0) / CLOCKS_PER_SEC);
    return 0;
}
