/* Exhaustive border search for the character-bordered Szollosi array (search/tbsz.py, solve_border_exh).

   Candidates: n quaternary vectors v_i of length dim = 2q (exponents e, entry i^e), the normalised quaternary
   vectors of W = V_T + V_T. The border needs k = 2t rows u_a and k columns c_b among them, rows pairwise
   orthogonal, columns pairwise orthogonal, and every u_a Mc^* c_b / 2q a fourth root of unity. The caller passes
   P_i = 2 v_i Mc^* (Gaussian integers), so the test is sum_l P_il i^(e_jl) in {M, iM, -M, -iM} with M = 4q.

   Input (little endian): int32 n, dim, k, M, nfirst; n*dim uint8 exponents; n*dim pairs int32 (re, im) of P;
   nfirst int32 indices of the rows tried first (all n when nfirst = 0).
   Output: "SOLUTION rows r_1 .. r_k cols c_1 .. c_k" (at most maxsols of them), then "DONE ...".
   The search is exhaustive: rows depth-first with the common compatible columns intersected as bitsets, then a
   k-clique of orthogonal columns among those. "solutions 0" proves that no border exists for this Mc.

       cc -O3 -march=native -fopenmp border_exh.c -o border_exh;  border_exh FILE [maxsols]                        */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static int n, dim, k, M, nfirst, W;
static uint64_t *B, *O;
static long long nodes = 0, sols = 0, maxsols = 1;
static int R[128], Cc[128];

static int popc(const uint64_t *a) {
    int s = 0;
    for (int w = 0; w < W; w++) s += __builtin_popcountll(a[w]);
    return s;
}

/* a clique of size need - depth in the orthogonality graph among the set bits of S (indices increasing) */
static int col_clique(const uint64_t *S, int depth, int need) {
    if (depth == need) return 1;
    if (popc(S) < need - depth) return 0;
    uint64_t *T = malloc(sizeof(uint64_t) * W);
    for (int w = 0; w < W; w++) {
        uint64_t x = S[w];
        while (x) {
            int b = __builtin_ctzll(x);
            x &= x - 1;
            int j = w * 64 + b;
            const uint64_t *Oj = O + (size_t)j * W;
            for (int v = 0; v < W; v++) T[v] = v < w ? 0 : (S[v] & Oj[v]);
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
    uint64_t *nc = malloc(sizeof(uint64_t) * W), *ncom = malloc(sizeof(uint64_t) * W);
    int stop = 0;
    for (int w = 0; w < W && !stop; w++) {
        uint64_t x = cand[w];
        while (x && !stop) {
            int b = __builtin_ctzll(x);
            x &= x - 1;
            int r = w * 64 + b;
            const uint64_t *Br = B + (size_t)r * W, *Or = O + (size_t)r * W;
            int cnt = 0;
            for (int v = 0; v < W; v++) { ncom[v] = com[v] & Br[v]; cnt += __builtin_popcountll(ncom[v]); }
            if (cnt < k) continue;
            for (int v = 0; v < W; v++) nc[v] = v < w ? 0 : (cand[v] & Or[v]);
            nc[w] &= b == 63 ? 0 : (~0ULL << (b + 1));
            if (popc(nc) < k - depth - 1) continue;
            R[depth] = r;
            if (row_dfs(depth + 1, nc, ncom)) stop = 1;
        }
    }
    free(nc);
    free(ncom);
    return stop;
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: border_exh FILE [maxsols]\n"); return 1; }
    if (argc >= 3) maxsols = atoll(argv[2]);
    FILE *f = fopen(argv[1], "rb");
    if (!f) return 2;
    int hdr[5];
    if (fread(hdr, 4, 5, f) != 5) return 3;
    n = hdr[0]; dim = hdr[1]; k = hdr[2]; M = hdr[3]; nfirst = hdr[4];
    W = (n + 63) / 64;
    uint8_t *E = malloc((size_t)n * dim);
    int32_t *P = malloc(sizeof(int32_t) * 2 * (size_t)n * dim);
    if (fread(E, 1, (size_t)n * dim, f) != (size_t)n * dim) return 4;
    if (fread(P, 4, 2 * (size_t)n * dim, f) != 2 * (size_t)n * dim) return 5;
    int *firsts = malloc(sizeof(int) * (nfirst ? nfirst : n));
    if (nfirst) {
        if (fread(firsts, 4, nfirst, f) != (size_t)nfirst) return 6;
    } else {
        nfirst = n;
        for (int i = 0; i < n; i++) firsts[i] = i;
    }
    fclose(f);
    B = calloc((size_t)n * W, 8);
    O = calloc((size_t)n * W, 8);
    if (!B || !O) { fprintf(stderr, "out of memory\n"); return 7; }
    clock_t t0 = clock();
    /* i^e as (re, im) */
    static const int RE[4] = {1, 0, -1, 0}, IM[4] = {0, 1, 0, -1};
#pragma omp parallel for schedule(dynamic, 64)
    for (int i = 0; i < n; i++) {
        const int32_t *Pi = P + 2 * (size_t)i * dim;
        const uint8_t *Ei = E + (size_t)i * dim;
        uint64_t *Bi = B + (size_t)i * W, *Oi = O + (size_t)i * W;
        for (int j = 0; j < n; j++) {
            const uint8_t *Ej = E + (size_t)j * dim;
            long long sr = 0, si = 0;
            int cr = 0, ci = 0;                       /* sum_l v_il conj(v_jl): exponent differences */
            for (int l = 0; l < dim; l++) {
                int e = Ej[l];
                long long pr = Pi[2 * l], pi = Pi[2 * l + 1];
                sr += pr * RE[e] - pi * IM[e];
                si += pr * IM[e] + pi * RE[e];
                int d = (Ei[l] - Ej[l]) & 3;
                cr += RE[d];
                ci += IM[d];
            }
            if ((sr == 0 && (si == M || si == -M)) || (si == 0 && (sr == M || sr == -M)))
                Bi[j >> 6] |= 1ULL << (j & 63);
            if (cr == 0 && ci == 0) Oi[j >> 6] |= 1ULL << (j & 63);
        }
    }
    double tprep = (double)(clock() - t0) / CLOCKS_PER_SEC;
    uint64_t *nc = calloc(W, 8), *com = calloc(W, 8);
    int stop = 0;
    for (int fi = 0; fi < nfirst && !stop; fi++) {
        int r = firsts[fi];
        const uint64_t *Br = B + (size_t)r * W, *Or = O + (size_t)r * W;
        int cnt = 0;
        for (int v = 0; v < W; v++) { com[v] = Br[v]; cnt += __builtin_popcountll(com[v]); }
        if (cnt < k) continue;
        /* a first row from the list, the others over all candidates with increasing index */
        memcpy(nc, Or, sizeof(uint64_t) * W);
        R[0] = r;
        if (row_dfs(1, nc, com)) stop = 1;
    }
    printf("DONE n %d k %d firsts %d nodes %lld solutions %lld prep %.1f seconds %.1f\n", n, k, nfirst, nodes, sols,
           tprep, (double)(clock() - t0) / CLOCKS_PER_SEC);
    return 0;
}
