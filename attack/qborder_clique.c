/* Border problem of the Z_q-bordered propus (propus.qborder_solve): N = 4q rows r in {+-1}^N (first entry +1),
   pairwise orthogonal, with r_i Mc r_j^T = +-N for all i, j. Vertices: the r with r Mc r^T = +-N; edges: orthogonal
   pairs with cross term +-N; a solution is an N-clique. For each start vertex s (random order) the subgraph on N(s)
   gets bitset adjacency and is searched for an (N - 1)-clique with greedy-colouring bounds (Tomita's MCQ).
   Input on stdin: N, then the N x N matrix Mc. Args: time limit (s), seed, max starts.
   Output: "SOLVED ..." and the N rows, or "TIMEOUT ..." / "NONE ..." (NONE: every start exhausted, no clique). */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdint.h>
static int N, M[32][32];
static uint32_t *V; static int nv;
static int **adj; static int *deg;
static double t0, tlim;
static double now(void) { struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts); return ts.tv_sec + 1e-9 * ts.tv_nsec; }
static inline int sgn(uint32_t r, int k) { return (r >> k) & 1 ? -1 : 1; }
/* local graph */
static int ln; static size_t lw; static uint64_t *LA; static int *lmap;
static int target, found; static int best[64], cur[64]; static long nodes;
#define ROW(i) (LA + (size_t)(i) * lw)
static int popc(const uint64_t *b) { int c = 0; for (int w = 0; w < lw; w++) c += __builtin_popcountll(b[w]); return c; }
static void expand(uint64_t *P, int depth) {
    if (found) return;
    if ((++nodes & 0x3ff) == 0 && now() - t0 > tlim) { found = -1; return; }
    /* greedy coloring of P: order vertices by color classes */
    int np_ = popc(P); if (depth + np_ < target) return;
    int *order = malloc(sizeof(int) * np_), *col = malloc(sizeof(int) * np_);
    uint64_t *U = malloc(8 * lw), *Q = malloc(8 * lw);
    memcpy(U, P, 8 * lw);
    int k = 0, color = 0;
    while (k < np_) {
        color++; memcpy(Q, U, 8 * lw);
        for (int w = 0; w < lw; w++) while (Q[w]) {
            int b = __builtin_ctzll(Q[w]); int v = w * 64 + b;
            Q[w] &= Q[w] - 1;
            U[v >> 6] &= ~(1ULL << (v & 63));
            /* remove neighbours of v from Q */
            const uint64_t *r = ROW(v); for (int x = 0; x < lw; x++) Q[x] &= ~r[x];
            order[k] = v; col[k] = color; k++;
        }
    }
    uint64_t *NP = malloc(8 * lw);
    for (int i = np_ - 1; i >= 0 && !found; i--) {
        if (depth + col[i] < target) break;
        int v = order[i];
        cur[depth] = v;
        if (depth + 1 == target) { memcpy(best, cur, sizeof(int) * target); found = 1; break; }
        const uint64_t *r = ROW(v); for (int x = 0; x < lw; x++) NP[x] = P[x] & r[x];
        expand(NP, depth + 1);
        P[v >> 6] &= ~(1ULL << (v & 63));
    }
    free(order); free(col); free(U); free(Q); free(NP);
}
int main(int argc, char **argv) {
    tlim = argc > 1 ? atof(argv[1]) : 600; unsigned seed = argc > 2 ? atoi(argv[2]) : 1; long maxstarts = argc > 3 ? atol(argv[3]) : 1000000;
    if (scanf("%d", &N) != 1) return 1;
    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) if (scanf("%d", &M[i][j]) != 1) return 1;
    t0 = now();
    uint32_t total = 1u << (N - 1);
    V = malloc(sizeof(uint32_t) * total); nv = 0;
    for (uint32_t x = 0; x < total; x++) {
        uint32_t r = x << 1; long q = 0;
        for (int i = 0; i < N; i++) { long s = 0; for (int j = 0; j < N; j++) if (M[i][j]) s += M[i][j] * sgn(r, j); q += sgn(r, i) * s; }
        if (q == N || q == -N) V[nv++] = r;
    }
    int8_t *U = malloc((size_t)nv * N);
    for (int a = 0; a < nv; a++) for (int i = 0; i < N; i++) { int s = 0; for (int j = 0; j < N; j++) s += M[i][j] * sgn(V[a], j); U[(size_t)a * N + i] = (int8_t)s; }
    adj = malloc(sizeof(int *) * nv); deg = calloc(nv, sizeof(int));
    int cap = 4096; int *buf = malloc(sizeof(int) * cap); long edges = 0;
    for (int a = 0; a < nv; a++) {
        int k = 0; const int8_t *u = U + (size_t)a * N;
        for (int b = 0; b < nv; b++) {
            if (b == a || __builtin_popcount(V[a] ^ V[b]) * 2 != N) continue;
            int s = 0; uint32_t rb = V[b];
            for (int i = 0; i < N; i++) s += (rb >> i) & 1 ? -u[i] : u[i];
            if (s == N || s == -N) { if (k == cap) { cap *= 2; buf = realloc(buf, sizeof(int) * cap); } buf[k++] = b; }
        }
        adj[a] = malloc(sizeof(int) * (k ? k : 1)); memcpy(adj[a], buf, sizeof(int) * k); deg[a] = k; edges += k;
    }
    fprintf(stderr, "vertices %d, mean degree %.1f [%.1fs]\n", nv, (double)edges / nv, now() - t0);
    srand(seed);
    int *perm = malloc(sizeof(int) * nv); for (int a = 0; a < nv; a++) perm[a] = a;
    for (int a = nv - 1; a > 0; a--) { int b = rand() % (a + 1); int t = perm[a]; perm[a] = perm[b]; perm[b] = t; }
    int *pos = malloc(sizeof(int) * nv);
    target = N - 1; found = 0;
    long starts = 0;
    for (int si = 0; si < nv && !found && starts < maxstarts; si++, starts++) {
        int s = perm[si];
        ln = deg[s]; lw = (ln + 63) / 64; lmap = adj[s];
        if (ln < target) continue;
        for (int a = 0; a < nv; a++) pos[a] = -1;
        for (int i = 0; i < ln; i++) pos[lmap[i]] = i;
        LA = calloc((size_t)ln * lw, 8);
        for (int i = 0; i < ln; i++) { int u = lmap[i]; for (int j = 0; j < deg[u]; j++) { int p = pos[adj[u][j]]; if (p >= 0) ROW(i)[p >> 6] |= 1ULL << (p & 63); } }
        uint64_t *P = calloc(lw, 8); for (int i = 0; i < ln; i++) P[i >> 6] |= 1ULL << (i & 63);
        expand(P, 0);
        free(P); free(LA);
        if (found == 1) {
            printf("SOLVED after %ld starts, %ld nodes\n", starts + 1, nodes);
            for (int i = 0; i < N; i++) printf("%d ", sgn(V[s], i));
            printf("\n");
            for (int d = 0; d < target; d++) { uint32_t r = V[lmap[best[d]]]; for (int i = 0; i < N; i++) printf("%d ", sgn(r, i)); printf("\n"); }
        }
        if (found == -1) break;
        /* exclude s from later local graphs: a clique through s was not found */
    }
    if (found != 1) printf("%s after %ld starts, %ld nodes\n", found == -1 ? "TIMEOUT" : "NONE", starts, nodes);
    fprintf(stderr, "[%.1fs]\n", now() - t0);
    return 0;
}
