/* Border problem of the Z_q-bordered propus (propus.qborder_solve_c), N = 4q <= 32: N rows r in {+-1}^N (first entry
   +1), pairwise orthogonal, with r_i Mc r_j^T = +-N for all i, j. Vertices: the r with r Mc r^T = +-N (bit k of r set
   when entry k is -1); edges: orthogonal pairs with cross term +-N; a solution is an N-clique. No global adjacency:
   for random adjacent start pairs (s, t), C = N(s) & N(t) is found by scanning, its adjacency is built as bitsets
   and an (N - 2)-clique is searched with greedy-colouring bounds (Tomita's MCQ). q = 5 (91,648 vertices) is solved in
   about a second; q = 7 has 17.4 million vertices and |C| near 15,000.
   Input on stdin: N, then the N x N matrix Mc. Args: time limit (s), seed, largest |C| accepted.
   Output: "SOLVED ..." and the N rows, or "TIMEOUT ...". */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdint.h>
static int N, M[32][32];
static uint32_t *V; static long nv;
static double t0, tlim;
static double now(void) { struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts); return ts.tv_sec + 1e-9 * ts.tv_nsec; }
static inline int sgn(uint32_t r, int k) { return (r >> k) & 1 ? -1 : 1; }
static void mcrow(uint32_t r, int8_t *u) { for (int i = 0; i < N; i++) { int s = 0; for (int j = 0; j < N; j++) if (M[i][j]) s += M[i][j] * sgn(r, j); u[i] = (int8_t)s; } }
static inline int adjq(const int8_t *u, uint32_t a, uint32_t b) {   /* u = Mc a */
    if (__builtin_popcount(a ^ b) * 2 != N) return 0;
    int s = 0; for (int i = 0; i < N; i++) s += (b >> i) & 1 ? -u[i] : u[i];
    return s == N || s == -N;
}
static int ln; static size_t lw; static uint64_t *LA;
static int target, found; static int best[64], cur[64]; static long nodes;
#define ROW(i) (LA + (size_t)(i) * lw)
static int popc(const uint64_t *b) { int c = 0; for (size_t w = 0; w < lw; w++) c += __builtin_popcountll(b[w]); return c; }
static void expand(uint64_t *P, int depth) {
    if (found) return;
    if ((++nodes & 0x3ff) == 0 && now() - t0 > tlim) { found = -1; return; }
    int np_ = popc(P); if (depth + np_ < target) return;
    int *order = malloc(sizeof(int) * np_), *col = malloc(sizeof(int) * np_);
    uint64_t *U = malloc(8 * lw), *Q = malloc(8 * lw), *NP = malloc(8 * lw);
    memcpy(U, P, 8 * lw);
    int k = 0, color = 0;
    while (k < np_) {
        color++; memcpy(Q, U, 8 * lw);
        for (size_t w = 0; w < lw; w++) while (Q[w]) {
            int b = __builtin_ctzll(Q[w]); int v = (int)(w * 64 + b);
            Q[w] &= Q[w] - 1; U[v >> 6] &= ~(1ULL << (v & 63));
            const uint64_t *r = ROW(v); for (size_t x = 0; x < lw; x++) Q[x] &= ~r[x];
            order[k] = v; col[k] = color; k++;
        }
    }
    for (int i = np_ - 1; i >= 0 && !found; i--) {
        if (depth + col[i] < target) break;
        int v = order[i]; cur[depth] = v;
        if (depth + 1 == target) { memcpy(best, cur, sizeof(int) * target); found = 1; break; }
        const uint64_t *r = ROW(v); for (size_t x = 0; x < lw; x++) NP[x] = P[x] & r[x];
        expand(NP, depth + 1);
        P[v >> 6] &= ~(1ULL << (v & 63));
    }
    free(order); free(col); free(U); free(Q); free(NP);
}
int main(int argc, char **argv) {
    tlim = argc > 1 ? atof(argv[1]) : 600; srand(argc > 2 ? atoi(argv[2]) : 1); long maxC = argc > 3 ? atol(argv[3]) : 60000;
    if (scanf("%d", &N) != 1) return 1;
    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) if (scanf("%d", &M[i][j]) != 1) return 1;
    t0 = now();
    uint64_t total = 1ULL << (N - 1); long cap = 1 << 20; V = malloc(sizeof(uint32_t) * cap); nv = 0;
    int8_t u[32];
    for (uint64_t x = 0; x < total; x++) {
        uint32_t r = (uint32_t)(x << 1); long q = 0;
        mcrow(r, u); for (int i = 0; i < N; i++) q += sgn(r, i) * u[i];
        if (q == N || q == -N) { if (nv == cap) { cap *= 2; V = realloc(V, sizeof(uint32_t) * cap); } V[nv++] = r; }
    }
    fprintf(stderr, "vertices %ld [%.1fs]\n", nv, now() - t0);
    uint32_t *Ns = malloc(sizeof(uint32_t) * nv), *C = malloc(sizeof(uint32_t) * nv);
    long starts = 0; target = N - 2; found = 0;
    while (!found && now() - t0 < tlim) {
        uint32_t s = V[((long)rand() * RAND_MAX + rand()) % nv]; int8_t us[32]; mcrow(s, us);
        long ns = 0; for (long a = 0; a < nv; a++) if (adjq(us, s, V[a])) Ns[ns++] = V[a];
        if (ns < N - 1) continue;
        for (int tries = 0; tries < 20 && !found && now() - t0 < tlim; tries++) {
            uint32_t t = Ns[((long)rand() * RAND_MAX + rand()) % ns]; int8_t ut[32]; mcrow(t, ut);
            long nc = 0; for (long a = 0; a < ns; a++) if (Ns[a] != t && adjq(ut, t, Ns[a])) C[nc++] = Ns[a];
            starts++;
            if (nc < N - 2 || nc > maxC) { if (starts % 50 == 0) fprintf(stderr, "start %ld |N(s)| %ld |C| %ld [%.0fs]\n", starts, ns, nc, now() - t0); continue; }
            ln = (int)nc; lw = (ln + 63) / 64; LA = calloc((size_t)ln * lw, 8);
            for (int i = 0; i < ln; i++) { int8_t ui[32]; mcrow(C[i], ui); for (int j = i + 1; j < ln; j++) if (adjq(ui, C[i], C[j])) { ROW(i)[j >> 6] |= 1ULL << (j & 63); ROW(j)[i >> 6] |= 1ULL << (i & 63); } }
            uint64_t *P = calloc(lw, 8); for (int i = 0; i < ln; i++) P[i >> 6] |= 1ULL << (i & 63);
            nodes = 0; expand(P, 0);
            if (starts % 10 == 0 || found) fprintf(stderr, "start %ld |N(s)| %ld |C| %ld nodes %ld found %d [%.0fs]\n", starts, ns, nc, nodes, found, now() - t0);
            if (found == 1) {
                printf("SOLVED after %ld starts\n", starts);
                uint32_t rows[64]; rows[0] = s; rows[1] = t; for (int d = 0; d < target; d++) rows[d + 2] = C[best[d]];
                for (int d = 0; d < N; d++) { for (int i = 0; i < N; i++) printf("%d ", sgn(rows[d], i)); printf("\n"); }
            }
            free(P); free(LA);
            if (found == -1) found = 0;      /* time limit inside this start: continue with the next if time remains */
        }
    }
    if (found != 1) printf("TIMEOUT after %ld starts\n", starts);
    return 0;
}
