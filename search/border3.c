// Border solver for the bordered Goethals-Seidel family with Q = Z_3 and T = all characters (search/bordered.py).
// Given X (12 x 12 Hadamard) and W = Mc/2 (12 x 12, entries 0/+-1, 3 nonzeros per column; Mc is the Goethals-Seidel
// array over Z_3 of the four compressions), find the rows c_0..c_11 of Ch (+-1, mutually orthogonal) such that
// sum_k W[k][j] c_k lies in Omega = { X^T s / 2 : s in {+-1}^12 } (entries +-1, +-3) for every column j. Then
// Rh = -X Ch^T W / 6 is a +-1 Hadamard matrix and X Ch^T + Rh Mc^T = 0, which is the border condition.
// Depth-first search over the rows of Ch, ordered so that the columns of W close early.
// Input on stdin: m (= 12), then X (m*m ints), then W (m*m ints); argv[1] = node limit.
// Prints "FOUND r_0 .. r_11" (rows as 12-bit masks, bit i set = entry -1) or "NONE nodes" / "LIMIT nodes".
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#define M 12
static int m, X[M][M], Wm[M][M];
static uint8_t *omega;                   // bitmap over 4^m keys
static int ord[M], pos[M];
static int ecnt[M], eidx[M][M];          // closing edges at level i
static int supk[M][3], sups[M][3];       // edge j: rows and signs
static int chosen[M];
static long long nodes = 0, maxnodes;
static int orthlist[1<<M][1000], orthn[1<<M];
static int vecv(int mask, int i) { return (mask >> i & 1) ? -1 : 1; }
static uint32_t key_of(const int *v) { uint32_t k = 0; for (int i = 0; i < m; i++) { int c = (v[i] + 3) / 2; k = k * 4 + c; } return k; }
static int found = 0;
static int check_edge(int j) {
    int v[M];
    for (int i = 0; i < m; i++) { v[i] = 0; for (int t = 0; t < 3; t++) v[i] += sups[j][t] * vecv(chosen[supk[j][t]], i); }
    for (int i = 0; i < m; i++) if (v[i] != 1 && v[i] != -1 && v[i] != 3 && v[i] != -3) return 0;
    uint32_t k = key_of(v); return omega[k >> 3] >> (k & 7) & 1;
}
static void rec(int lev) {
    if (found) return;
    if (++nodes > maxnodes) return;
    if (lev == m) { found = 1; return; }
    int r = ord[lev];
    // candidates: orthogonal to all chosen rows so far
    int base = -1;
    for (int l = 0; l < lev; l++) { base = chosen[ord[l]]; break; }
    int ncand; const int *cand; static int all[1<<M]; static int init = 0;
    if (!init) { for (int i = 0; i < (1<<m); i++) all[i] = i; init = 1; }
    if (base < 0) { cand = all; ncand = 1 << (m - 1); }      // first row: entry 0 is +1 (global sign)
    else { cand = orthlist[base]; ncand = orthn[base]; }
    for (int ci = 0; ci < ncand && !found; ci++) {
        int c = cand[ci], ok = 1;
        for (int l = 1; l < lev; l++) if (__builtin_popcount(c ^ chosen[ord[l]]) != m / 2) { ok = 0; break; }
        if (!ok) continue;
        chosen[r] = c;
        for (int e = 0; e < ecnt[lev]; e++) if (!check_edge(eidx[lev][e])) { ok = 0; break; }
        if (ok) rec(lev + 1);
    }
}
int main(int argc, char **argv) {
    maxnodes = argc > 1 ? atoll(argv[1]) : 2000000000LL;
    if (scanf("%d", &m) != 1 || m != M) return 2;
    for (int i = 0; i < m; i++) for (int j = 0; j < m; j++) if (scanf("%d", &X[i][j]) != 1) return 2;
    for (int i = 0; i < m; i++) for (int j = 0; j < m; j++) if (scanf("%d", &Wm[i][j]) != 1) return 2;
    omega = calloc((1u << (2 * m)) / 8 + 1, 1);
    int cnt = 0;
    for (int s = 0; s < (1 << m); s++) {
        int y[M], ok = 1;
        for (int j = 0; j < m; j++) { int t = 0; for (int i = 0; i < m; i++) t += X[i][j] * vecv(s, i); if (t % 2) { ok = 0; break; } y[j] = t / 2; if (y[j] != 1 && y[j] != -1 && y[j] != 3 && y[j] != -3) { ok = 0; break; } }
        if (!ok) continue;
        uint32_t k = key_of(y); omega[k >> 3] |= 1 << (k & 7); cnt++;
    }
    fprintf(stderr, "omega %d\n", cnt);
    for (int a = 0; a < (1 << m); a++) { orthn[a] = 0; for (int b = 0; b < (1 << m); b++) if (__builtin_popcount(a ^ b) == m / 2) orthlist[a][orthn[a]++] = b; }
    for (int j = 0; j < m; j++) { int t = 0; for (int k = 0; k < m; k++) if (Wm[k][j]) { if (t == 3) return 3; supk[j][t] = k; sups[j][t] = Wm[k][j]; t++; } if (t != 3) return 3; }
    // greedy order: maximize closed edges
    int used[M] = {0};
    for (int lev = 0; lev < m; lev++) {
        int best = -1, bs = -1;
        for (int r = 0; r < m; r++) if (!used[r]) {
            int sc = 0;
            for (int j = 0; j < m; j++) { int in = 0, rest = 1; for (int t = 0; t < 3; t++) { if (supk[j][t] == r) in = 1; else if (!used[supk[j][t]]) rest = 0; } if (in && rest) sc += 100; else if (in) { int u = 0; for (int t = 0; t < 3; t++) u += used[supk[j][t]]; sc += u; } }
            if (sc > bs) { bs = sc; best = r; }
        }
        ord[lev] = best; used[best] = 1; pos[best] = lev;
    }
    for (int lev = 0; lev < m; lev++) ecnt[lev] = 0;
    for (int j = 0; j < m; j++) { int mx = 0; for (int t = 0; t < 3; t++) if (pos[supk[j][t]] > mx) mx = pos[supk[j][t]]; eidx[mx][ecnt[mx]++] = j; }
    rec(0);
    if (found) { printf("FOUND"); for (int r = 0; r < m; r++) printf(" %d", chosen[r]); printf("\n"); }
    else printf("%s %lld\n", nodes > maxnodes ? "LIMIT" : "NONE", nodes);
    return 0;
}
