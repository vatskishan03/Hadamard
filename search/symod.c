/* Full symmetric orthogonal designs in Balonin (propus) form over Z_n with K variables.
 * Blocks A (symmetric), B (= C), D are circulants with entries +-x_i (every entry one variable), and with
 * c_ij(s) = sum_g [a_i(g) a_j(g+s) + 2 b_i(g) b_j(g+s) + d_i(g) d_j(g+s)] the design conditions are
 *   c_ii(s) = 0 and c_ij(s) + c_ji(s) = 0 for s != 0 (enough when every plug-in is symmetric),
 * "strong" pairs additionally have c_ij(s) = 0 for every s != 0 (enough for arbitrary plug-ins).
 * Weights s_i = c_ii(0) = (#x_i in A) + 2 (#x_i in B) + (#x_i in D); s_1 + .. + s_K = 4n.
 * usage: symod n K
 */
#include <stdio.h>
#include <stdlib.h>

static int n, K, na, ent[64];            /* entries: A half (na values), B (n), D (n); value = +-(var+1) */
static long long found;

static void expand(int *A, int *B, int *D) {
    for (int g = 0; g < n; g++) { int h = g <= n - g ? g : n - g; A[g] = ent[h]; }
    for (int g = 0; g < n; g++) { B[g] = ent[na + g]; D[g] = ent[na + n + g]; }
}

static int check(int print) {
    int A[32], B[32], D[32];
    expand(A, B, D);
    int w[8] = {0};
    for (int g = 0; g < n; g++) { w[abs(A[g]) - 1] += 1; w[abs(B[g]) - 1] += 2; w[abs(D[g]) - 1] += 1; }
    for (int i = 0; i < K; i++) if (!w[i]) return 0;           /* every variable used */
    int strong[8][8] = {{0}};
    for (int i = 0; i < K; i++) for (int j = 0; j < K; j++) strong[i][j] = 1;
    for (int s = 1; s < n; s++) {
        int c[8][8] = {{0}};
        for (int g = 0; g < n; g++) {
            int t = (g + s) % n;
            int X[3][2] = {{A[g], A[t]}, {B[g], B[t]}, {D[g], D[t]}}, wt[3] = {1, 2, 1};
            for (int k = 0; k < 3; k++) {
                int i = abs(X[k][0]) - 1, j = abs(X[k][1]) - 1;
                int sg = (X[k][0] > 0 ? 1 : -1) * (X[k][1] > 0 ? 1 : -1);
                c[i][j] += wt[k] * sg;
            }
        }
        for (int i = 0; i < K; i++) {
            if (c[i][i]) return 0;
            for (int j = i + 1; j < K; j++) {
                if (c[i][j] + c[j][i]) return 0;
                if (c[i][j]) strong[i][j] = strong[j][i] = 0;
            }
        }
    }
    found++;
    if (print) {
        printf("weights");
        for (int i = 0; i < K; i++) printf(" %d", w[i]);
        printf("  strong pairs:");
        for (int i = 0; i < K; i++) for (int j = i + 1; j < K; j++) if (strong[i][j]) printf(" %d%d", i + 1, j + 1);
        printf("  A");
        for (int g = 0; g < n; g++) printf(" %+d", A[g]);
        printf("  B"); for (int g = 0; g < n; g++) printf(" %+d", B[g]);
        printf("  D"); for (int g = 0; g < n; g++) printf(" %+d", D[g]);
        printf("\n");
    }
    return 1;
}

static long long seenw[64][64][64];
static void rec(int pos, int tot, int maxvar) {
    if (pos == tot) {
        int A[32], B[32], D[32]; expand(A, B, D);
        if (check(0)) {
            int w[8] = {0};
            for (int g = 0; g < n; g++) { w[abs(A[g]) - 1] += 1; w[abs(B[g]) - 1] += 2; w[abs(D[g]) - 1] += 1; }
            long long *cnt = &seenw[w[0]][K > 1 ? w[1] : 0][K > 2 ? w[2] : 0];
            if ((*cnt)++ < 2) check(1);
        }
        return;
    }
    /* symmetry breaking: variables introduced in order (first use of var i+1 after var i); first entry +x1 */
    for (int v = 1; v <= maxvar + 1 && v <= K; v++) for (int sg = -1; sg <= 1; sg += 2) {
        if (pos == 0 && sg < 0) continue;
        if (v == maxvar + 1 && sg < 0) continue;               /* sign of a variable's first occurrence fixed */
        ent[pos] = sg * v;
        rec(pos + 1, tot, v > maxvar ? v : maxvar);
    }
}

int main(int argc, char **argv) {
    n = atoi(argv[1]); K = atoi(argv[2]);
    na = n / 2 + 1;
    rec(0, na + 2 * n, 0);
    printf("total designs %lld (n=%d, K=%d)\n", found, n, K);
    return 0;
}
