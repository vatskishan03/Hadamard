// Complete-graph templates: n x n Hermitian Psi = D + 2 Phi over Z[w], Phi_uv sixth roots of unity (u != v), D = diag(+1 x a, -1 x (n-a)),
// Psi^2 = (1 + 4(n-1)) I. Gauge: Phi[0][j] = 1. Row-by-row backtracking, checking (Psi^2)_{uv} once rows u, v are complete.
// usage: kn n a [maxprint]
#include <stdio.h>
#include <stdlib.h>
typedef struct { int x, y; } Z;                       // x + y w
static Z U6[6];                                        // 1, -w^2 = 1 + w, w, -1, w^2 = -1 - w, -w
static inline Z zmul(Z a, Z b) { Z r = { a.x*b.x - a.y*b.y, a.x*b.y + a.y*b.x - a.y*b.y }; return r; }
static inline Z zconj(Z a) { Z r = { a.x - a.y, -a.y }; return r; }
int n, A, D[16]; int ph[16][16]; long long nodes = 0, sols = 0; int maxprint;
static Z ent(int u, int v) { return ph[u][v] < 0 ? (Z){0,0} : U6[ph[u][v]]; }
static int pairok(int u, int v) {      // (Psi^2)_{uv} = 2(Du + Dv) Phi_uv + 4 sum_{k != u,v} Phi_uk Phi_kv = 0
  Z s = {0, 0};
  for (int k = 0; k < n; k++) if (k != u && k != v) { Z t = zmul(ent(u, k), ent(k, v)); s.x += t.x; s.y += t.y; }
  Z f = ent(u, v); int c = D[u] + D[v];
  return 2 * c * f.x + 4 * s.x == 0 && 2 * c * f.y + 4 * s.y == 0;
}
static void row(int u, int j) {
  nodes++;
  if (j == n) {          // row u complete (entries below diagonal known from earlier rows)
    for (int v = 0; v < u; v++) if (!pairok(v, u)) return;
    if (u + 1 == n) {
      sols++;
      if (sols <= maxprint) { printf("SOL"); for (int a = 0; a < n; a++) for (int b = a + 1; b < n; b++) printf(" %d", ph[a][b]); printf("\n"); }
      return;
    }
    row(u + 1, u + 2); return;
  }
  for (int c = 0; c < 6; c++) {
    ph[u][j] = c; ph[j][u] = (6 - c) % 6;   // conj of w-power: unit index k means e^{i pi k/3}; conj is -k
    row(u, j + 1);
  }
  ph[u][j] = ph[j][u] = -1;
}
int main(int argc, char **argv) {
  n = atoi(argv[1]); A = atoi(argv[2]); maxprint = argc > 3 ? atoi(argv[3]) : 5;
  // U6[k] = e^{i pi k / 3} in basis (1, w), w = e^{2 pi i/3}: e^{i pi/3} = 1 + w
  Z t[6] = { {1,0}, {1,1}, {0,1}, {-1,0}, {-1,-1}, {0,-1} }; for (int k = 0; k < 6; k++) U6[k] = t[k];
  for (int i = 0; i < n; i++) { D[i] = i < A ? 1 : -1; for (int j = 0; j < n; j++) ph[i][j] = -1; }
  for (int j = 1; j < n; j++) { ph[0][j] = 0; ph[j][0] = 0; }
  if (n == 2) { printf("n=2\n"); return 0; }
  row(1, 2);
  printf("n=%d a=%d nodes=%lld sols=%lld\n", n, A, nodes, sols);
  return 0;
}
