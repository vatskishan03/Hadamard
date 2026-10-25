// Templates on a given support graph: Psi = D + 2 Phi over Z[w], Phi_uv sixth root on edges, 0 off the support;
// Psi^2 = (1 + 4 deg) I (graph deg-regular). Gauge: Phi on a BFS spanning tree = 1. usage: kg n a "u-v u-v ..." [maxprint]
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef struct { int x, y; } Z;
static Z U6[6] = { {1,0}, {1,1}, {0,1}, {-1,0}, {-1,-1}, {0,-1} };
static inline Z zmul(Z a, Z b) { Z r = { a.x*b.x - a.y*b.y, a.x*b.y + a.y*b.x - a.y*b.y }; return r; }
int n, A, D[24], adj[24][24], ph[24][24], tree[24][24]; long long nodes = 0, sols = 0, maxp;
static Z ent(int u, int v) { return (!adj[u][v] || ph[u][v] < 0) ? (Z){0,0} : U6[ph[u][v]]; }
static int pairok(int u, int v) {
  Z s = {0,0};
  for (int k = 0; k < n; k++) if (k != u && k != v) { Z t = zmul(ent(u,k), ent(k,v)); s.x += t.x; s.y += t.y; }
  Z f = ent(u, v); int c = D[u] + D[v];
  return 2*c*f.x + 4*s.x == 0 && 2*c*f.y + 4*s.y == 0;
}
static void row(int u, int j) {
  nodes++;
  if (j == n) {
    for (int v = 0; v < u; v++) if (!pairok(v, u)) return;
    if (u + 1 == n) { for (int v = 0; v < n; v++) if (!pairok(v, v) && 0) return; sols++; if (sols <= maxp) { printf("SOL"); for (int a = 0; a < n; a++) for (int b = a+1; b < n; b++) if (adj[a][b]) printf(" %d", ph[a][b]); printf("\n"); } return; }
    row(u + 1, u + 2); return;
  }
  if (!adj[u][j] || tree[u][j]) { row(u, j + 1); return; }
  for (int c = 0; c < 6; c++) { ph[u][j] = c; ph[j][u] = (6 - c) % 6; row(u, j + 1); }
  ph[u][j] = ph[j][u] = -1;
}
int main(int argc, char **argv) {
  n = atoi(argv[1]); A = atoi(argv[2]); maxp = argc > 4 ? atoll(argv[4]) : 5;
  memset(adj, 0, sizeof adj); memset(tree, 0, sizeof tree);
  char *s = strdup(argv[3]); for (char *t = strtok(s, " "); t; t = strtok(NULL, " ")) { int a, b; sscanf(t, "%d-%d", &a, &b); adj[a][b] = adj[b][a] = 1; }
  for (int i = 0; i < n; i++) { D[i] = i < A ? 1 : -1; for (int j = 0; j < n; j++) ph[i][j] = -1; }
  // BFS tree from 0; tree edges get phase 0
  int seen[24] = {0}, q[24], h = 0, tl = 0; q[tl++] = 0; seen[0] = 1;
  while (h < tl) { int u = q[h++]; for (int v = 0; v < n; v++) if (adj[u][v] && !seen[v]) { seen[v] = 1; q[tl++] = v; tree[u][v] = tree[v][u] = 1; ph[u][v] = ph[v][u] = 0; } }
  row(0, 1);
  printf("n=%d a=%d nodes=%lld sols=%lld\n", n, A, nodes, sols);
  return 0;
}
