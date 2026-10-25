// BH(n,6): n x n matrices of sixth roots with X X^* = n I; normalised first row and column = 1; prints all (rows 2..n).
#include <stdio.h>
#include <stdlib.h>
typedef struct { int x, y; } Z;
static Z U6[6] = { {1,0}, {1,1}, {0,1}, {-1,0}, {-1,-1}, {0,-1} };
static inline Z zmul(Z a, Z b) { Z r = { a.x*b.x - a.y*b.y, a.x*b.y + a.y*b.x - a.y*b.y }; return r; }
int n, X[12][12]; long long sols = 0, maxp;
static int orth(int r, int s) { Z t = {0,0}; for (int j = 0; j < n; j++) { Z a = zmul(U6[X[r][j]], U6[(6 - X[s][j]) % 6]); t.x += a.x; t.y += a.y; } return t.x == 0 && t.y == 0; }
static void go(int r, int j) {
  if (j == n) { for (int s = 0; s < r; s++) if (!orth(r, s)) return; if (r + 1 == n) { sols++; if (sols <= maxp) { for (int a = 1; a < n; a++) for (int b = 1; b < n; b++) printf("%d", X[a][b]); printf("\n"); } return; } X[r+1][0] = 0; go(r + 1, 1); return; }
  for (int c = 0; c < 6; c++) { X[r][j] = c; go(r, j + 1); }
}
int main(int argc, char **argv) { n = atoi(argv[1]); maxp = argc > 2 ? atoll(argv[2]) : 10; for (int j = 0; j < n; j++) X[0][j] = 0; X[1][0] = 0; go(1, 1); fprintf(stderr, "n=%d normalised BH(n,6): %lld\n", n, sols); return 0; }
