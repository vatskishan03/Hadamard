// Decoupled bipartite cyclic flag structures (NOTES.md §96) on interval supports: row R1 has atoms at offsets mu*j,
// j = 1..d, row R2 at offsets -mu*(d+j), j = 1..d (as at p = 17, where every decoupled solution is of this form up to
// multipliers). Each atom gets a B-class y in {0,1} and a phase c in Z3; coefficient kappa = sigma(offset) w^c,
// sigma = chi * flip. Row condition (decoupled, one orientation): for every h != 0,
//   F(h) = sum over atom pairs i != j of the same class with delta_i - delta_j = h of kappa_i conj(kappa_j) = 0.
// Cross condition: for every h, sum over (i in R1, j in R2), same class, delta_i - delta_j = h, of kappa_i conj(kappa'_j w^g) = 0.
// usage: intcyc p mu flip [maxprint]
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef struct { int x, y; } Z;
static inline Z zmul(Z a, Z b) { Z r = { a.x*b.x - a.y*b.y, a.x*b.y + a.y*b.x - a.y*b.y }; return r; }
static inline Z zconj(Z a) { Z r = { a.x - a.y, -a.y }; return r; }
static const Z WP[3] = { {1,0}, {0,1}, {-1,-1} };
int p, d, mu, flip, sig[1024];
typedef struct { int cls[64]; int ph[64]; } Row;
Row *rows[2]; long nrows[2], cap[2];
int del[2][64]; int CL[64], PH[64];
int hn[1024], hi_[1024][128], hj_[1024][128], hlast[1024]; int hlist_at[64][1024], nh_at[64];
static Z kap(int dl, int c) { Z w = WP[c]; Z r = { w.x * sig[dl], w.y * sig[dl] }; return r; }
static int R;   // which row is being enumerated
static void store() {
  if (nrows[R] == cap[R]) { cap[R] = cap[R] ? 2 * cap[R] : 1 << 16; rows[R] = realloc(rows[R], cap[R] * sizeof(Row)); }
  Row *r = &rows[R][nrows[R]++]; memcpy(r->cls, CL, sizeof(int) * d); memcpy(r->ph, PH, sizeof(int) * d);
}
static void bt(int k) {
  if (k == d) { store(); return; }
  for (int y = 0; y < 2; y++) for (int c = 0; c < 3; c++) {
    if (k == 0 && (c || y)) continue;          // gauge: first atom phase 0, class 0 (class swap symmetry)
    CL[k] = y; PH[k] = c; int ok = 1;
    for (int t = 0; t < nh_at[k] && ok; t++) {
      int h = hlist_at[k][t]; Z F = {0, 0};
      for (int q = 0; q < hn[h]; q++) {
        int i = hi_[h][q], j = hj_[h][q];
        if (CL[i] != CL[j]) continue;
        Z t2 = zmul(kap(del[R][i], PH[i]), zconj(kap(del[R][j], PH[j]))); F.x += t2.x; F.y += t2.y;
      }
      if (F.x || F.y) ok = 0;
    }
    if (ok) bt(k + 1);
  }
}
static void enumrow(int r) {
  R = r;
  for (int h = 0; h < p; h++) { hn[h] = 0; hlast[h] = -1; }
  for (int i = 0; i < d; i++) for (int j = 0; j < d; j++) if (i != j) {
    int h = ((del[r][i] - del[r][j]) % p + p) % p; hi_[h][hn[h]] = i; hj_[h][hn[h]] = j; hn[h]++;
    int mx = i > j ? i : j; if (mx > hlast[h]) hlast[h] = mx;
  }
  for (int t = 0; t < d; t++) nh_at[t] = 0;
  for (int h = 1; h < p; h++) if (hn[h]) { int t = hlast[h]; hlist_at[t][nh_at[t]++] = h; }
  bt(0);
}
int main(int argc, char **argv) {
  p = atoi(argv[1]); mu = atoi(argv[2]); flip = atoi(argv[3]); int maxprint = argc > 4 ? atoi(argv[4]) : 5;
  d = (p - 1) / 4;
  int qr[1024] = {0}; for (int x = 1; x < p; x++) qr[x * x % p] = 1;
  for (int x = 1; x < p; x++) sig[x] = (qr[x] ? 1 : -1) * (flip ? -1 : 1);
  // alternating order 1, d, 2, d-1, ... so that the long shifts close early
  for (int k = 0; k < d; k++) {
    int j = (k % 2 == 0) ? k / 2 : d - 1 - k / 2;
    del[0][k] = (mu * (j + 1)) % p; del[1][k] = ((p - mu * (d + j + 1) % p) % p + p) % p;
  }
  enumrow(0); enumrow(1);
  fprintf(stderr, "p=%d mu=%d flip=%d rows: %ld %ld\n", p, mu, flip, nrows[0], nrows[1]);
  long nsol = 0;
  // the class labels of row 2 may be swapped relative to row 1: row 2 was gauged with its first atom in class 0
  for (long a = 0; a < nrows[0]; a++) for (long b = 0; b < nrows[1]; b++) for (int sw = 0; sw < 2; sw++) {
    Row *r1 = &rows[0][a], *r2 = &rows[1][b];
    int n1 = 0; for (int i = 0; i < d; i++) { n1 += r1->cls[i] == 0; n1 += (r2->cls[i] ^ sw) == 0; }
    if (n1 != d) continue;                       // B-degree balance: each B class gets d atoms
    for (int g = 0; g < 3; g++) {
      static Z X[1024]; for (int h = 0; h < p; h++) X[h].x = X[h].y = 0;
      for (int i = 0; i < d; i++) for (int j = 0; j < d; j++) {
        if (r1->cls[i] != (r2->cls[j] ^ sw)) continue;
        int h = ((del[0][i] - del[1][j]) % p + p) % p;
        Z t2 = zmul(kap(del[0][i], r1->ph[i]), zconj(kap(del[1][j], (r2->ph[j] + g) % 3))); X[h].x += t2.x; X[h].y += t2.y;
      }
      int ok = 1; for (int h = 0; h < p; h++) if (X[h].x || X[h].y) { ok = 0; break; }
      if (ok) {
        nsol++;
        if (nsol <= maxprint) {
          printf("SOL g=%d R1:", g); for (int i = 0; i < d; i++) printf(" (%d,c%d,ph%d)", del[0][i], 2 * r1->cls[i], r1->ph[i]);
          printf(" R2:"); for (int i = 0; i < d; i++) printf(" (%d,c%d,ph%d)", del[1][i], 2 * (r2->cls[i] ^ sw), (r2->ph[i] + g) % 3);
          printf("\n"); fflush(stdout);
        }
      }
    }
  }
  printf("p=%d mu=%d flip=%d solutions=%ld\n", p, mu, flip, nsol);
  return 0;
}
