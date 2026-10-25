// Bipartite cyclic flag structures for symmetric conference matrices of order 9p + 1 (p = 1 mod 4 prime; search/NOTES.md
// section 96). The flag matrix is Psi = D + 2 Phi on the 4p flags (point, direction of AG(2, 3)); for a structure invariant
// under Z_p with edges only between D = +1 and D = -1 flags, Psi^2 = p I is X X* = d I, d = (p - 1)/4, for a 4 x 4 matrix
// X over Z[w][Z_p] with a real structure. A row (class A1 or A2) has d atoms: an offset in Z_p* (one of each +-pair over
// both rows), a B-column (class B1/B2 and orientation), a cube-root phase; the coefficient is sigma(offset) w^phase,
// sigma = chi * (+-1 by the sign convention flip). Row conditions: same-column self-correlation F and opposite-orientation
// correlation G vanish off zero. Stage 1: rows by a depth-first search over the +-pairs (offset sign, B-class) pruned
// by the number of differences that occur once (each needs a partner), only masks minimal under multiplication by
// squares; orientations and phases by backtracking with each difference checked when its last atom is set. Stage 2:
// R1 with every image mu * R2 of a canonical row whose mask is the complement; cross conditions; B-degree d.
//
//   bipcyc p flip [decoupled]     decoupled = 1: every atom has orientation 0 (then G = 0 identically)
// Output lines `SOL g=.. R1: (offset,cCOLUMN,phPHASE) ... R2: ...` are read by build_bip.py.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef struct { int x, y; } Z;
static inline Z zmul(Z a, Z b) { Z r = { a.x*b.x - a.y*b.y, a.x*b.y + a.y*b.x - a.y*b.y }; return r; }
static inline Z zconj(Z a) { Z r = { a.x - a.y, -a.y }; return r; }
static const Z WP[3] = { {1,0}, {0,1}, {-1,-1} };
int p, d, m, flip; int sig[512];
typedef struct { int16_t del[40]; int8_t col[40]; int8_t ph[40]; uint64_t pairmask; } Row;
Row *rows; long nrows = 0, caprows = 0;
static inline Z kap(int del, int c) { Z w = WP[c]; Z r = { w.x*sig[del], w.y*sig[del] }; return r; }
// per-survivor data
int D_[40], CL[40]; int ORI[40], PH[40];
int npr[40]; int prI[40][1600], prJ[40][1600], prH[40][1600];   // pairs closing at atom k: (i, j, h), i,j <= k
static int closecheck(int k) {
  // recompute F,G for every h that has a pair closing at k: need all pairs of that h (all close at <= k? no)
  // simpler: hlast[h] = max atom index over pairs of h; check h when k == hlast[h]
  return 1;
}
int hlast[512]; int hn[512]; int hi_[512][256], hj_[512][256];
long nsolrow = 0;
static void store(uint64_t pm) {
  if (nrows == caprows) { caprows = caprows ? 2 * caprows : 1 << 16; rows = realloc(rows, caprows * sizeof(Row)); }
  Row *r = &rows[nrows++]; r->pairmask = pm;
  for (int i = 0; i < d; i++) { r->del[i] = D_[i]; r->col[i] = 2 * CL[i] + ORI[i]; r->ph[i] = PH[i]; }
}
int hlist_at[40][512]; int nh_at[40];
int onlyo0 = 0;
static void bt(int k, uint64_t pm) {
  if (k == d) { store(pm); return; }
  for (int o = 0; o < (onlyo0 ? 1 : 2); o++) for (int c = 0; c < 3; c++) {
    if (k == 0 && c) continue;
    ORI[k] = o; PH[k] = c;
    int ok = 1;
    for (int t = 0; t < nh_at[k] && ok; t++) {
      int h = hlist_at[k][t]; Z F = {0,0}, G = {0,0};
      for (int q = 0; q < hn[h]; q++) {
        int i = hi_[h][q], j = hj_[h][q];
        Z ki = kap(D_[i], PH[i]), kj = kap(D_[j], PH[j]);
        if (ORI[i] == ORI[j]) { Z t2 = zmul(ki, zconj(kj)); F.x += t2.x; F.y += t2.y; }
        else { Z t2 = zmul(ki, kj); G.x += t2.x; G.y += t2.y; }
      }
      if (F.x || F.y || G.x || G.y) ok = 0;
    }
    if (ok) bt(k + 1, pm);
  }
}
int sq[512]; int nsq = 0;
static uint64_t mulmask(uint64_t pm, int mu);

static uint64_t mulmask(uint64_t pm, int mu) {
  uint64_t r = 0;
  for (int j = 1; j <= m; j++) if (pm >> (j - 1) & 1) { int x = j * mu % p; if (x > m) x = p - x; r |= 1ull << (x - 1); }
  return r;
}

int cntd[512], n1 = 0; long leaves = 0, survivors = 0;
int selD[40], selC[40];
static inline void addd(int h, int dlt) { int o = cntd[h]; cntd[h] += dlt; if (o == 1) n1--; if (cntd[h] == 1) n1++; }
static void leaf(uint64_t pm) {
  leaves++;
  for (int t = 0; t < nsq; t++) if (mulmask(pm, sq[t]) < pm) return;     // canonical masks only
  survivors++;
  for (int i = 0; i < d; i++) { D_[i] = selD[i]; CL[i] = selC[i]; }
  for (int h = 0; h < p; h++) { hn[h] = 0; hlast[h] = -1; }
  for (int i = 0; i < d; i++) for (int j = 0; j < d; j++) if (i != j && CL[i] == CL[j]) {
    int h = ((D_[i] - D_[j]) % p + p) % p; hi_[h][hn[h]] = i; hj_[h][hn[h]] = j; hn[h]++;
    int mx = i > j ? i : j; if (mx > hlast[h]) hlast[h] = mx;
  }
  for (int t = 0; t < d; t++) nh_at[t] = 0;
  for (int h = 1; h < p; h++) if (hn[h]) { int t = hlast[h]; hlist_at[t][nh_at[t]++] = h; }
  bt(0, pm);
}
static void dfs1(int j, int k, uint64_t pm) {
  int r = d - k;
  if (r == 0) { if (n1 == 0) leaf(pm); return; }
  if (m - j + 1 < r) return;
  if (n1 > 2 * r * k + r * (r - 1)) return;
  for (int sgn = 0; sgn < 2; sgn++) for (int y = 0; y < 2; y++) {
    int del = sgn ? p - j : j;
    for (int i = 0; i < k; i++) if (selC[i] == y) { addd(((del - selD[i]) % p + p) % p, 1); addd(((selD[i] - del) % p + p) % p, 1); }
    selD[k] = del; selC[k] = y;
    dfs1(j + 1, k + 1, pm | (1ull << (j - 1)));
    for (int i = 0; i < k; i++) if (selC[i] == y) { addd(((del - selD[i]) % p + p) % p, -1); addd(((selD[i] - del) % p + p) % p, -1); }
  }
  dfs1(j + 1, k, pm);
}

int main(int argc, char **argv) {
  p = atoi(argv[1]); flip = argc > 2 ? atoi(argv[2]) : 0; onlyo0 = argc > 3 ? atoi(argv[3]) : 0; d = (p - 1) / 4; m = (p - 1) / 2;
  int qr[512] = {0}; for (int x = 1; x < p; x++) qr[x * x % p] = 1;
  for (int x = 1; x < p; x++) sig[x] = (qr[x] ? 1 : -1) * (flip ? -1 : 1);
  for (int x = 1; x < p; x++) if (qr[x]) sq[nsq++] = x;
  long cand = 0, surv = 0, masks = 0;
  memset(cntd, 0, sizeof(cntd)); n1 = 0;
  dfs1(1, 0, 0);
  cand = leaves; surv = survivors;
  fprintf(stderr, "p=%d flip=%d canonical masks=%ld candidates=%ld survivors=%ld canonical rows=%ld\n", p, flip, masks, cand, surv, nrows);
  fflush(stderr);
  // join: R1 canonical row, R2 = mu * (canonical row) for all squares mu, complementary masks
  long nsol = 0; uint64_t full = (m == 64) ? ~0ull : ((1ull << m) - 1);
  // index canonical rows by mask (rows are stored grouped by mask in DFS order; build an index of first/last)
  // sort indices by mask
  long *ord = malloc(nrows * sizeof(long)); for (long i = 0; i < nrows; i++) ord[i] = i;
  int cmpm(const void *x, const void *y) { uint64_t a = rows[*(const long *)x].pairmask, b = rows[*(const long *)y].pairmask; return a < b ? -1 : a > b; }
  qsort(ord, nrows, sizeof(long), cmpm);
  for (long a = 0; a < nrows; a++) {
    Row *r1 = &rows[a];
    uint64_t M2 = full ^ r1->pairmask;
    // canonical form of M2 and all mu with mu * canon = M2
    uint64_t cm = M2; for (int t = 0; t < nsq; t++) { uint64_t x = mulmask(M2, sq[t]); if (x < cm) cm = x; }
    // binary search rows with mask cm
    long lo = 0, hi = nrows;
    while (lo < hi) { long md = (lo + hi) / 2; if (rows[ord[md]].pairmask < cm) lo = md + 1; else hi = md; }
    for (int t = 0; t < nsq; t++) {
      int mu = sq[t];
      if (mulmask(cm, mu) != M2) continue;
      for (long q = lo; q < nrows && rows[ord[q]].pairmask == cm; q++) {
        Row *rb = &rows[ord[q]];
        Row r2; for (int i = 0; i < d; i++) { r2.del[i] = rb->del[i] * mu % p; r2.col[i] = rb->col[i]; r2.ph[i] = rb->ph[i]; }
        int nb1 = 0; for (int i = 0; i < d; i++) { nb1 += (r1->col[i] >> 1) == 0; nb1 += (r2.col[i] >> 1) == 0; }
        if (nb1 != d) continue;
        for (int g = 0; g < 3; g++) {
          static Z X[512], Y[512];
          for (int h = 0; h < p; h++) { X[h].x = X[h].y = 0; Y[h].x = Y[h].y = 0; }
          for (int i = 0; i < d; i++) for (int j = 0; j < d; j++) {
            int h = ((r1->del[i] - r2.del[j]) % p + p) % p;
            Z ki = kap(r1->del[i], r1->ph[i]), kj = kap(r2.del[j], (r2.ph[j] + g) % 3);
            if (r1->col[i] == r2.col[j]) { Z t2 = zmul(ki, zconj(kj)); X[h].x += t2.x; X[h].y += t2.y; }
            if (r2.col[j] == (r1->col[i] ^ 1)) { Z t2 = zmul(ki, kj); Y[h].x += t2.x; Y[h].y += t2.y; }
          }
          int ok = 1; for (int h = 0; h < p; h++) if (X[h].x || X[h].y || Y[h].x || Y[h].y) { ok = 0; break; }
          if (ok) {
            nsol++;
            if (nsol <= 10) {
              printf("SOL g=%d R1:", g); for (int i = 0; i < d; i++) printf(" (%d,c%d,ph%d)", r1->del[i], r1->col[i], r1->ph[i]);
              printf(" R2:"); for (int i = 0; i < d; i++) printf(" (%d,c%d,ph%d)", r2.del[i], r2.col[i], (r2.ph[i] + g) % 3);
              printf("\n"); fflush(stdout);
            }
          }
        }
      }
    }
  }
  printf("p=%d flip=%d canonical-rows=%ld solutions(R1 canonical)=%ld\n", p, flip, nrows, nsol);
  return 0;
}
