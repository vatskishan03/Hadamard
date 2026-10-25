// Cyclic flag structures for conference matrices of order 9p+1 invariant under a multiplier group K = <mu> of odd order k:
// the block for difference mu^j s is the block for s with phase c + j*tau (tau = 0, or 1, 2 when 3 | k). Enumerates every
// assignment of (u, v, c) to the orbit representatives (general, non-bipartite) and every tau; checks Psi^2 = p I exactly
// (fs3 formalism). usage: mulcyc p k flip [maxprint]
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef struct { int x, y; } Z;
static inline Z zmul(Z a, Z b) { Z r = { a.x*b.x - a.y*b.y, a.x*b.y + a.y*b.x - a.y*b.y }; return r; }
static inline Z zconj(Z a) { Z r = { a.x - a.y, -a.y }; return r; }
static const Z WP[3] = { {1,0}, {0,1}, {-1,-1} };
int p, k, m, d, flip, mu; int sig[1024];
int vec[8][2], negidx[8], cls[8], Dsign[8], posrep[4];
int bu[1024], bv[1024], bc[1024];          // blocks[s], s = 1..m (fs3 convention)
typedef struct { int g, col; Z v; } Ent;
Ent ent[8][4096]; int nent[8];
static void push(int r, int g, int c, Z v) { ent[r][nent[r]].g = ((g % p) + p) % p; ent[r][nent[r]].col = c; ent[r][nent[r]].v = v; nent[r]++; }
static int check(void) {
  for (int a = 0; a < 8; a++) { nent[a] = 0; Z dv = { Dsign[a], 0 }; push(a, 0, a, dv); }
  for (int s = 1; s <= m; s++) {
    int u = bu[s], v = bv[s], c = bc[s];
    Z v1 = { WP[(3 - c) % 3].x * 2 * sig[s], WP[(3 - c) % 3].y * 2 * sig[s] }, v2 = { WP[c].x * 2 * sig[s], WP[c].y * 2 * sig[s] };
    push(u, s, negidx[v], v1); push(negidx[u], s, v, v2); push(negidx[v], -s, u, v2); push(v, -s, negidx[u], v1);
  }
  static Z T[1024];
  for (int a = 0; a < 8; a++) for (int b = a; b < 8; b++) {
    for (int h = 0; h < p; h++) T[h].x = T[h].y = 0;
    for (int i = 0; i < nent[a]; i++) for (int j = 0; j < nent[b]; j++) if (ent[a][i].col == ent[b][j].col) {
      int h = ((ent[a][i].g - ent[b][j].g) % p + p) % p; Z t = zmul(ent[a][i].v, zconj(ent[b][j].v)); T[h].x += t.x; T[h].y += t.y;
    }
    for (int h = 0; h < p; h++) { int tx = (h == 0 && a == b) ? p : 0; if (T[h].x != tx || T[h].y) return 0; }
  }
  return 1;
}
int main(int argc, char **argv) {
  p = atoi(argv[1]); k = atoi(argv[2]); flip = atoi(argv[3]); int maxp = argc > 4 ? atoi(argv[4]) : 5;
  m = (p - 1) / 2; d = (p - 1) / 4;
  int kk = 0; for (int a = 0; a < 3; a++) for (int b = 0; b < 3; b++) if (a || b) { vec[kk][0] = a; vec[kk][1] = b; kk++; }
  for (int i = 0; i < 8; i++) for (int j = 0; j < 8; j++) if ((vec[i][0] + vec[j][0]) % 3 == 0 && (vec[i][1] + vec[j][1]) % 3 == 0) negidx[i] = j;
  int reps[4][2] = { {0,1}, {1,0}, {1,1}, {1,2} };
  for (int c = 0; c < 4; c++) for (int i = 0; i < 8; i++) if (vec[i][0] == reps[c][0] && vec[i][1] == reps[c][1]) { posrep[c] = i; cls[i] = cls[negidx[i]] = c; }
  for (int i = 0; i < 8; i++) { int z = (vec[i][0] == 0) + (vec[i][1] == 0); Dsign[i] = (z == 1) ? 1 : -1; }
  int qr[1024] = {0}; for (int x = 1; x < p; x++) qr[x * x % p] = 1;
  for (int s = 1; s < p; s++) sig[s] = (qr[s] ? -1 : 1) * (flip ? -1 : 1);
  // generator of the subgroup of order k
  int g = 0; for (int c = 2; c < p && !g; c++) { int o = 1, x = c; while (x != 1) { x = x * c % p; o++; } if (o == p - 1) g = c; }
  mu = 1; for (int i = 0; i < (p - 1) / k; i++) mu = mu * g % p;
  // orbits of +-pairs: representatives s (1..m); for each s in 1..m: (rep index, j with s = +-mu^j rep, sign: whether mu^j rep = s or -s)
  int orb[1024], oj[1024], osg[1024], nreps = 0, rep[1024];
  for (int s = 1; s <= m; s++) orb[s] = -1;
  for (int s = 1; s <= m; s++) if (orb[s] < 0) {
    rep[nreps] = s; int x = s;
    for (int j = 0; j < k; j++) { int t = x <= m ? x : p - x; if (orb[t] < 0) { orb[t] = nreps; oj[t] = j; osg[t] = (x <= m) ? 1 : -1; } x = x * mu % p; }
    nreps++;
  }
  int ntau = (k % 3 == 0) ? 3 : 1;
  fprintf(stderr, "p=%d k=%d mu=%d orbits of +-pairs=%d tau options=%d\n", p, k, mu, nreps, ntau);
  if (nreps > 4) { fprintf(stderr, "too many orbits for brute force\n"); return 1; }
  long long tot = 1; for (int i = 0; i < nreps; i++) tot *= 96; tot *= ntau;
  long long sols = 0;
  for (long long code = 0; code < tot; code++) {
    long long t = code; int tau = t % ntau; t /= ntau;
    int U[8], V[8], Cc[8];
    for (int i = 0; i < nreps; i++) { int o = t % 96; t /= 96; U[i] = posrep[o / 24]; V[i] = (o / 3) % 8; Cc[i] = o % 3; }
    for (int s = 1; s <= m; s++) {
      int r = orb[s], j = oj[s]; int c = (Cc[r] + j * tau) % 3;
      if (osg[s] > 0) { bu[s] = U[r]; bv[s] = V[r]; bc[s] = c; }
      else { bu[s] = V[r]; bv[s] = U[r]; bc[s] = c; }   // transpose: hyper(v, u, c)
      // canonical (u in posrep): (u, v, c) ~ (-u, -v, -c)
      if (cls[bu[s]] >= 0 && bu[s] != posrep[cls[bu[s]]]) { bu[s] = negidx[bu[s]]; bv[s] = negidx[bv[s]]; bc[s] = (3 - bc[s]) % 3; }
    }
    if (check()) {
      sols++;
      if (sols <= maxp) { printf("SOL:"); for (int s = 1; s <= m; s++) printf(" [%d:(%d,%d)(%d,%d)c%d]", s, vec[bu[s]][0], vec[bu[s]][1], vec[bv[s]][0], vec[bv[s]][1], bc[s]); printf("\n"); fflush(stdout); }
    }
  }
  printf("p=%d k=%d flip=%d solutions=%lld (of %lld)\n", p, k, flip, sols, tot);
  return 0;
}
