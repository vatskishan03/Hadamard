/* odmx.c: exhaustive search for skew-type OD(4t; 1, m, 4t - m - 1) in Goethals-Seidel form, invariant under a
 * multiplier group H of Z_t, with a given b part (see odm_sa.c for the conditions).
 *
 * Blocks: A1 = a e0 + b u1 + c w1 (w1 skew), Ai = b ui + c wi.  The b part must be H-invariant (checked) and satisfy
 * sum_i PAF(ui) = m delta with u1 skew (checked).  Each c part is H-invariant: w1 has one sign per pair of opposite
 * H-orbits, wi one sign per H-orbit.  For every row-sum vector (0, r2, r3, r4) with r2^2 + r3^2 + r4^2 = 4t - m - 1
 * and sum_i rowsum(ui) ri = 0, the four lists of c parts with these row sums are formed; each entry is the vector
 * (PAF(wi)(s), Xi(s)) over representatives s of the pair orbits, Xi(s) = sum_a ui(a) [wi(a + s) + wi(a - s)].
 * A solution is a choice with the four vectors summing to zero (meet in the middle: lists 1 + 2 hashed, 3 + 4
 * looked up).
 *
 * Usage: odmx t "h1,h2,..." "i:pos:sign,..." [maxsol]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define MAXT 512
static int t, nH, Hel[MAXT];
static int U[4][MAXT];
static int orb[MAXT], norb, osz[MAXT], orep[MAXT];
static int prep[MAXT], npr;          /* representatives of pair orbits (s and -s identified) */
static int D;                        /* vector dimension = 2 npr */

typedef struct { int16_t v[64]; int32_t idx; } Ent;
static int *lists[4]; static int lcnt[4];       /* assignments as bitmasks (index into per-block var list) */
static int nvar[4], var_orb[4][MAXT];           /* per block: variable -> orbit id (block 1: orbit of the + side) */

static int W[MAXT];
static void build_w(int i, uint64_t mask) {
    memset(W, 0, sizeof W);
    for (int k = 0; k < nvar[i]; k++) {
        int sg = (mask >> k) & 1 ? -1 : 1, o = var_orb[i][k];
        for (int g = 0; g < t; g++) if (orb[g] == o) W[g] = sg;
        if (i == 0) {     /* skew: the opposite orbit gets -sg */
            int r = orep[o], nr = (t - r) % t, oo = orb[nr];
            for (int g = 0; g < t; g++) if (orb[g] == oo) W[g] = -sg;
        }
    }
}
static int rowsum(void) { int s = 0; for (int g = 0; g < t; g++) s += W[g]; return s; }
static void vec(int i, int16_t *v) {
    for (int k = 0; k < npr; k++) {
        int s = prep[k];
        long p = 0, x = 0;
        for (int g = 0; g < t; g++) p += W[g] * W[(g + s) % t];
        for (int a = 0; a < t; a++) if (U[i][a]) x += U[i][a] * (W[(a + s) % t] + W[(a - s + t) % t]);
        v[k] = (int16_t)p; v[npr + k] = (int16_t)x;
    }
}

static uint64_t hv(const int16_t *v) {
    uint64_t x = 1469598103934665603ULL;
    for (int k = 0; k < D; k++) { x ^= (uint16_t)v[k]; x *= 1099511628211ULL; }
    return x;
}

int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage: odmx t h1,h2,... template [maxsol]\n"); return 1; }
    t = atoi(argv[1]);
    long maxsol = argc > 4 ? atol(argv[4]) : 10;
    int gens[16], ng = 0;
    for (char *tok = strtok(strdup(argv[2]), ","); tok; tok = strtok(0, ",")) gens[ng++] = atoi(tok);
    nH = 1; Hel[0] = 1;
    for (int ch = 1; ch;) {
        ch = 0;
        for (int i = 0; i < nH; i++) for (int j = 0; j < ng; j++) {
            int x = (int)((long)Hel[i] * gens[j] % t), seen = 0;
            for (int k = 0; k < nH; k++) if (Hel[k] == x) seen = 1;
            if (!seen) { Hel[nH++] = x; ch = 1; }
        }
    }
    memset(U, 0, sizeof U);
    for (char *tok = strtok(strdup(argv[3]), ","); tok; tok = strtok(0, ",")) {
        int i, p, s;
        if (sscanf(tok, "%d:%d:%d", &i, &p, &s) != 3) { fprintf(stderr, "bad template\n"); return 1; }
        U[i - 1][((p % t) + t) % t] = s;
    }
    /* orbits of H on Z_t */
    for (int g = 0; g < t; g++) orb[g] = -1;
    norb = 0;
    for (int g = 0; g < t; g++) {
        if (orb[g] >= 0) continue;
        orep[norb] = g; osz[norb] = 0;
        for (int k = 0; k < nH; k++) { int x = (int)((long)Hel[k] * g % t); if (orb[x] < 0) { orb[x] = norb; osz[norb]++; } }
        norb++;
    }
    /* admissibility: no orbit (except {0}) contains both s and -s */
    for (int g = 1; g < t; g++) if (orb[g] == orb[(t - g) % t]) { fprintf(stderr, "H maps %d to -%d: not admissible\n", g, g); return 1; }
    /* template H-invariance and b-part checks */
    int m = 0;
    for (int i = 0; i < 4; i++) for (int g = 0; g < t; g++) {
        if (U[i][g]) m++;
        for (int k = 0; k < nH; k++) if (U[i][(int)((long)Hel[k] * g % t)] != U[i][g]) { fprintf(stderr, "template not H-invariant\n"); return 1; }
    }
    for (int s = 1; s < t; s++) {
        long z = 0;
        for (int i = 0; i < 4; i++) for (int g = 0; g < t; g++) z += U[i][g] * U[i][(g + s) % t];
        if (z) { fprintf(stderr, "b part PAF sum at %d = %ld\n", s, z); return 1; }
    }
    if (U[0][0]) { fprintf(stderr, "u1(0) must be 0\n"); return 1; }
    for (int g = 1; g < t; g++) if (U[0][g] != -U[0][t - g]) { fprintf(stderr, "u1 not skew\n"); return 1; }
    /* pair orbit representatives */
    npr = 0;
    int seenp[MAXT]; memset(seenp, 0, sizeof seenp);
    for (int o = 1; o < norb; o++) {
        int r = orep[o], oo = orb[(t - r) % t];
        if (seenp[o] || seenp[oo]) continue;
        seenp[o] = seenp[oo] = 1;
        prep[npr++] = r;
    }
    D = 2 * npr;
    if (D > 64) { fprintf(stderr, "too many pair orbits (%d)\n", npr); return 1; }
    /* variables per block */
    for (int i = 0; i < 4; i++) {
        nvar[i] = 0;
        for (int o = 0; o < norb; o++) {
            int r = orep[o];
            if (r == 0 && (i == 0 || U[i][0])) continue;   /* position 0: a (block 1) or b */
            if (U[i][r]) continue;                /* b positions */
            if (i == 0) {                          /* skew: one variable per opposite pair, take the orbit with smaller rep */
                int oo = orb[(t - r) % t];
                if (orep[oo] < r) continue;
            }
            var_orb[i][nvar[i]++] = o;
        }
    }
    int ru[4];
    for (int i = 0; i < 4; i++) { ru[i] = 0; for (int g = 0; g < t; g++) ru[i] += U[i][g]; }
    int cw = 4 * t - m - 1;
    fprintf(stderr, "t=%d |H|=%d orbits %d pair-orbits %d m=%d c-weight %d vars %d %d %d %d, b row sums %d %d %d %d\n",
            t, nH, norb, npr, m, cw, nvar[0], nvar[1], nvar[2], nvar[3], ru[0], ru[1], ru[2], ru[3]);
    for (int i = 0; i < 4; i++) if (nvar[i] > 26) { fprintf(stderr, "block %d has %d variables: too many\n", i + 1, nvar[i]); return 1; }
    long found = 0;
    /* free counts for parity */
    int nfree[4];
    for (int i = 0; i < 4; i++) { nfree[i] = 0; for (int g = 0; g < t; g++) nfree[i] += !(U[i][g] || (i == 0 && g == 0)); }
    for (int r2 = -nfree[1]; r2 <= nfree[1]; r2++) for (int r3 = -nfree[2]; r3 <= nfree[2]; r3++) for (int r4 = -nfree[3]; r4 <= nfree[3]; r4++) {
        if (r2 * r2 + r3 * r3 + r4 * r4 != cw) continue;
        if (ru[1] * r2 + ru[2] * r3 + ru[3] * r4) continue;
        if ((nfree[1] - r2) & 1 || (nfree[2] - r3) & 1 || (nfree[3] - r4) & 1) continue;
        int tg[4] = {0, r2, r3, r4};
        Ent *L[4]; int n[4];
        for (int i = 0; i < 4; i++) {
            long cap = 1L << nvar[i];
            L[i] = malloc(sizeof(Ent) * (size_t)(cap < 4000000 ? cap : 4000000));
            n[i] = 0;
            for (long mask = 0; mask < cap; mask++) {
                build_w(i, mask);
                if (rowsum() != tg[i]) continue;
                if (n[i] >= 4000000) { fprintf(stderr, "list %d too long\n", i); return 1; }
                vec(i, L[i][n[i]].v); L[i][n[i]].idx = (int32_t)mask; n[i]++;
            }
        }
        fprintf(stderr, "rows (0, %d, %d, %d): lists %d %d %d %d\n", r2, r3, r4, n[0], n[1], n[2], n[3]);
        if (!n[0] || !n[1] || !n[2] || !n[3]) { for (int i = 0; i < 4; i++) free(L[i]); continue; }
        /* hash sums of lists 0 and 1 */
        long N01 = (long)n[0] * n[1];
        int hb = 1; while ((1L << hb) < 2 * N01) hb++;
        long hsz = 1L << hb;
        int64_t *hk = malloc(sizeof(int64_t) * hsz);
        int32_t *ha = malloc(sizeof(int32_t) * hsz), *hbb = malloc(sizeof(int32_t) * hsz);
        for (long k = 0; k < hsz; k++) ha[k] = -1;
        int16_t sv[64];
        for (int a = 0; a < n[0]; a++) for (int b = 0; b < n[1]; b++) {
            for (int k = 0; k < D; k++) sv[k] = L[0][a].v[k] + L[1][b].v[k];
            uint64_t h = hv(sv) & (hsz - 1);
            while (ha[h] >= 0) h = (h + 1) & (hsz - 1);
            hk[h] = (int64_t)hv(sv); ha[h] = a; hbb[h] = b;
        }
        for (int c = 0; c < n[2] && found < maxsol; c++) for (int d = 0; d < n[3] && found < maxsol; d++) {
            for (int k = 0; k < D; k++) sv[k] = -(L[2][c].v[k] + L[3][d].v[k]);
            uint64_t full = hv(sv), h = full & (hsz - 1);
            while (ha[h] >= 0) {
                if (hk[h] == (int64_t)full) {
                    int a = ha[h], b = hbb[h], ok = 1;
                    for (int k = 0; k < D; k++) if (L[0][a].v[k] + L[1][b].v[k] != sv[k]) { ok = 0; break; }
                    if (ok) {
                        found++;
                        printf("SOL t=%d m=%d rows (0, %d, %d, %d)\n", t, m, r2, r3, r4);
                        int idx[4] = {L[0][a].idx, L[1][b].idx, L[2][c].idx, L[3][d].idx};
                        for (int i = 0; i < 4; i++) {
                            build_w(i, (uint64_t)idx[i]);
                            for (int g = 0; g < t; g++) {
                                int code = (i == 0 && g == 0) ? 1 : U[i][g] ? 2 * U[i][g] : 3 * W[g];
                                printf("%d%c", code, g == t - 1 ? '\n' : ' ');
                            }
                        }
                        fflush(stdout);
                    }
                }
                h = (h + 1) & (hsz - 1);
            }
        }
        free(hk); free(ha); free(hbb);
        for (int i = 0; i < 4; i++) free(L[i]);
        if (found >= maxsol) break;
    }
    fprintf(stderr, "solutions %ld\n", found);
    return 0;
}
