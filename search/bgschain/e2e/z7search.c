// Exhaustive search for Z_7-Singer-invariant borders at k = 8 (Q = Z_2^3), spike-at-0 compressions with amplitudes d.
// Generators e (32 entries: block i, slot 0 = fixed position (eta_i), slots 1..7 = cyclic positions): e is a Legendre quadruple:
//   blocks eps_i in {+-1}^7 with sum = +-1, eta_i = -sum, sum_i PAF_{eps_i}(u) = -4 (u = 1..3).
// E = rows {T^a e_j}, F = columns {T^a f_j'} (j, j' = 1..4).  Conditions:
//   orth: <e, T^u e'> = 0 for all u (e != e'), same for f's;  compat: |f^T T^u (Gamma (x) I) e| = 32 for all u.
// usage: z7search d1 d2 d3 d4 [maxsol]
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

typedef struct { int8_t v[32]; } vec;
static int d[4];
static int G[4][4];
static int S70[70][7]; static int nS = 0;
static int PAF[70][3];
static vec *C; static int nC = 0;
static int8_t (*W)[7][32]; // W[e][u][k] = (T^u (Gamma x I) e)[k]
static int *mm;            // minimal multiplier-variant index
static long long nsol = 0; static long long maxsol = 1;

static void shiftv(const int8_t *v, int u, int8_t *out) {
    for (int i = 0; i < 4; i++) { out[i * 8] = v[i * 8]; for (int b = 0; b < 7; b++) out[i * 8 + 1 + b] = v[i * 8 + 1 + ((b + u) % 7)]; }
}
static void multv(const int8_t *v, int mu, int8_t *out) { // position t -> mu*t
    for (int i = 0; i < 4; i++) { out[i * 8] = v[i * 8]; for (int t = 0; t < 7; t++) out[i * 8 + 1 + (mu * t) % 7] = v[i * 8 + 1 + t]; }
}
static int cmpv(const int8_t *a, const int8_t *b) { return memcmp(a, b, 32); }
static void canon_shift_sign(const int8_t *v, int8_t *best) { // lexicographically smallest of 14 variants (as unsigned bytes: map -1->0,1->1 for memcmp)
    int8_t tmp[32], t2[32]; int first = 1;
    for (int a = 0; a < 7; a++) { shiftv(v, a, tmp); for (int s = 0; s < 2; s++) {
        for (int k = 0; k < 32; k++) t2[k] = s ? -tmp[k] : tmp[k];
        // compare with +1 > -1: use (x+1) as unsigned
        int8_t u1[32], u2[32]; for (int k = 0; k < 32; k++) u1[k] = t2[k] + 1;
        if (first) { memcpy(best, t2, 32); first = 0; }
        else { for (int k = 0; k < 32; k++) u2[k] = best[k] + 1; if (memcmp(u1, u2, 32) < 0) memcpy(best, t2, 32); }
    } }
}
static int is_canon(const int8_t *v) { int8_t b[32]; canon_shift_sign(v, b); return memcmp(b, v, 32) == 0; }

static int cmp_vec(const void *a, const void *b) { // sort by (x+1) bytes
    const int8_t *x = ((const vec *)a)->v, *y = ((const vec *)b)->v;
    for (int k = 0; k < 32; k++) { if (x[k] != y[k]) return x[k] < y[k] ? -1 : 1; }
    return 0;
}

static int orth(const int8_t *e, const int8_t *g) {
    for (int u = 0; u < 7; u++) {
        int8_t t[32]; shiftv(g, u, t); int s = 0; for (int k = 0; k < 32; k++) s += e[k] * t[k];
        if (s != 0) return 0;
    }
    return 1;
}
static int compat(const int8_t *f, int ei) {
    for (int u = 0; u < 7; u++) { int s = 0; const int8_t *w = W[ei][u]; for (int k = 0; k < 32; k++) s += f[k] * w[k]; if (s != 32 && s != -32) return 0; }
    return 1;
}

static int E[4]; static int firstE[4], firstF[4]; static int have = 0;
static int Fsel[4];

static void fsearch(int depth, int *list, int n) {
    if (depth == 4) { nsol++; if (!have) { memcpy(firstE, E, sizeof(E)); memcpy(firstF, Fsel, sizeof(Fsel)); have = 1; } return; }
    for (int a = 0; a < n; a++) {
        int fi = list[a]; Fsel[depth] = fi;
        int *nl = malloc(sizeof(int) * (n + 1)); int m = 0;
        for (int b = a + 1; b < n; b++) if (orth(C[fi].v, C[list[b]].v)) nl[m++] = list[b];
        if (m >= 3 - depth) fsearch(depth + 1, nl, m);
        free(nl);
        if (nsol >= maxsol) return;
    }
}

static long long nE1 = 0, nE1ok = 0, nE2 = 0, nE3 = 0, nE4 = 0;
static void esearch(int depth, int *elist, int ne, int *flist, int nf, int e1) {
    // depth = number of E gens chosen so far (>=1)
    if (depth == 4) { nE4++; fsearch(0, flist, nf); return; }
    for (int a = 0; a < ne; a++) {
        int ei = elist[a];
        // new f list
        int *nfl = malloc(sizeof(int) * (nf + 1)); int m = 0;
        for (int b = 0; b < nf; b++) if (compat(C[flist[b]].v, ei)) nfl[m++] = flist[b];
        if (m < 4) { free(nfl); continue; }
        E[depth] = ei;
        int *nel = malloc(sizeof(int) * (ne + 1)); int me = 0;
        for (int b = a + 1; b < ne; b++) if (orth(C[ei].v, C[elist[b]].v)) nel[me++] = elist[b];
        if (depth == 1) nE2++; else if (depth == 2) nE3++;
        if (me >= 3 - depth) esearch(depth + 1, nel, me, nfl, m, e1);
        free(nel); free(nfl);
        if (nsol >= maxsol) return;
    }
}

int main(int argc, char **argv) {
    for (int i = 0; i < 4; i++) d[i] = atoi(argv[1 + i]);
    if (argc > 5) maxsol = atoll(argv[5]);
    // GS array gs4(d)
    int g[4][4] = {{d[0], d[1], d[2], d[3]}, {-d[1], d[0], d[3], -d[2]}, {-d[2], -d[3], d[0], d[1]}, {-d[3], d[2], -d[1], d[0]}};
    memcpy(G, g, sizeof(g));
    // sequences with |sum| = 1
    for (int x = 0; x < 128; x++) { int s[7], sum = 0; for (int t = 0; t < 7; t++) { s[t] = (x >> t) & 1 ? -1 : 1; sum += s[t]; } if (sum == 1 || sum == -1) { memcpy(S70[nS], s, sizeof(s)); nS++; } }
    for (int a = 0; a < nS; a++) for (int u = 1; u <= 3; u++) { int p = 0; for (int t = 0; t < 7; t++) p += S70[a][t] * S70[a][(t + u) % 7]; PAF[a][u - 1] = p; }
    C = malloc(sizeof(vec) * 3000000);
    for (int a = 0; a < nS; a++) for (int b = 0; b < nS; b++) {
        int p1[3]; for (int u = 0; u < 3; u++) p1[u] = PAF[a][u] + PAF[b][u];
        for (int c = 0; c < nS; c++) {
            int p2[3]; for (int u = 0; u < 3; u++) p2[u] = p1[u] + PAF[c][u];
            for (int e = 0; e < nS; e++) {
                if (p2[0] + PAF[e][0] != -4 || p2[1] + PAF[e][1] != -4 || p2[2] + PAF[e][2] != -4) continue;
                int idx[4] = {a, b, c, e}; vec v;
                for (int i = 0; i < 4; i++) { int sum = 0; for (int t = 0; t < 7; t++) { v.v[i * 8 + 1 + t] = S70[idx[i]][t]; sum += S70[idx[i]][t]; } v.v[i * 8] = -sum; }
                if (!is_canon(v.v)) continue;
                C[nC++] = v;
            }
        }
    }
    fprintf(stderr, "canonical Legendre-quadruple generators: %d\n", nC);
    qsort(C, nC, sizeof(vec), cmp_vec);
    // W images
    W = malloc((size_t)nC * 7 * 32);
    for (int e = 0; e < nC; e++) {
        int8_t w[32];
        for (int j = 0; j < 4; j++) for (int p = 0; p < 8; p++) { int s = 0; for (int i = 0; i < 4; i++) s += G[j][i] * C[e].v[i * 8 + p]; w[j * 8 + p] = (int8_t)s; }
        for (int u = 0; u < 7; u++) shiftv(w, u, W[e][u]);
    }
    // multiplier-minimal index
    mm = malloc(sizeof(int) * nC);
    for (int e = 0; e < nC; e++) {
        int best = e;
        for (int mu = 2; mu < 7; mu++) {
            int8_t t[32], c[32]; multv(C[e].v, mu, t); canon_shift_sign(t, c);
            vec key; memcpy(key.v, c, 32);
            vec *r = bsearch(&key, C, nC, sizeof(vec), cmp_vec);
            if (r) { int id = (int)(r - C); if (id < best) best = id; } else { fprintf(stderr, "multiplier image not in list!\n"); }
        }
        mm[e] = best;
    }
    long long ncanon = 0; for (int e = 0; e < nC; e++) if (mm[e] == e) ncanon++;
    fprintf(stderr, "multiplier-canonical e1 candidates: %lld\n", ncanon);
    int *all = malloc(sizeof(int) * nC); for (int i = 0; i < nC; i++) all[i] = i;
    for (int e1 = 0; e1 < nC; e1++) {
        if (mm[e1] != e1) continue;
        nE1++;
        // F list compatible with e1
        int *fl = malloc(sizeof(int) * (nC / 4 + 16)); int nf = 0; int cap = nC / 4 + 16;
        for (int f = 0; f < nC; f++) if (compat(C[f].v, e1)) { if (nf < cap) fl[nf++] = f; }
        if (nf < 4) { free(fl); continue; }
        nE1ok++;
        // E candidates: orth to e1, mm(e) > e1 (or >= and distinct), other indices
        int *el = malloc(sizeof(int) * (nC / 8 + 16)); int ne = 0; int capE = nC / 8 + 16;
        for (int e = 0; e < nC; e++) { if (e == e1 || mm[e] < e1) continue; if (orth(C[e1].v, C[e].v)) { if (ne < capE) el[ne++] = e; } }
        E[0] = e1;
        if (ne >= 3) esearch(1, el, ne, fl, nf, e1);
        free(el); free(fl);
        if (nE1 % 500 == 0) fprintf(stderr, "e1 %d/%d: with>=4 compat F: %lld, E2 nodes %lld, E3 %lld, E4 %lld, solutions %lld\n", e1, nC, nE1ok, nE2, nE3, nE4, nsol);
        if (nsol >= maxsol) break;
    }
    printf("d=(%d,%d,%d,%d): e1 candidates %lld, with >=4 compat F %lld, E2 %lld E3 %lld E4 %lld, solutions %lld\n", d[0], d[1], d[2], d[3], nE1, nE1ok, nE2, nE3, nE4, nsol);
    if (have) {
        printf("SOL\n");
        for (int j = 0; j < 4; j++) { for (int k = 0; k < 32; k++) printf("%d ", C[firstE[j]].v[k]); printf("\n"); }
        for (int j = 0; j < 4; j++) { for (int k = 0; k < 32; k++) printf("%d ", C[firstF[j]].v[k]); printf("\n"); }
    }
    return 0;
}
