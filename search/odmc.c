/* odmc.c: the list of H-invariant c parts of one Goethals-Seidel block (see odmx.c) with a given b part and row sum,
 * built by a depth-first search over the orbit variables with incremental vectors.
 *
 * Block 1 (skew): one variable per pair of opposite orbits (values v on O, -v on -O), position 0 holds a.
 * Blocks 2-4: one variable per orbit not carrying b.  Each list entry is the int16 vector
 *   (PAF(w)(s_k) for the pair-orbit representatives s_k, X(s_k) = sum_a u(a)[w(a + s_k) + w(a - s_k)])
 * in the representative order of odmx.c, and the mask of the variable values (bit j set: variable j = -1).
 *
 * Usage:  odmc t "h1,h2" block "i:pos:sign,..." rowsum out.lst      (block = 1..4, the entries "block:pos:sign" of
 *                                                                    the template give its b part; block 1 is skew)
 *         odmc t "h1,h2" block "template" rowsum --decode mask      (prints the block as codes)
 * File format: int32 n, int32 D, int16 vec[n][D], uint64 mask[n].
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define MAXT 1024
#define MAXV 64
#define MAXR 40
static int t, nH, Hel[MAXT], orb[MAXT], norb, orep[MAXT], osz[MAXT];
static int prep[MAXR], npr;
static int U[MAXT];
static int blk;                    /* block index 1..4; block 1 is the skew one */
static int nv, vorb[MAXV][2], vcoef[MAXV][2], vlen[MAXV], vsize[MAXV];
static int Dg[MAXV][MAXR], Cx[MAXV][MAXV][MAXR], Lv[MAXV][MAXR];
static int target;
static FILE *out; static long cnt;
static int16_t *vbuf; static uint64_t *mbuf; static long cap;

static void flush_entry(int *P, int *X, uint64_t mask) {
    if (cnt == cap) {
        cap = cap ? 2 * cap : 1 << 16;
        vbuf = realloc(vbuf, sizeof(int16_t) * 2 * npr * cap);
        mbuf = realloc(mbuf, sizeof(uint64_t) * cap);
        if (!vbuf || !mbuf) { fprintf(stderr, "out of memory\n"); exit(1); }
    }
    for (int k = 0; k < npr; k++) { vbuf[cnt * 2 * npr + k] = (int16_t)P[k]; vbuf[cnt * 2 * npr + npr + k] = (int16_t)X[k]; }
    mbuf[cnt] = mask;
    cnt++;
}

static int val[MAXV];
static int suffix_size[MAXV + 1];

static void dfs(int j, int rs, int *P, int *X, uint64_t mask) {
    if (j == nv) { if (rs == target) flush_entry(P, X, mask); return; }
    int rem = target - rs;
    if (rem > suffix_size[j] || rem < -suffix_size[j]) return;
    int P2[MAXR], X2[MAXR];
    for (int s = 0; s < 2; s++) {
        int v = s ? -1 : 1;
        for (int k = 0; k < npr; k++) {
            int d = Dg[j][k];
            for (int i = 0; i < j; i++) d += v * val[i] * Cx[j][i][k];
            P2[k] = P[k] + d;
            X2[k] = X[k] + v * Lv[j][k];
        }
        val[j] = v;
        dfs(j + 1, rs + v * vsize[j], P2, X2, s ? (mask | (1ULL << j)) : mask);
    }
}

int main(int argc, char **argv) {
    if (argc < 7) { fprintf(stderr, "usage: odmc t gens block template rowsum out|--decode mask\n"); return 1; }
    t = atoi(argv[1]);
    int gens[16], ng = 0;
    for (char *tok = strtok(strdup(argv[2]), ","); tok; tok = strtok(0, ",")) gens[ng++] = atoi(tok);
    blk = atoi(argv[3]);
    target = atoi(argv[5]);
    int decode = !strcmp(argv[6], "--decode");
    nH = 1; Hel[0] = 1;
    for (int ch = 1; ch;) {
        ch = 0;
        for (int i = 0; i < nH; i++) for (int j = 0; j < ng; j++) {
            int x = (int)((long)Hel[i] * gens[j] % t), seen = 0;
            for (int k = 0; k < nH; k++) if (Hel[k] == x) seen = 1;
            if (!seen) { Hel[nH++] = x; ch = 1; }
        }
    }
    /* b part of this block: template entries "i:pos:sign" with i == blk (block 2 means the given block's entries are
       passed with i = 2) */
    memset(U, 0, sizeof U);
    for (char *tok = strtok(strdup(argv[4]), ","); tok; tok = strtok(0, ",")) {
        int i, p, s;
        if (sscanf(tok, "%d:%d:%d", &i, &p, &s) != 3) { fprintf(stderr, "bad template\n"); return 1; }
        if (i == blk) U[((p % t) + t) % t] = s;
    }
    for (int g = 0; g < t; g++) orb[g] = -1;
    norb = 0;
    for (int g = 0; g < t; g++) {
        if (orb[g] >= 0) continue;
        orep[norb] = g; osz[norb] = 0;
        for (int k = 0; k < nH; k++) { int x = (int)((long)Hel[k] * g % t); if (orb[x] < 0) { orb[x] = norb; osz[norb]++; } }
        norb++;
    }
    for (int g = 1; g < t; g++) if (orb[g] == orb[(t - g) % t]) { fprintf(stderr, "H not admissible\n"); return 1; }
    /* pair-orbit representatives, same order as odmx.c */
    npr = 0;
    int seenp[MAXT]; memset(seenp, 0, sizeof seenp);
    for (int o = 1; o < norb; o++) {
        int r = orep[o], oo = orb[(t - r) % t];
        if (seenp[o] || seenp[oo]) continue;
        seenp[o] = seenp[oo] = 1;
        if (npr >= MAXR) { fprintf(stderr, "too many pair orbits\n"); return 1; }
        prep[npr++] = r;
    }
    /* variables, same order as odmx.c (orbits in index order; skew block: orbit with the smaller representative) */
    nv = 0;
    for (int o = 0; o < norb; o++) {
        int r = orep[o];
        if (r == 0 && (blk == 1 || U[0])) continue;
        if (U[r]) continue;
        if (blk == 1) {
            int oo = orb[(t - r) % t];
            if (orep[oo] < r) continue;
            vorb[nv][0] = o; vcoef[nv][0] = 1; vorb[nv][1] = oo; vcoef[nv][1] = -1; vlen[nv] = 2; vsize[nv] = 0;
        } else {
            vorb[nv][0] = o; vcoef[nv][0] = 1; vlen[nv] = 1; vsize[nv] = osz[o];
        }
        if (++nv > MAXV) { fprintf(stderr, "too many variables\n"); return 1; }
    }
    if (decode) {
        uint64_t mask = strtoull(argv[7], 0, 10);
        int W[MAXT]; memset(W, 0, sizeof W);
        for (int j = 0; j < nv; j++) {
            int v = (mask >> j) & 1 ? -1 : 1;
            for (int l = 0; l < vlen[j]; l++) for (int g = 0; g < t; g++) if (orb[g] == vorb[j][l]) W[g] = v * vcoef[j][l];
        }
        for (int g = 0; g < t; g++) {
            int code = (blk == 1 && g == 0) ? 1 : U[g] ? 2 * U[g] : 3 * W[g];
            printf("%d%c", code, g == t - 1 ? '\n' : ' ');
        }
        return 0;
    }
    /* T[o1][o2][k] = #{x in o1 : x + s_k in o2}; build per-variable tables */
    static int T[MAXT][MAXT / 4 + 1];   /* not used: compute on the fly below */
    (void)T;
    int *Tk = calloc((size_t)norb * norb * npr, sizeof(int));
    for (int x = 0; x < t; x++) for (int k = 0; k < npr; k++) Tk[((size_t)orb[x] * norb + orb[(x + prep[k]) % t]) * npr + k]++;
    for (int j = 0; j < nv; j++) {
        for (int k = 0; k < npr; k++) {
            int d = 0;
            for (int a = 0; a < vlen[j]; a++) for (int b = 0; b < vlen[j]; b++)
                d += vcoef[j][a] * vcoef[j][b] * Tk[((size_t)vorb[j][a] * norb + vorb[j][b]) * npr + k];
            Dg[j][k] = d;
            int l = 0;
            for (int a = 0; a < vlen[j]; a++)
                for (int g = 0; g < t; g++) if (U[g]) {
                    l += U[g] * vcoef[j][a] * ((orb[(g + prep[k]) % t] == vorb[j][a]) + (orb[(g - prep[k] % t + t) % t] == vorb[j][a]));
                }
            Lv[j][k] = l;
        }
        for (int i = 0; i < j; i++) for (int k = 0; k < npr; k++) {
            int c = 0;
            for (int a = 0; a < vlen[j]; a++) for (int b = 0; b < vlen[i]; b++) {
                int o1 = vorb[j][a], o2 = vorb[i][b];
                c += vcoef[j][a] * vcoef[i][b] * (Tk[((size_t)o1 * norb + o2) * npr + k] + Tk[((size_t)o2 * norb + o1) * npr + k]);
            }
            Cx[j][i][k] = c;
        }
    }
    suffix_size[nv] = 0;
    for (int j = nv - 1; j >= 0; j--) suffix_size[j] = suffix_size[j + 1] + vsize[j];
    int P0[MAXR] = {0}, X0[MAXR] = {0};
    dfs(0, 0, P0, X0, 0);
    out = fopen(argv[6], "wb");
    int32_t hdr[2] = {(int32_t)cnt, 2 * npr};
    fwrite(hdr, 4, 2, out);
    fwrite(vbuf, sizeof(int16_t), (size_t)cnt * 2 * npr, out);
    fwrite(mbuf, sizeof(uint64_t), (size_t)cnt, out);
    fclose(out);
    fprintf(stderr, "block %d rowsum %d: %d variables, %ld entries\n", blk, target, nv, cnt);
    return 0;
}
