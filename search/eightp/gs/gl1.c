// gl1.c : level 1 of the REAL Goethals-Seidel-type plug-ins (four FREE blocks), any abelian K (shift tables), any odd r  (search/NOTES.md sections 86 and 92).
// Each block is a list of (s,e) configs (defect weight nw exactly), either EXPLICIT (arrays S,E) or GENERATED (uniform over distinct (s,e) configs of the given weight,
// regenerated deterministically from (seed, index); optional exact row-sum filter R = sum u in rok (rejection sampling), optional canonical block sign).
// The level-1 system is  sum_X PAF_u(k) = target(k) for k in the reps (u = s + r e); join (blk0 x blk2) against (blk1 x blk3), class-restricted hash join,
// blocks 0,1 are the OUTER (small) lists, 2,3 the INNER (large) lists.  Every output instance is checked exactly on all features.
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#define MAXN 640
static inline uint64_t sm64(uint64_t *x) { uint64_t z = (*x += 0x9E3779B97F4A7C15ULL); z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL; z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL; return z ^ (z >> 31); }
static inline uint64_t xs(uint64_t *s) { uint64_t x = *s; x ^= x >> 12; x ^= x << 25; x ^= x >> 27; *s = x; return x * 0x2545F4914F6CDD1DULL; }

typedef struct {
    int n, r, nrep; const int32_t *shift;
    int nu; const int *w, *forced, *poff, *plist; const double *dpw; int nfree; const int *freeidx; int fw;
    uint64_t coef[96];
} Ctx;
typedef struct { int nw; long L; uint64_t seed; const int8_t *S, *E; const uint8_t *rok; int RM; int cansign; int mw; const uint8_t *sok, *eok; int SM, EM; int fs, fe; } Blk;

#define NONE_FIX (-100000)
static double TP_[MAXN + 1][2 * MAXN + 2], TM_[MAXN + 1][2 * MAXN + 2], TD_[MAXN + 1][2 * MAXN + 2];
static void subset_table(const int *wt, int nit, double (*T)[2 * MAXN + 2], int W) {
    for (int s = 0; s <= W; s++) T[nit][s] = (s == 0) ? 1.0 : 0.0;
    for (int j = nit - 1; j >= 0; j--) for (int s = 0; s <= W; s++) T[j][s] = T[j + 1][s] + (s >= wt[j] ? T[j + 1][s - wt[j]] : 0.0);
}
// uniformly random subset of items with total weight s (table T); plus[j] = 1 if chosen
static void subset_sample(const int *wt, int nit, double (*T)[2 * MAXN + 2], int s, uint64_t *rs, uint8_t *plus) {
    for (int j = 0; j < nit; j++) {
        double tk = s >= wt[j] ? T[j + 1][s - wt[j]] : 0.0, tot = T[j][s];
        double u = (double)(xs(rs) >> 11) * (1.0 / 9007199254740992.0);
        if (u * tot < tk) { plus[j] = 1; s -= wt[j]; } else plus[j] = 0;
    }
}

static int gen_one(const Ctx *c, int nw, uint64_t *rs, int8_t *s, int8_t *e, int cansign, int mw, int fs, int fe) {
    int n = c->n; int stride = n + 1;
    memset(s, 0, n); memset(e, 0, n);
    uint8_t defect[MAXN]; memset(defect, 0, sizeof(defect));
    for (int i = 0; i < c->nu; i++) if (c->forced[i]) defect[i] = 1;
    int rem = nw - c->fw;
    for (int ii = 0; ii < c->nfree; ii++) {
        int i = c->freeidx[ii]; int wi = c->w[i];
        double skip = 2.0 * c->dpw[(ii + 1) * stride + rem];
        double take = rem >= wi ? 4.0 * c->dpw[(ii + 1) * stride + rem - wi] : 0.0;
        double p = take / (take + skip + 1e-300);
        double u = (double)(xs(rs) >> 11) * (1.0 / 9007199254740992.0);
        if (u < p) { defect[i] = 1; rem -= wi; }
    }
    // m-type defects (s = -2e) of total orbit weight exactly mw (mw < 0: types uniform i.i.d.)
    uint8_t mtype[MAXN]; memset(mtype, 0, sizeof(mtype));
    int dl[MAXN], nd = 0; for (int i = 0; i < c->nu; i++) if (defect[i]) dl[nd++] = i;
    if (mw >= 0) {
        static double cnt[MAXN + 1][MAXN * 4 + 2];
        int W = mw;
        for (int w = 0; w <= W; w++) cnt[nd][w] = (w == 0) ? 1.0 : 0.0;
        for (int j = nd - 1; j >= 0; j--) { int wj = c->w[dl[j]]; for (int w = 0; w <= W; w++) cnt[j][w] = cnt[j + 1][w] + (w >= wj ? cnt[j + 1][w - wj] : 0.0); }
        if (cnt[0][W] <= 0) return 0;
        int rm = W;
        for (int j = 0; j < nd; j++) { int wj = c->w[dl[j]]; double tk = rm >= wj ? cnt[j + 1][rm - wj] : 0.0; double tot = cnt[j][rm];
            double u = (double)(xs(rs) >> 11) * (1.0 / 9007199254740992.0);
            if (u * tot < tk) { mtype[dl[j]] = 1; rm -= wj; } }
    } else {
        for (int j = 0; j < nd; j++) mtype[dl[j]] = (xs(rs) & 1) ? 1 : 0;
    }
    int sg = 0;
    if (fs == NONE_FIX || fe == NONE_FIX) {
        for (int i = 0; i < c->nu; i++) {
            uint64_t b = xs(rs); int be = (b & 1) ? 1 : -1;
            int sv, ev;
            if (defect[i]) { ev = be; sv = mtype[i] ? -2 * be : 0; } else { ev = 0; sv = be; }
            for (int p = c->poff[i]; p < c->poff[i + 1]; p++) { s[c->plist[p]] = (int8_t)sv; e[c->plist[p]] = (int8_t)ev; }
            if (sg == 0) { int v = sv != 0 ? sv : ev; sg = v; }
        }
    } else {
        // exact class sampler: sum s = fs, sum e = fe.   a = sum_{p-defects} e w, b = sum_{m-defects} e w, a + b = fe; sum_{nondef} s w = fs + 2 b
        int Pw[MAXN], Mw[MAXN], Dw[MAXN], Pi[MAXN], Mi[MAXN], Di[MAXN], np_ = 0, nm_ = 0, ndd = 0;
        for (int i = 0; i < c->nu; i++) {
            if (defect[i]) { if (mtype[i]) { Mw[nm_] = c->w[i]; Mi[nm_++] = i; } else { Pw[np_] = c->w[i]; Pi[np_++] = i; } }
            else { Dw[ndd] = c->w[i]; Di[ndd++] = i; }
        }
        int WP = 0, WM = 0, WD = 0; for (int j = 0; j < np_; j++) WP += Pw[j]; for (int j = 0; j < nm_; j++) WM += Mw[j]; for (int j = 0; j < ndd; j++) WD += Dw[j];
        subset_table(Pw, np_, TP_, WP); subset_table(Mw, nm_, TM_, WM); subset_table(Dw, ndd, TD_, WD);
        double wts[2 * MAXN + 2]; int bv[2 * MAXN + 2]; int nb = 0; double tot = 0;
        for (int b = -WM; b <= WM; b += 2) {
            int a = fe - b, sgd = fs + 2 * b;
            if (((WP + a) & 1) || ((WM + b) & 1) || ((WD + sgd) & 1)) continue;
            int pp = (WP + a) / 2, pm = (WM + b) / 2, pd = (WD + sgd) / 2;
            if (pp < 0 || pp > WP || pm < 0 || pm > WM || pd < 0 || pd > WD) continue;
            double w = TP_[0][pp] * TM_[0][pm] * TD_[0][pd];
            if (w > 0) { wts[nb] = w; bv[nb++] = b; tot += w; }
        }
        if (nb == 0) return 0;
        double u = (double)(xs(rs) >> 11) * (1.0 / 9007199254740992.0) * tot; int pick = nb - 1;
        for (int k = 0; k < nb; k++) { if (u < wts[k]) { pick = k; break; } u -= wts[k]; }
        int b = bv[pick], a = fe - b, sgd = fs + 2 * b;
        uint8_t plusP[MAXN], plusM[MAXN], plusD[MAXN];
        subset_sample(Pw, np_, TP_, (WP + a) / 2, rs, plusP); subset_sample(Mw, nm_, TM_, (WM + b) / 2, rs, plusM); subset_sample(Dw, ndd, TD_, (WD + sgd) / 2, rs, plusD);
        for (int j = 0; j < np_; j++) { int i = Pi[j]; int ev = plusP[j] ? 1 : -1; for (int p = c->poff[i]; p < c->poff[i + 1]; p++) { s[c->plist[p]] = 0; e[c->plist[p]] = (int8_t)ev; } }
        for (int j = 0; j < nm_; j++) { int i = Mi[j]; int ev = plusM[j] ? 1 : -1; for (int p = c->poff[i]; p < c->poff[i + 1]; p++) { s[c->plist[p]] = (int8_t)(-2 * ev); e[c->plist[p]] = (int8_t)ev; } }
        for (int j = 0; j < ndd; j++) { int i = Di[j]; int sv = plusD[j] ? 1 : -1; for (int p = c->poff[i]; p < c->poff[i + 1]; p++) { s[c->plist[p]] = (int8_t)sv; e[c->plist[p]] = 0; } }
        for (int i = 0; i < c->nu && sg == 0; i++) { int z = c->plist[c->poff[i]]; int v = s[z] != 0 ? s[z] : e[z]; if (v) sg = v; }
    }
    if (cansign && sg < 0) for (int z = 0; z < n; z++) { s[z] = -s[z]; e[z] = -e[z]; }
    return 1;
}
static void gen_cfg(const Ctx *c, const Blk *b, uint64_t j, int8_t *s, int8_t *e) {
    int n = c->n;
    if (b->S) { memcpy(s, b->S + (size_t)j * n, n); memcpy(e, b->E + (size_t)j * n, n); return; }
    uint64_t st = b->seed ^ (j * 0xD1B54A32D192ED03ULL); uint64_t r0 = sm64(&st); uint64_t rs = r0 | 1;
    for (int tries = 0; tries < 20000; tries++) {
        if (!gen_one(c, b->nw, &rs, s, e, b->cansign, b->mw, b->fs, b->fe)) continue;
        if (!b->rok && !b->sok && !b->eok) return;
        long sg_ = 0, eg_ = 0; for (int z = 0; z < n; z++) { sg_ += s[z]; eg_ += e[z]; }
        long R = sg_ + (long)c->r * eg_;
        if (b->rok && !(R >= -b->RM && R <= b->RM && b->rok[R + b->RM])) continue;
        if (b->sok && !(sg_ >= -b->SM && sg_ <= b->SM && b->sok[sg_ + b->SM])) continue;
        if (b->eok && !(eg_ >= -b->EM && eg_ <= b->EM && b->eok[eg_ + b->EM])) continue;
        return;
    }
}

static void feats(const Ctx *c, const int8_t *s, const int8_t *e, int32_t *f) {
    int n = c->n, r = c->r; int u[MAXN];
    for (int z = 0; z < n; z++) u[z] = s[z] + r * e[z];
    for (int i = 0; i < c->nrep; i++) { const int32_t *sh = c->shift + (size_t)i * n; int a = 0; for (int z = 0; z < n; z++) a += u[z] * u[sh[z]]; f[i] = a; }
}
static inline uint64_t hsh(const Ctx *c, const int32_t *f) { uint64_t h = 0; for (int i = 0; i < c->nrep; i++) h += c->coef[i] * (uint64_t)(int64_t)f[i]; return h; }

typedef struct { uint64_t key; uint32_t a, c; } ent;
typedef struct { uint64_t ch; uint32_t j; } cj_t;
static int cj_cmp(const void *a, const void *bb) { const cj_t *x = a, *y = bb; if (x->ch != y->ch) return x->ch < y->ch ? -1 : 1; return x->j < y->j ? -1 : (x->j > y->j); }

// stats[0] = pairs enumerated on side 0, stats[1] = candidates probed, stats[2] = hash-matching candidates
long gl1_run(int n, int r, int nrep, const int32_t *shift, const int32_t *target,
             int nu, const int *w, const int *forced, const int *poff, const int *plist, const double *dpw, int nfree, const int *freeidx, int fw,
             const int *nw, const int *mwv, const long *L, const uint64_t *seeds, const int8_t **SS, const int8_t **EE, const uint8_t **rok, int RM, const int *cansign, const uint8_t **sok, const uint8_t **eok, int SM, int EM, const int *fsv, const int *fev,
             int lgM, long maxout, int8_t *oS, int8_t *oE, double *stats, uint64_t hseed, long stop_after)
{
    Ctx cx = { n, r, nrep, shift, nu, w, forced, poff, plist, dpw, nfree, freeidx, fw, {0} };
    uint64_t st = hseed * 6364136223846793005ULL + 1; for (int i = 0; i < 96; i++) cx.coef[i] = sm64(&st) | 1;
    Blk B[4];
    for (int b = 0; b < 4; b++) { B[b].nw = nw[b]; B[b].L = L[b]; B[b].seed = seeds[b]; B[b].S = SS[b]; B[b].E = EE[b]; B[b].rok = rok[b]; B[b].RM = RM; B[b].cansign = cansign[b]; B[b].mw = mwv[b]; B[b].sok = sok[b]; B[b].eok = eok[b]; B[b].SM = SM; B[b].EM = EM; B[b].fs = fsv[b]; B[b].fe = fev[b]; }
    uint64_t M = 1ULL << lgM, mm = M - 1;
    int32_t f[96];
    uint64_t *h[4]; uint32_t *map_[4]; long Lw[4];
    int8_t s_[MAXN], e_[MAXN];
    for (int b = 0; b < 4; b++) {
        h[b] = malloc(8L * (L[b] + 1)); map_[b] = NULL; Lw[b] = L[b];
        if (B[b].S) {
            for (long j = 0; j < L[b]; j++) { gen_cfg(&cx, &B[b], j, s_, e_); feats(&cx, s_, e_, f); h[b][j] = hsh(&cx, f); }
        } else {
            // generated list: drop duplicate configs (identical (s,e)); keep the first occurrence of each
            cj_t *cj = malloc(sizeof(cj_t) * (L[b] + 1)); uint64_t *hh = malloc(8L * (L[b] + 1));
            for (long j = 0; j < L[b]; j++) {
                gen_cfg(&cx, &B[b], j, s_, e_); feats(&cx, s_, e_, f); hh[j] = hsh(&cx, f);
                uint64_t c = 0x9E3779B97F4A7C15ULL; for (int z = 0; z < n; z++) { c = (c ^ (uint64_t)(uint8_t)s_[z]) * 0x100000001B3ULL; c ^= c >> 29; c = (c ^ (uint64_t)(uint8_t)e_[z]) * 0x100000001B3ULL; c ^= c >> 31; }
                cj[j].ch = c; cj[j].j = (uint32_t)j;
            }
            qsort(cj, L[b], sizeof(cj_t), cj_cmp);
            map_[b] = malloc(4L * (L[b] + 1)); long m = 0;
            for (long k = 0; k < L[b]; k++) if (k == 0 || cj[k].ch != cj[k - 1].ch) { map_[b][m] = cj[k].j; h[b][m] = hh[cj[k].j]; m++; }
            Lw[b] = m; free(cj); free(hh);
        }
    }
    uint64_t ht = hsh(&cx, target);
    // inner lists (blocks 2,3): sort by class
    long *bs[2]; uint64_t *sh_[2]; uint32_t *io[2];
    for (int m = 0; m < 2; m++) {
        int b = 2 + m; long Lb = Lw[b];
        bs[m] = calloc(M + 1, sizeof(long));
        for (long i = 0; i < Lb; i++) bs[m][(h[b][i] & mm) + 1]++;
        for (uint64_t q = 0; q < M; q++) bs[m][q + 1] += bs[m][q];
        sh_[m] = malloc(8L * (Lb + 1)); io[m] = malloc(4L * (Lb + 1));
        long *fl = malloc(sizeof(long) * M); memcpy(fl, bs[m], sizeof(long) * M);
        for (long i = 0; i < Lb; i++) { long p = fl[h[b][i] & mm]++; sh_[m][p] = h[b][i]; io[m][p] = (uint32_t)i; }
        free(fl); free(h[b]);
    }
    long found = 0; size_t cap = 1 << 10; ent *tab = malloc(sizeof(ent) * cap);
    size_t bmcap = 1 << 10; uint64_t *bm = malloc(8 * bmcap);
    int8_t s4[4][MAXN], e4[4][MAXN]; int32_t fs[4][96];
    double st0 = 0, st1 = 0, st2 = 0;
    for (uint64_t cls = 0; cls < M; cls++) {
        long npairs = 0;
        for (long a = 0; a < Lw[0]; a++) { uint64_t need = (cls - h[0][a]) & mm; npairs += bs[0][need + 1] - bs[0][need]; }
        st0 += npairs;
        if (!npairs) continue;
        size_t nc = 1; while (nc < (size_t)npairs * 2) nc <<= 1;
        if (nc > cap) { cap = nc; tab = realloc(tab, sizeof(ent) * cap); }
        uint64_t cm = nc - 1;
        for (size_t i = 0; i < nc; i++) tab[i].a = 0xffffffffu;
        int bmb = 4; while (((size_t)1 << bmb) < (size_t)npairs * 16 && bmb < 28) bmb++;
        size_t bmw = ((size_t)1 << bmb) >> 6; if (bmw < 1) bmw = 1; if (bmw > bmcap) { bmcap = bmw; bm = realloc(bm, 8 * bmcap); }
        memset(bm, 0, 8 * bmw); uint64_t bmm = ((uint64_t)1 << bmb) - 1;
        for (long a = 0; a < Lw[0]; a++) {
            uint64_t need = (cls - h[0][a]) & mm; uint64_t ha = h[0][a];
            for (long p = bs[0][need]; p < bs[0][need + 1]; p++) {
                uint64_t key = ha + sh_[0][p]; uint64_t i = ((key * 0x9E3779B97F4A7C15ULL) >> 24) & cm;
                { uint64_t bi = (key * 0xC2B2AE3D27D4EB4FULL) >> 40 & bmm; bm[bi >> 6] |= 1ULL << (bi & 63); }
                while (tab[i].a != 0xffffffffu) i = (i + 1) & cm;
                tab[i].key = key; tab[i].a = (uint32_t)a; tab[i].c = io[0][p];
            }
        }
        uint64_t rt = (ht - cls) & mm;
        for (long b = 0; b < Lw[1]; b++) {
            uint64_t need = (rt - h[1][b]) & mm; uint64_t hb = ht - h[1][b];
            for (long p = bs[1][need]; p < bs[1][need + 1]; p++) {
                uint64_t key = hb - sh_[1][p]; st1 += 1;
                { uint64_t bi = (key * 0xC2B2AE3D27D4EB4FULL) >> 40 & bmm; if (!((bm[bi >> 6] >> (bi & 63)) & 1)) continue; }
                uint64_t i = ((key * 0x9E3779B97F4A7C15ULL) >> 24) & cm;
                while (tab[i].a != 0xffffffffu) {
                    if (tab[i].key == key) {
                        st2 += 1;
                        uint32_t a = tab[i].a, cj = tab[i].c, dj = io[1][p];
                        gen_cfg(&cx, &B[0], map_[0] ? map_[0][a] : a, s4[0], e4[0]); gen_cfg(&cx, &B[1], map_[1] ? map_[1][b] : (uint64_t)b, s4[1], e4[1]);
                        gen_cfg(&cx, &B[2], map_[2] ? map_[2][cj] : cj, s4[2], e4[2]); gen_cfg(&cx, &B[3], map_[3] ? map_[3][dj] : dj, s4[3], e4[3]);
                        for (int k = 0; k < 4; k++) feats(&cx, s4[k], e4[k], fs[k]);
                        int ok = 1; for (int k = 0; k < nrep; k++) if (fs[0][k] + fs[1][k] + fs[2][k] + fs[3][k] != target[k]) { ok = 0; break; }
                        if (ok) {
                            if (found < maxout) for (int k = 0; k < 4; k++) { memcpy(oS + ((size_t)found * 4 + k) * n, s4[k], n); memcpy(oE + ((size_t)found * 4 + k) * n, e4[k], n); }
                            found++;
                            if (stop_after > 0 && found >= stop_after) { found = stop_after; goto done; }
                        }
                    }
                    i = (i + 1) & cm;
                }
            }
        }
    }
done:
    if (stats) { stats[0] = st0; stats[1] = st1; stats[2] = st2; stats[3] = (double)Lw[0]; stats[4] = (double)Lw[1]; stats[5] = (double)Lw[2]; stats[6] = (double)Lw[3]; }
    free(bm); free(tab); free(h[0]); free(h[1]);
    for (int m = 0; m < 2; m++) { free(bs[m]); free(sh_[m]); free(io[m]); }
    for (int b = 0; b < 4; b++) if (map_[b]) free(map_[b]);
    return found;
}

void gl1_gen(int n, int r, int nu, const int *w, const int *forced, const int *poff, const int *plist, const double *dpw, int nfree, const int *freeidx, int fw,
             int nw, int mw, uint64_t seed, uint64_t j, const uint8_t *rok, int RM, int cansign, const uint8_t *sok, const uint8_t *eok, int SM, int EM, int fs, int fe, int8_t *s, int8_t *e)
{
    Ctx cx = { n, r, 0, NULL, nu, w, forced, poff, plist, dpw, nfree, freeidx, fw, {0} };
    Blk b = { nw, 0, seed, NULL, NULL, rok, RM, cansign, mw, sok, eok, SM, EM, fs, fe };
    gen_cfg(&cx, &b, j, s, e);
}
