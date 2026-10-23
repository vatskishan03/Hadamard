// qp.c : two-level exact engine for the QUATERNARY skew Paley plug-in (round 36, plug lens).
// G = K x F_q, blocks U: G -> {1,i,-1,-i}, H_K-invariant, U = s(x)delta + t(x)Q + e(x)J, values x0 (f=0), x1 (squares), x2 (non-squares) per K-orbit.
// level 1: w = x0 + (q-1)e  (Gaussian integers):  sum_X PAF(w_X) = 4|K|q delta_K            (class (W,R) restricted hash join, A x B | C x D)
// level 2: t given s,e:     sum_X (s s* + q t t*) = 4|K|q delta, sum_X (t s* - s t*) = 0   (scaled by 2: exact Gaussian integers; 2x2 meet in the middle)
// block A skew: A(-g) = i conj A(g), A(0) in {1,-i}; blocks B,C,D normalised x0(0)=1.  Every filter is a necessary condition.
// usage: qp shape.txt key=value ...   keys: maxprint first l2 psd wcap sample seed part nparts
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
typedef long long ll; typedef unsigned long long u64; typedef unsigned int u32;
#define MAXNU 20
static int SB = 6, SBM = 63;
static int nu, q, Ksz, ell[MAXNU], neg_[MAXNU], Nt[MAXNU][MAXNU][MAXNU];
static double Cre[MAXNU][MAXNU], Cim[MAXNU][MAXNU];
static int nJ, Jl[MAXNU];
static u64 cf1[MAXNU][2], cfI[MAXNU][2], cfQ[MAXNU][2];
static u64 al[MAXNU][MAXNU], be[MAXNU][MAXNU], ga[MAXNU][MAXNU], de_[MAXNU][MAXNU];
static u64 T1, T2;
static const int UR[4] = {1, 0, -1, 0}, UI[4] = {0, 1, 0, -1};
static int opt_maxprint = 5, opt_first = 0, opt_l2 = 1, opt_psd = 1, opt_wcap = -1, opt_part = 0, opt_nparts = 1, opt_verbose = 0;
static double opt_sample[3] = {1, 1, 1}; static u64 opt_seed = 1;
static ll opt_l2cap = 400000000LL; static int opt_real = 0; static int opt_listonly = 0; static int opt_l2c = 0; static int opt_fastf = 1; static int opt_bE = 0; static int opt_d1 = 0; static double opt_l2p = 1.0; static ll n_l2hit = 0; static double t_l2 = 0;

typedef struct { int type, x0, e2r, e2i, wr, wi, c, s2r, s2i, nt, t2r[4], t2i[4]; } Sym;
static Sym sy[36]; static int nsy = 0;
static void mksyms(void) {
    static const int ER[4] = {1, -1, -1, 1}, EI[4] = {1, 1, -1, -1};
    nsy = 0;
    for (int a = 0; a < 4; a++) { Sym *s = &sy[nsy++]; memset(s, 0, sizeof(Sym)); s->type = 0; s->x0 = a; s->wr = UR[a]; s->wi = UI[a]; s->c = 0; s->nt = 4; for (int b = 0; b < 4; b++) { s->t2r[b] = 2 * UR[b]; s->t2i[b] = 2 * UI[b]; } }
    for (int a = 0; a < 4; a++) for (int b = 0; b < 4; b++) { Sym *s = &sy[nsy++]; memset(s, 0, sizeof(Sym)); s->type = 1; s->x0 = a; s->e2r = 2 * UR[b]; s->e2i = 2 * UI[b]; s->wr = UR[a] + (q - 1) * UR[b]; s->wi = UI[a] + (q - 1) * UI[b]; s->c = (q - 1) + 2 * UR[(a - b) & 3]; s->nt = 1; }
    for (int a = 0; a < 4; a++) for (int m = 0; m < 4; m++) { Sym *s = &sy[nsy++]; memset(s, 0, sizeof(Sym)); s->type = 2; s->x0 = a; s->e2r = ER[m]; s->e2i = EI[m]; s->wr = UR[a] + (q - 1) / 2 * ER[m]; s->wi = UI[a] + (q - 1) / 2 * EI[m]; s->c = (q - 1) / 2 + UR[a] * ER[m] + UI[a] * EI[m]; s->nt = 2; s->t2r[0] = -EI[m]; s->t2i[0] = ER[m]; s->t2r[1] = EI[m]; s->t2i[1] = -ER[m]; }
    for (int i = 0; i < nsy; i++) { sy[i].s2r = 2 * UR[sy[i].x0] - sy[i].e2r; sy[i].s2i = 2 * UI[sy[i].x0] - sy[i].e2i; }
    if (opt_real) {   // real alphabet: x0, x1, x2 in {+-1}: types E and F only, x0 even exponent, E: t = +-1 (2t = +-2), F: x1 = +-1
        int k = 0;
        for (int i = 0; i < nsy; i++) {
            Sym t = sy[i]; if (t.x0 & 1) continue; if (t.type == 2) continue;
            if (t.type == 1 && (t.e2i != 0)) continue;     // F: x1 = i^b must be real: e2 = 2 i^b real
            if (t.type == 0) { t.nt = 2; t.t2r[0] = 2; t.t2i[0] = 0; t.t2r[1] = -2; t.t2i[1] = 0; }
            sy[k++] = t;
        }
        nsy = k;
    }
}
static u64 rs_state;
static inline u64 xs(void) { u64 x = rs_state; x ^= x >> 12; x ^= x << 25; x ^= x >> 27; rs_state = x; return x * 0x2545F4914F6CDD1DULL; }
static inline u64 sm64(u64 *x) { u64 z = (*x += 0x9E3779B97F4A7C15ULL); z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL; z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL; return z ^ (z >> 31); }

// ---------------------------------------------------------------- level-1 hash / vector
static u64 hash1(const int *wr, const int *wi) {
    u64 h = 0;
    for (int jj = 0; jj < nJ; jj++) {
        int j = Jl[jj]; ll re = 0, im = 0;
        for (int a = 0; a < nu; a++) for (int b = 0; b < nu; b++) { int n = Nt[j][a][b]; if (!n) continue; re += (ll)n * ((ll)wr[a] * wr[b] + (ll)wi[a] * wi[b]); im += (ll)n * ((ll)wi[a] * wr[b] - (ll)wr[a] * wi[b]); }
        h += cf1[jj][0] * (u64)re + cf1[jj][1] * (u64)im;
    }
    return h;
}
static void vec1(const int *wr, const int *wi, ll *v) {  // v[2*jj], v[2*jj+1]
    for (int jj = 0; jj < nJ; jj++) {
        int j = Jl[jj]; ll re = 0, im = 0;
        for (int a = 0; a < nu; a++) for (int b = 0; b < nu; b++) { int n = Nt[j][a][b]; if (!n) continue; re += (ll)n * ((ll)wr[a] * wr[b] + (ll)wi[a] * wi[b]); im += (ll)n * ((ll)wi[a] * wr[b] - (ll)wr[a] * wi[b]); }
        v[2 * jj] = re; v[2 * jj + 1] = im;
    }
}
// expand a config to per-orbit level-1 data.  kind 0 = A (skew), 1 = free.
// cfg: 6 bits per orbit = symbol index (A: only orbit 0 and the pair representatives o < neg[o] are used)
static void expand1(int kind, u64 cfg, int *wr, int *wi, int *ty, int *cw) {
    for (int a = 0; a < nu; a++) {
        if (kind == 1 || a == 0 || a < neg_[a]) { int s = (cfg >> (SB * a)) & SBM; wr[a] = sy[s].wr; wi[a] = sy[s].wi; if (ty) ty[a] = sy[s].type; if (cw) cw[a] = sy[s].c; }
    }
    if (kind == 0) for (int a = 1; a < nu; a++) if (a < neg_[a]) { int p = neg_[a]; if (opt_real) { wr[p] = -wr[a]; wi[p] = -wi[a]; } else { wr[p] = wi[a]; wi[p] = wr[a]; } if (ty) ty[p] = ty[a]; if (cw) cw[p] = cw[a]; }
}
// ---------------------------------------------------------------- lists
typedef struct { u64 key, cfg; int W, Rr, Ri; int pad; } El;
typedef struct { El *v; long n, cap; } List;
static void lpush(List *L, El *e) { if (L->n == L->cap) { L->cap = L->cap ? L->cap * 2 : 1024; L->v = realloc(L->v, sizeof(El) * L->cap); if (!L->v) { fprintf(stderr, "OOM list\n"); exit(2); } } L->v[L->n++] = *e; }
static int psd_ok(const int *wr, const int *wi) {
    double lim = 4.0 * Ksz * q + 1e-6;
    for (int j = 0; j < nu; j++) {
        double re = 0, im = 0;
        for (int a = 0; a < nu; a++) { re += wr[a] * Cre[j][a] - wi[a] * Cim[j][a]; im += wr[a] * Cim[j][a] + wi[a] * Cre[j][a]; }
        if (re * re + im * im > lim) return 0;
    }
    return 1;
}
static void leaf(List *L, int kind, u64 cfg, int W, const int *wr, const int *wi, double psamp) {
    if (psamp < 1.0) { double u = (double)(xs() >> 11) * (1.0 / 9007199254740992.0); if (u >= psamp) return; }
    if (opt_psd && !psd_ok(wr, wi)) return;
    El e; e.cfg = cfg; e.W = W; e.pad = 0; ll Rr = 0, Ri = 0;
    for (int a = 0; a < nu; a++) { Rr += (ll)ell[a] * wr[a]; Ri += (ll)ell[a] * wi[a]; }
    e.Rr = (int)Rr; e.Ri = (int)Ri; e.key = hash1(wr, wi);
    lpush(L, &e);
}
static int gwcap;  // cap on defect weight of a block
static double rthin_free = 1.0, rthin_skew = 1.0;
static inline int thin_skip(double r) { return r < 1.0 && ((double)(xs() >> 11) * (1.0 / 9007199254740992.0) >= r); }
static void rec_free(List *L, int o, int W, u64 cfg, int *wr, int *wi, double ps) {
    if (o == nu) { leaf(L, 1, cfg, W, wr, wi, 1.0); return; }
    for (int s = 0; s < nsy; s++) {
        if (o == 0 && sy[s].x0 != 0) continue;
        if (opt_d1 && ell[o] > 1 && sy[s].type != 0) continue;
        if (o > 0 && thin_skip(rthin_free)) continue;
        int W2 = W + ell[o] * sy[s].c; if (W2 > gwcap) continue;
        wr[o] = sy[s].wr; wi[o] = sy[s].wi;
        rec_free(L, o + 1, W2, cfg | ((u64)s << (SB * o)), wr, wi, ps);
    }
}
static void rec_skew(List *L, int o, int W, u64 cfg, int *wr, int *wi, double ps) {
    while (o < nu && o > 0 && !(o < neg_[o])) o++;
    if (o == nu) { leaf(L, 0, cfg, W, wr, wi, 1.0); return; }
    for (int s = 0; s < nsy; s++) {
        if (o > 0 && opt_d1 && ell[o] > 1 && sy[s].type != 0) continue;
        if (o > 0 && thin_skip(rthin_skew)) continue;
        if (o == 0 && opt_real) { if (!(sy[s].type == 0 && sy[s].x0 == 0)) continue; }
        else if (o == 0) { if (sy[s].type != 2) continue; if (!((sy[s].e2r == 1 && sy[s].e2i == 1) || (sy[s].e2r == -1 && sy[s].e2i == -1))) continue; if (!(sy[s].x0 == 0 || sy[s].x0 == 3)) continue; }
        int W2 = W + (o == 0 ? 1 : 2) * ell[o] * sy[s].c; if (W2 > gwcap) continue;
        wr[o] = sy[s].wr; wi[o] = sy[s].wi;
        if (o > 0) { int p = neg_[o]; if (opt_real) { wr[p] = -wr[o]; wi[p] = -wi[o]; } else { wr[p] = wi[o]; wi[p] = wr[o]; } }
        rec_skew(L, o + 1, W2, cfg | ((u64)s << (SB * o)), wr, wi, ps);
    }
}
static int cmp_cls(const void *x, const void *y) { const El *a = x, *b = y; if (a->W != b->W) return a->W < b->W ? -1 : 1; if (a->Rr != b->Rr) return a->Rr < b->Rr ? -1 : 1; if (a->Ri != b->Ri) return a->Ri < b->Ri ? -1 : 1; return 0; }
typedef struct { int W, Rr, Ri; ll nrm; long lo, hi; } Cls;
static Cls *mkcls(List *L, long *nc) { Cls *c = malloc(sizeof(Cls) * (L->n + 1)); long k = 0; for (long i = 0; i < L->n;) { long j = i; while (j < L->n && cmp_cls(&L->v[i], &L->v[j]) == 0) j++; c[k].W = L->v[i].W; c[k].Rr = L->v[i].Rr; c[k].Ri = L->v[i].Ri; c[k].nrm = (ll)c[k].Rr * c[k].Rr + (ll)c[k].Ri * c[k].Ri; c[k].lo = i; c[k].hi = j; k++; i = j; } *nc = k; return c; }

// ---------------------------------------------------------------- level 2
typedef struct { int sr[MAXNU], si[MAXNU], er[MAXNU], ei[MAXNU]; int nslot, slotorb[MAXNU], nopt[MAXNU]; int optr[MAXNU][4], opti[MAXNU][4]; } Blk2;
static void blk_setup(int kind, u64 cfg, Blk2 *b) {
    b->nslot = 0;
    for (int a = 0; a < nu; a++) {
        if (kind == 1 || a == 0 || a < neg_[a]) {
            int s = (cfg >> (SB * a)) & SBM; b->sr[a] = sy[s].s2r; b->si[a] = sy[s].s2i; b->er[a] = sy[s].e2r; b->ei[a] = sy[s].e2i;
            int k = b->nslot++; b->slotorb[k] = a; b->nopt[k] = sy[s].nt; for (int c = 0; c < sy[s].nt; c++) { b->optr[k][c] = sy[s].t2r[c]; b->opti[k][c] = sy[s].t2i[c]; }
        }
    }
    if (kind == 0) for (int a = 1; a < nu; a++) if (a < neg_[a]) { int p = neg_[a]; if (opt_real) { b->sr[p] = -b->sr[a]; b->si[p] = -b->si[a]; b->er[p] = -b->er[a]; b->ei[p] = -b->ei[a]; } else { b->sr[p] = b->si[a]; b->si[p] = b->sr[a]; b->er[p] = b->ei[a]; b->ei[p] = b->er[a]; } }
}
static void blk_tarr(int kind, const Blk2 *b, const int *dig, int *tr, int *ti) {
    for (int k = 0; k < b->nslot; k++) { int a = b->slotorb[k]; tr[a] = b->optr[k][dig[k]]; ti[a] = b->opti[k][dig[k]]; if (kind == 0 && a > 0) { int p = neg_[a]; if (opt_real) { tr[p] = tr[a]; ti[p] = ti[a]; } else { tr[p] = -ti[a]; ti[p] = -tr[a]; } } }
}
static u64 hash2(const int *sr, const int *si, const int *tr, const int *ti) {
    u64 h = 0;
    for (int a = 0; a < nu; a++) for (int b = 0; b < nu; b++) {
        ll ssr = (ll)sr[a] * sr[b] + (ll)si[a] * si[b], ssi = (ll)si[a] * sr[b] - (ll)sr[a] * si[b];
        ll ttr = (ll)tr[a] * tr[b] + (ll)ti[a] * ti[b], tti = (ll)ti[a] * tr[b] - (ll)tr[a] * ti[b];
        ll tsr = (ll)tr[a] * sr[b] + (ll)ti[a] * si[b], tsi = (ll)ti[a] * sr[b] - (ll)tr[a] * si[b];
        ll str_ = (ll)sr[a] * tr[b] + (ll)si[a] * ti[b], sti = (ll)si[a] * tr[b] - (ll)sr[a] * ti[b];
        h += al[a][b] * (u64)(ssr + (ll)q * ttr) + be[a][b] * (u64)(ssi + (ll)q * tti) + ga[a][b] * (u64)(tsr - str_) + de_[a][b] * (u64)(tsi - sti);
    }
    return h;
}
static void vec2(const int *sr, const int *si, const int *tr, const int *ti, ll *v) {  // 4 per jj : I re, I im, Q re, Q im
    for (int jj = 0; jj < nJ; jj++) {
        int j = Jl[jj]; ll Ir = 0, Ii = 0, Qr = 0, Qi = 0;
        for (int a = 0; a < nu; a++) for (int b = 0; b < nu; b++) {
            int n = Nt[j][a][b]; if (!n) continue;
            ll ssr = (ll)sr[a] * sr[b] + (ll)si[a] * si[b], ssi = (ll)si[a] * sr[b] - (ll)sr[a] * si[b];
            ll ttr = (ll)tr[a] * tr[b] + (ll)ti[a] * ti[b], tti = (ll)ti[a] * tr[b] - (ll)tr[a] * ti[b];
            ll tsr = (ll)tr[a] * sr[b] + (ll)ti[a] * si[b], tsi = (ll)ti[a] * sr[b] - (ll)tr[a] * si[b];
            ll str_ = (ll)sr[a] * tr[b] + (ll)si[a] * ti[b], sti = (ll)si[a] * tr[b] - (ll)sr[a] * ti[b];
            Ir += n * (ssr + (ll)q * ttr); Ii += n * (ssi + (ll)q * tti); Qr += n * (tsr - str_); Qi += n * (tsi - sti);
        }
        v[4 * jj] = Ir; v[4 * jj + 1] = Ii; v[4 * jj + 2] = Qr; v[4 * jj + 3] = Qi;
    }
}
typedef struct { u64 h; u32 idx; } Fe;
static int cmp_n(const void *a, const void *b) { const ll *x = a, *y = b; return (x[0] > y[0]) - (x[0] < y[0]); }
static void digits_of(const Blk2 *b, u32 idx, int *dig) { for (int k = 0; k < b->nslot; k++) { dig[k] = idx & 3; idx >>= 2; } }
static Fe *fibre(int kind, const Blk2 *b, long *nout) {
    long tot = 1; for (int k = 0; k < b->nslot; k++) tot *= b->nopt[k];
    Fe *F = malloc(sizeof(Fe) * tot); long m = 0;
    int dig[MAXNU] = {0}, tr[MAXNU], ti[MAXNU];
    for (;;) {
        blk_tarr(kind, b, dig, tr, ti);
        u32 idx = 0; for (int k = b->nslot - 1; k >= 0; k--) idx = (idx << 2) | dig[k];
        F[m].h = hash2(b->sr, b->si, tr, ti); F[m].idx = idx; m++;
        int k = 0; while (k < b->nslot) { dig[k]++; if (dig[k] < b->nopt[k]) break; dig[k] = 0; k++; }
        if (k == b->nslot) break;
    }
    *nout = m; return F;
}
// global counters
static ll n_l1 = 0, n_l2sol = 0, n_l2inst = 0, n_l2big = 0, n_l2work = 0; static ll n_l1ok_exact = 0; static int printed = 0; static int stop_flag = 0;
static void print_solution(const int kinds[4], const u64 cfgs[4], const Blk2 B2[4], const u32 idx[4]) {
    printf("SOL");
    for (int X = 0; X < 4; X++) {
        int dig[MAXNU], tr[MAXNU], ti[MAXNU]; digits_of(&B2[X], idx[X], dig); blk_tarr(kinds[X], &B2[X], dig, tr, ti);
        char s0[MAXNU + 1], s1[MAXNU + 1], s2[MAXNU + 1];
        for (int a = 0; a < nu; a++) {
            // x0 = (s2+e2)/2, x1 = (e2+t2)/2, x2 = (e2-t2)/2
            int v0r = (B2[X].sr[a] + B2[X].er[a]) / 2, v0i = (B2[X].si[a] + B2[X].ei[a]) / 2;
            int v1r = (B2[X].er[a] + tr[a]) / 2, v1i = (B2[X].ei[a] + ti[a]) / 2;
            int v2r = (B2[X].er[a] - tr[a]) / 2, v2i = (B2[X].ei[a] - ti[a]) / 2;
            int ex[3];
            int vr[3] = {v0r, v1r, v2r}, vi[3] = {v0i, v1i, v2i};
            for (int c = 0; c < 3; c++) { ex[c] = -1; for (int u = 0; u < 4; u++) if (UR[u] == vr[c] && UI[u] == vi[c]) ex[c] = u; if (ex[c] < 0) { fprintf(stderr, "BAD unit at block %d orbit %d comp %d (%d,%d)\n", X, a, c, vr[c], vi[c]); ex[c] = 9; } }
            s0[a] = '0' + ex[0]; s1[a] = '0' + ex[1]; s2[a] = '0' + ex[2];
        }
        s0[nu] = s1[nu] = s2[nu] = 0;
        printf(" %s|%s|%s", s0, s1, s2);
    }
    printf("\n"); fflush(stdout);
}
static int exact2(const int kinds[4], const Blk2 B2[4], const u32 idx[4]) {
    ll tot[4 * MAXNU]; memset(tot, 0, sizeof(tot));
    for (int X = 0; X < 4; X++) {
        int dig[MAXNU], tr[MAXNU], ti[MAXNU]; digits_of(&B2[X], idx[X], dig); blk_tarr(kinds[X], &B2[X], dig, tr, ti);
        ll v[4 * MAXNU]; vec2(B2[X].sr, B2[X].si, tr, ti, v); for (int i = 0; i < 4 * nJ; i++) tot[i] += v[i];
    }
    for (int i = 0; i < 4 * nJ; i++) { ll want = (i == 0) ? 16LL * Ksz * q : 0; if (tot[i] != want) return 0; }
    return 1;
}
typedef struct { u64 key; u32 i, j; } Ent;
static void level2(const u64 cfg[4]) {
    ll sol0 = n_l2sol; clock_t tc0 = clock();
    int kinds[4] = {0, 1, 1, 1};
    Blk2 B2[4]; Fe *F[4]; long n[4];
    for (int X = 0; X < 4; X++) { blk_setup(kinds[X], cfg[X], &B2[X]); F[X] = fibre(kinds[X], &B2[X], &n[X]); }
    n_l2inst++;
    int ord[4] = {0, 1, 2, 3};
    for (int i = 0; i < 4; i++) for (int j = i + 1; j < 4; j++) if (n[ord[j]] < n[ord[i]]) { int t = ord[i]; ord[i] = ord[j]; ord[j] = t; }
    int X0 = ord[0], X1 = ord[1], X2 = ord[2], X3 = ord[3];
    ll pl = n[X0] * n[X3], pr = n[X1] * n[X2];
    int ta0, ta1, tb0, tb1; ll ps, pp;
    if (pl <= pr) { ta0 = X0; ta1 = X3; tb0 = X1; tb1 = X2; ps = pl; pp = pr; } else { ta0 = X1; ta1 = X2; tb0 = X0; tb1 = X3; ps = pr; pp = pl; }
    if (ps > 60000000LL || pp > opt_l2cap) { n_l2big++; for (int X = 0; X < 4; X++) free(F[X]); return; }
    n_l2work += ps + pp;
    size_t tsz = 1; while (tsz < (size_t)ps * 2) tsz <<= 1;
    Ent *tab = calloc(tsz, sizeof(Ent)); u64 *bm; size_t bmbits = tsz * 4; bm = calloc(bmbits / 64 + 1, 8);
    for (long i = 0; i < n[ta0]; i++) for (long j = 0; j < n[ta1]; j++) {
        u64 k = F[ta0][i].h + F[ta1][j].h; size_t h = (k * 0x9E3779B97F4A7C15ULL) >> 20; h &= tsz - 1;
        while (tab[h].i | tab[h].j) h = (h + 1) & (tsz - 1);
        tab[h].key = k; tab[h].i = (u32)i + 1; tab[h].j = (u32)j + 1;
        u64 bi = (k * 0xC2B2AE3D27D4EB4FULL) >> 20; bi &= bmbits - 1; bm[bi >> 6] |= 1ULL << (bi & 63);
    }
    for (long i = 0; i < n[tb0]; i++) {
        u64 hi = T2 - F[tb0][i].h;
        for (long j = 0; j < n[tb1]; j++) {
            u64 k = hi - F[tb1][j].h;
            { u64 bi = (k * 0xC2B2AE3D27D4EB4FULL) >> 20; bi &= bmbits - 1; if (!((bm[bi >> 6] >> (bi & 63)) & 1)) continue; }
            size_t h = (k * 0x9E3779B97F4A7C15ULL) >> 20; h &= tsz - 1;
            while (tab[h].i | tab[h].j) {
                if (tab[h].key == k) {
                    u32 idx[4]; idx[ta0] = F[ta0][tab[h].i - 1].idx; idx[ta1] = F[ta1][tab[h].j - 1].idx; idx[tb0] = F[tb0][i].idx; idx[tb1] = F[tb1][j].idx;
                    if (exact2(kinds, B2, idx)) {
                        n_l2sol++;
                        if (printed < opt_maxprint) { printed++; print_solution(kinds, cfg, B2, idx); }
                        if (opt_first && printed >= opt_maxprint) stop_flag = 1;
                    }
                }
                h = (h + 1) & (tsz - 1);
            }
        }
    }
    free(tab); free(bm); for (int X = 0; X < 4; X++) free(F[X]);
    if (n_l2sol > sol0) n_l2hit++;
    t_l2 += (double)(clock() - tc0) / CLOCKS_PER_SEC;
}

// ---------------------------------------------------------------- level 2 with DC classes (T = sum ell t; n = |T|^2, l = Im(T conj S)): class-bucketed 2x2 meet in the middle
typedef struct { u64 h; u32 idx; int n, l; } Fe2;
static int cmp_fe2(const void *a, const void *b) { const Fe2 *x = a, *y = b; if (x->n != y->n) return x->n < y->n ? -1 : 1; if (x->l != y->l) return x->l < y->l ? -1 : 1; return 0; }
static Fe2 *fibre2(int kind, const Blk2 *b, long *nout) {
    long tot = 1; for (int k = 0; k < b->nslot; k++) tot *= b->nopt[k];
    Fe2 *F = malloc(sizeof(Fe2) * tot); long m = 0;
    int dig[MAXNU] = {0}, tr[MAXNU], ti[MAXNU];
    ll S2r = 0, S2i = 0; for (int a = 0; a < nu; a++) { S2r += (ll)ell[a] * b->sr[a]; S2i += (ll)ell[a] * b->si[a]; }
    for (;;) {
        blk_tarr(kind, b, dig, tr, ti);
        u32 idx = 0; for (int k = b->nslot - 1; k >= 0; k--) idx = (idx << 2) | dig[k];
        ll Tr = 0, Ti = 0; for (int a = 0; a < nu; a++) { Tr += (ll)ell[a] * tr[a]; Ti += (ll)ell[a] * ti[a]; }
        F[m].h = hash2(b->sr, b->si, tr, ti); F[m].idx = idx; F[m].n = (int)(Tr * Tr + Ti * Ti); F[m].l = (int)(Ti * S2r - Tr * S2i); m++;
        int k = 0; while (k < b->nslot) { dig[k]++; if (dig[k] < b->nopt[k]) break; dig[k] = 0; k++; }
        if (k == b->nslot) break;
    }
    qsort(F, m, sizeof(Fe2), cmp_fe2);
    *nout = m; return F;
}

// fast fibre for FREE blocks (kind 1): the hash is a quadratic form in the t's -> depth-first accumulation
static u64 fD[MAXNU][4], fP[MAXNU][MAXNU][4][4];
static int f_nu_slots; static const Blk2 *f_b; static Fe2 *f_out; static long f_m; static ll f_S2r, f_S2i;
static int f_dig[MAXNU], f_tr[MAXNU], f_ti[MAXNU];
static void frec(int a, u64 h, ll Tr, ll Ti) {
    if (a == nu) {
        u32 idx = 0; for (int k = nu - 1; k >= 0; k--) idx = (idx << 2) | f_dig[k];
        Fe2 *e = &f_out[f_m++]; e->h = h; e->idx = idx; e->n = (int)(Tr * Tr + Ti * Ti); e->l = (int)(Ti * f_S2r - Tr * f_S2i); return;
    }
    for (int c = 0; c < f_b->nopt[a]; c++) {
        u64 hh = h + fD[a][c];
        for (int b = 0; b < a; b++) hh += fP[b][a][f_dig[b]][c];
        f_dig[a] = c;
        frec(a + 1, hh, Tr + (ll)ell[a] * f_b->optr[a][c], Ti + (ll)ell[a] * f_b->opti[a][c]);
    }
}
static Fe2 *fibre2_fast(const Blk2 *b, long *nout) {
    long tot = 1; for (int k = 0; k < nu; k++) tot *= b->nopt[k];
    Fe2 *F = malloc(sizeof(Fe2) * tot);
    const int *sr = b->sr, *si = b->si;
    u64 C0 = 0;
    for (int a = 0; a < nu; a++) for (int c = 0; c < nu; c++) { ll ssr = (ll)sr[a] * sr[c] + (ll)si[a] * si[c], ssi = (ll)si[a] * sr[c] - (ll)sr[a] * si[c]; C0 += al[a][c] * (u64)ssr + be[a][c] * (u64)ssi; }
    u64 Lr[MAXNU], Li[MAXNU];
    for (int a = 0; a < nu; a++) {
        u64 lr = 0, li = 0;
        for (int c = 0; c < nu; c++) {
            lr += ga[a][c] * (u64)(ll)sr[c] - de_[a][c] * (u64)(ll)si[c];
            li += ga[a][c] * (u64)(ll)si[c] + de_[a][c] * (u64)(ll)sr[c];
            lr += -ga[c][a] * (u64)(ll)sr[c] - de_[c][a] * (u64)(ll)si[c];
            li += -ga[c][a] * (u64)(ll)si[c] + de_[c][a] * (u64)(ll)sr[c];
        }
        Lr[a] = lr; Li[a] = li;
    }
    for (int a = 0; a < nu; a++) for (int c = 0; c < b->nopt[a]; c++) {
        ll tr = b->optr[a][c], ti = b->opti[a][c];
        fD[a][c] = (u64)tr * Lr[a] + (u64)ti * Li[a] + (u64)q * al[a][a] * (u64)(tr * tr + ti * ti);
    }
    for (int a = 0; a < nu; a++) for (int d = a + 1; d < nu; d++) for (int ca = 0; ca < b->nopt[a]; ca++) for (int cd = 0; cd < b->nopt[d]; cd++) {
        ll tra = b->optr[a][ca], tia = b->opti[a][ca], trd = b->optr[d][cd], tid = b->opti[d][cd];
        ll ttr = tra * trd + tia * tid, tti = tia * trd - tra * tid;
        fP[a][d][ca][cd] = (u64)q * (((al[a][d] + al[d][a]) * (u64)ttr) + ((be[a][d] - be[d][a]) * (u64)tti));
    }
    f_b = b; f_out = F; f_m = 0; f_S2r = 0; f_S2i = 0; for (int a = 0; a < nu; a++) { f_S2r += (ll)ell[a] * sr[a]; f_S2i += (ll)ell[a] * si[a]; }
    frec(0, C0, 0, 0);
    qsort(F, f_m, sizeof(Fe2), cmp_fe2);
    *nout = f_m; return F;
}
typedef struct { int n, l; long lo, hi; } C2;
static C2 *mkc2(Fe2 *F, long m, long *nc) { C2 *c = malloc(sizeof(C2) * (m + 1)); long k = 0; for (long i = 0; i < m;) { long j = i; while (j < m && F[j].n == F[i].n && F[j].l == F[i].l) j++; c[k].n = F[i].n; c[k].l = F[i].l; c[k].lo = i; c[k].hi = j; k++; i = j; } *nc = k; return c; }
static ll n_dcfail = 0, n_l2cprobe = 0;
static void level2c(const u64 cfg[4]) {
    ll sol0 = n_l2sol; clock_t tc0 = clock();
    int kinds[4] = {0, 1, 1, 1};
    Blk2 B2[4];
    ll sumS = 0;
    for (int X = 0; X < 4; X++) { blk_setup(kinds[X], cfg[X], &B2[X]); ll Sr = 0, Si = 0; for (int a = 0; a < nu; a++) { Sr += (ll)ell[a] * B2[X].sr[a]; Si += (ll)ell[a] * B2[X].si[a]; } sumS += Sr * Sr + Si * Si; }
    ll rest = 16LL * Ksz * q - sumS;
    if (rest < 0 || rest % q) { n_dcfail++; return; }
    ll Cn = rest / q;
    if (opt_verbose >= 2) { fprintf(stderr, "L2DBG Cn=%lld", Cn); for (int X = 0; X < 4; X++) { ll Sr = 0, Si = 0; int nE = 0; for (int a = 0; a < nu; a++) { Sr += (ll)ell[a] * B2[X].sr[a]; Si += (ll)ell[a] * B2[X].si[a]; if (B2[X].nopt[a] == 4) nE++; } fprintf(stderr, " | S2=(%lld,%lld) nE=%d", Sr, Si, nE); } fprintf(stderr, "\n"); }
    Fe2 *F[4]; long n[4]; C2 *cc[4]; long nc[4];
    for (int X = 0; X < 4; X++) { F[X] = (opt_fastf && kinds[X] == 1) ? fibre2_fast(&B2[X], &n[X]) : fibre2(kinds[X], &B2[X], &n[X]); cc[X] = mkc2(F[X], n[X], &nc[X]); }
    n_l2inst++;
    int ord[4] = {0, 1, 2, 3};
    for (int i = 0; i < 4; i++) for (int j = i + 1; j < 4; j++) if (n[ord[j]] < n[ord[i]]) { int t = ord[i]; ord[i] = ord[j]; ord[j] = t; }
    int a0 = ord[0], a1 = ord[1], b0 = ord[2], b1 = ord[3];
    ll ps = n[a0] * n[a1];
    if (ps > 60000000LL || n[b0] * n[b1] > 4000000000LL) { n_l2big++; for (int X = 0; X < 4; X++) { free(F[X]); free(cc[X]); } return; }
    n_l2work += ps;
    size_t tsz = 1; while (tsz < (size_t)ps * 2) tsz <<= 1;
    Ent *tab = calloc(tsz, sizeof(Ent)); size_t bmbits = tsz * 4; u64 *bm = calloc(bmbits / 64 + 1, 8);
    for (long i = 0; i < n[a0]; i++) for (long j = 0; j < n[a1]; j++) {
        u64 k = F[a0][i].h + F[a1][j].h; size_t h = (k * 0x9E3779B97F4A7C15ULL) >> 20; h &= tsz - 1;
        while (tab[h].i | tab[h].j) h = (h + 1) & (tsz - 1);
        tab[h].key = k; tab[h].i = (u32)i + 1; tab[h].j = (u32)j + 1;
        u64 bi = (k * 0xC2B2AE3D27D4EB4FULL) >> 20; bi &= bmbits - 1; bm[bi >> 6] |= 1ULL << (bi & 63);
    }
    // left class-pair set
    size_t ssz = 1; { size_t cnt = (size_t)nc[a0] * nc[a1]; while (ssz < 2 * cnt + 16) ssz <<= 1; }
    u64 *sset = malloc(8 * ssz); memset(sset, 0xff, 8 * ssz);
    for (long i = 0; i < nc[a0]; i++) for (long j = 0; j < nc[a1]; j++) {
        ll nn = (ll)cc[a0][i].n + cc[a1][j].n, ll_ = (ll)cc[a0][i].l + cc[a1][j].l;
        u64 key = ((u64)nn << 32) | (u32)(ll_ + (1LL << 30)); size_t h = (key * 0x9E3779B97F4A7C15ULL) >> 24; h &= ssz - 1; int f = 0;
        while (sset[h] != ~0ULL) { if (sset[h] == key) { f = 1; break; } h = (h + 1) & (ssz - 1); }
        if (!f) sset[h] = key;
    }
    for (long ci = 0; ci < nc[b0]; ci++) for (long cj = 0; cj < nc[b1]; cj++) {
        ll nn = Cn - cc[b0][ci].n - cc[b1][cj].n; ll l2_ = -((ll)cc[b0][ci].l + cc[b1][cj].l);
        if (nn < 0) continue;
        u64 key = ((u64)nn << 32) | (u32)(l2_ + (1LL << 30)); size_t h = (key * 0x9E3779B97F4A7C15ULL) >> 24; h &= ssz - 1; int f = 0;
        while (sset[h] != ~0ULL) { if (sset[h] == key) { f = 1; break; } h = (h + 1) & (ssz - 1); }
        if (!f) continue;
        for (long i = cc[b0][ci].lo; i < cc[b0][ci].hi; i++) {
            u64 hi = T2 - F[b0][i].h;
            for (long j = cc[b1][cj].lo; j < cc[b1][cj].hi; j++) {
                u64 k = hi - F[b1][j].h; n_l2cprobe++;
                { u64 bi = (k * 0xC2B2AE3D27D4EB4FULL) >> 20; bi &= bmbits - 1; if (!((bm[bi >> 6] >> (bi & 63)) & 1)) continue; }
                size_t hh = (k * 0x9E3779B97F4A7C15ULL) >> 20; hh &= tsz - 1;
                while (tab[hh].i | tab[hh].j) {
                    if (tab[hh].key == k) {
                        u32 idx[4]; idx[a0] = F[a0][tab[hh].i - 1].idx; idx[a1] = F[a1][tab[hh].j - 1].idx; idx[b0] = F[b0][i].idx; idx[b1] = F[b1][j].idx;
                        if (exact2(kinds, B2, idx)) {
                            n_l2sol++;
                            if (printed < opt_maxprint) { printed++; print_solution(kinds, cfg, B2, idx); }
                            if (opt_first && printed >= opt_maxprint) stop_flag = 1;
                        }
                    }
                    hh = (hh + 1) & (tsz - 1);
                }
            }
        }
    }
    free(tab); free(bm); free(sset); for (int X = 0; X < 4; X++) { free(F[X]); free(cc[X]); }
    if (n_l2sol > sol0) n_l2hit++;
    t_l2 += (double)(clock() - tc0) / CLOCKS_PER_SEC;
}
// ---------------------------------------------------------------- level-1 join
static List LA, LL;  // A list, shared B/C/D list
typedef struct { u64 key; u32 a, b; } Ent1;
static int exact1(const u64 cfg[4]) {
    ll tot[2 * MAXNU]; memset(tot, 0, sizeof(tot));
    for (int X = 0; X < 4; X++) { int wr[MAXNU], wi[MAXNU]; expand1(X == 0 ? 0 : 1, cfg[X], wr, wi, NULL, NULL); ll v[2 * MAXNU]; vec1(wr, wi, v); for (int i = 0; i < 2 * nJ; i++) tot[i] += v[i]; }
    for (int i = 0; i < 2 * nJ; i++) { ll want = (i == 0) ? 4LL * Ksz * q : 0; if (tot[i] != want) return 0; }
    return 1;
}

static ll n_pairs_left = 0, n_probe = 0, n_cand = 0;
static int opt_sym = 0; static long opt_resume = 0;
typedef struct { u64 key; } SKey;
static List LB_;
static void join1(void) {
    List *LBp = opt_bE ? &LB_ : &LL;
    long nA, nL, nB; Cls *cA = mkcls(&LA, &nA), *cL = mkcls(&LL, &nL), *cB = opt_bE ? mkcls(&LB_, &nB) : cL; if (!opt_bE) nB = nL;
    int Wtot = 4 * Ksz; ll Ntot = 4LL * Ksz * q;
    fprintf(stderr, "lists: |A|=%ld (%ld classes) |L|=%ld (%ld classes) |LB|=%ld (%ld classes)\n", LA.n, nA, LL.n, nL, LBp->n, nB);
    // (1) set of left sums required by right class pairs
    size_t rsz = 1; { size_t cnt = (size_t)nL * (size_t)nL; if (opt_sym) cnt = cnt / 2 + nL; while (rsz < 2 * cnt + 16) rsz <<= 1; if (rsz > ((size_t)1 << 28)) rsz = (size_t)1 << 28; }
    u64 *rset = malloc(8 * rsz); memset(rset, 0xff, 8 * rsz); int rfull = 0; size_t rcount = 0;
    for (long ic = 0; ic < nL; ic++) for (long id = (opt_sym ? ic : 0); id < nL; id++) {
        int Wr = Wtot - cL[ic].W - cL[id].W; ll Nr = Ntot - cL[ic].nrm - cL[id].nrm; if (Wr < 0 || Nr < 0) continue;
        u64 key = ((u64)Wr << 40) | (u64)Nr; size_t h = (key * 0x9E3779B97F4A7C15ULL) >> 24; h &= rsz - 1; int f = 0;
        while (rset[h] != ~0ULL) { if (rset[h] == key) { f = 1; break; } h = (h + 1) & (rsz - 1); }
        if (!f) { if (rcount * 2 > rsz) rfull = 1; else { rset[h] = key; rcount++; } }
    }
    fprintf(stderr, "right class pairs: %zu distinct required left sums%s\n", rcount, rfull ? " (FULL)" : "");
    // (2) compatible left class pairs
    char *compat = calloc((size_t)nA * nB, 1); ll *insA = calloc(nA, sizeof(ll));
    for (long ia = 0; ia < nA; ia++) for (long ib = 0; ib < nB; ib++) {
        int Wl = cA[ia].W + cB[ib].W; ll Nl = cA[ia].nrm + cB[ib].nrm; if (Wl > Wtot || Nl > Ntot) continue;
        int f = rfull; if (!f) { u64 key = ((u64)Wl << 40) | (u64)Nl; size_t h = (key * 0x9E3779B97F4A7C15ULL) >> 24; h &= rsz - 1; while (rset[h] != ~0ULL) { if (rset[h] == key) { f = 1; break; } h = (h + 1) & (rsz - 1); } }
        if (f) { compat[(size_t)ia * nB + ib] = 1; insA[ia] += cB[ib].hi - cB[ib].lo; }
    }
    long *clsA = malloc(sizeof(long) * (LA.n + 1)); for (long ia = 0; ia < nA; ia++) for (long x = cA[ia].lo; x < cA[ia].hi; x++) clsA[x] = ia;
    free(rset);
    ll cap = 16000000LL; long pidx = 0;
    for (long alo = 0; alo < LA.n && !stop_flag;) {
        long ahi = alo; ll acc = 0;
        while (ahi < LA.n && (acc + insA[clsA[ahi]] <= cap || ahi == alo)) { acc += insA[clsA[ahi]]; ahi++; }
        if (acc == 0) { alo = ahi; continue; }
        size_t nl = (size_t)acc; size_t tsz = 1; while (tsz < 2 * nl + 2) tsz <<= 1;
        Ent1 *tab = calloc(tsz, sizeof(Ent1)); size_t bmbits = tsz * 8; u64 *bm = calloc(bmbits / 64 + 1, 8);
        if (!tab || !bm) { fprintf(stderr, "OOM left table\n"); exit(2); }
        size_t ssz = 1; { size_t cnt = 0; for (long x = alo; x < ahi; x++) { cnt++; } cnt = (size_t)nA * (size_t)nL; while (ssz < 2 * cnt + 16) ssz <<= 1; if (ssz > ((size_t)1 << 26)) ssz = (size_t)1 << 26; }
        u64 *sset = malloc(8 * ssz); memset(sset, 0xff, 8 * ssz); size_t scount = 0; int sfull = 0;
        long lastcls = -1;
        for (long x = alo; x < ahi; x++) {
            long ia = clsA[x];
            if (ia != lastcls) { lastcls = ia;
                for (long ib = 0; ib < nB; ib++) if (compat[(size_t)ia * nB + ib]) {
                    int Wl = cA[ia].W + cB[ib].W; ll Nl = cA[ia].nrm + cB[ib].nrm; u64 key = ((u64)Wl << 40) | (u64)Nl; size_t h = (key * 0x9E3779B97F4A7C15ULL) >> 24; h &= ssz - 1; int f = 0;
                    while (sset[h] != ~0ULL) { if (sset[h] == key) { f = 1; break; } h = (h + 1) & (ssz - 1); }
                    if (!f) { if (scount * 2 > ssz) sfull = 1; else { sset[h] = key; scount++; } }
                } }
            for (long ib = 0; ib < nB; ib++) if (compat[(size_t)ia * nB + ib]) for (long y = cB[ib].lo; y < cB[ib].hi; y++) {
                u64 k = LA.v[x].key + LBp->v[y].key; size_t h = (k * 0x9E3779B97F4A7C15ULL) >> 20; h &= tsz - 1;
                while (tab[h].a | tab[h].b) h = (h + 1) & (tsz - 1);
                tab[h].key = k; tab[h].a = (u32)x + 1; tab[h].b = (u32)y + 1;
                u64 bi = (k * 0xC2B2AE3D27D4EB4FULL) >> 20; bi &= bmbits - 1; bm[bi >> 6] |= 1ULL << (bi & 63);
            }
        }
        n_pairs_left += nl;
        fprintf(stderr, "chunk A[%ld,%ld): left table %zu entries (compat-filtered), %zu distinct left class sums\n", alo, ahi, nl, scount);
        for (long ic = (alo == 0 ? opt_resume : 0); ic < nL && !stop_flag; ic++) {
            long id0 = opt_sym ? ic : 0;
            for (long id = id0; id < nL && !stop_flag; id++) {
                int Wr = Wtot - cL[ic].W - cL[id].W; ll Nr = Ntot - cL[ic].nrm - cL[id].nrm; if (Wr < 0 || Nr < 0) continue;
                u64 key = ((u64)Wr << 40) | (u64)Nr; size_t h = (key * 0x9E3779B97F4A7C15ULL) >> 24; h &= ssz - 1; int found = sfull;
                while (!found && sset[h] != ~0ULL) { if (sset[h] == key) { found = 1; break; } h = (h + 1) & (ssz - 1); }
                if (!found) continue;
                pidx++;
                if (opt_nparts > 1 && (pidx % opt_nparts) != opt_part) continue;
                for (long x = cL[ic].lo; x < cL[ic].hi && !stop_flag; x++) {
                    u64 hx = T1 - LL.v[x].key;
                    long y0 = (opt_sym && id == ic) ? x : cL[id].lo;
                    for (long y = y0; y < cL[id].hi; y++) {
                        u64 k = hx - LL.v[y].key; n_probe++;
                        { u64 bi = (k * 0xC2B2AE3D27D4EB4FULL) >> 20; bi &= bmbits - 1; if (!((bm[bi >> 6] >> (bi & 63)) & 1)) continue; }
                        size_t hh = (k * 0x9E3779B97F4A7C15ULL) >> 20; hh &= tsz - 1;
                        while (tab[hh].a | tab[hh].b) {
                            if (tab[hh].key == k) {
                                n_cand++;
                                u64 cfg[4] = {LA.v[tab[hh].a - 1].cfg, LBp->v[tab[hh].b - 1].cfg, LL.v[x].cfg, LL.v[y].cfg};
                                if (exact1(cfg)) {
                                    n_l1++;
                                    if (opt_l2) { if (opt_l2p >= 1.0 || (double)(xs() >> 11) * (1.0 / 9007199254740992.0) < opt_l2p) { if (opt_l2c) level2c(cfg); else level2(cfg); } }
                                    if (stop_flag) break;
                                }
                            }
                            hh = (hh + 1) & (tsz - 1);
                        }
                        if (stop_flag) break;
                    }
                }
            }
            if ((ic & 255) == 0 || ic == nL - 1) fprintf(stderr, "PROGRESS chunk %ld ic %ld/%ld l1=%lld l2inst=%lld l2sol=%lld probes=%lld t=%.0fs\n", alo, ic, nL, n_l1, n_l2inst, n_l2sol, n_probe, (double)clock() / CLOCKS_PER_SEC);
        }
        free(tab); free(bm); free(sset);
        alo = ahi;
    }
}
int main(int argc, char **argv) {
    FILE *f = fopen(argv[1], "r"); if (!f) { fprintf(stderr, "no shape\n"); return 1; }
    if (fscanf(f, "%d %d %d", &nu, &q, &Ksz) != 3) return 1;
    for (int i = 0; i < nu; i++) if (fscanf(f, "%d", &ell[i]) != 1) return 1;
    for (int i = 0; i < nu; i++) if (fscanf(f, "%d", &neg_[i]) != 1) return 1;
    for (int j = 0; j < nu; j++) for (int a = 0; a < nu; a++) for (int b = 0; b < nu; b++) if (fscanf(f, "%d", &Nt[j][a][b]) != 1) return 1;
    for (int j = 0; j < nu; j++) for (int a = 0; a < nu; a++) if (fscanf(f, "%lf %lf", &Cre[j][a], &Cim[j][a]) != 2) return 1;
    fclose(f);
    for (int i = 2; i < argc; i++) {
        char *eq = strchr(argv[i], '='); if (!eq) continue; *eq = 0; const char *k = argv[i], *v = eq + 1;
        if (!strcmp(k, "maxprint")) opt_maxprint = atoi(v); else if (!strcmp(k, "first")) opt_first = atoi(v); else if (!strcmp(k, "l2")) opt_l2 = atoi(v); else if (!strcmp(k, "psd")) opt_psd = atoi(v);
        else if (!strcmp(k, "wcap")) opt_wcap = atoi(v); else if (!strcmp(k, "part")) opt_part = atoi(v); else if (!strcmp(k, "nparts")) opt_nparts = atoi(v); else if (!strcmp(k, "seed")) opt_seed = strtoull(v, 0, 10);
        else if (!strcmp(k, "pA")) opt_sample[0] = atof(v); else if (!strcmp(k, "pL")) opt_sample[1] = atof(v); else if (!strcmp(k, "verbose")) opt_verbose = atoi(v); else if (!strcmp(k, "sym")) opt_sym = atoi(v); else if (!strcmp(k, "listonly")) opt_listonly = atoi(v); else if (!strcmp(k, "l2p")) opt_l2p = atof(v); else if (!strcmp(k, "l2c")) opt_l2c = atoi(v); else if (!strcmp(k, "fastf")) opt_fastf = atoi(v); else if (!strcmp(k, "bE")) opt_bE = atoi(v); else if (!strcmp(k, "real")) opt_real = atoi(v); else if (!strcmp(k, "d1")) opt_d1 = atoi(v); else if (!strcmp(k, "resume")) opt_resume = atol(v); else if (!strcmp(k, "l2cap")) opt_l2cap = atoll(v);
    }
    rs_state = opt_seed * 0x9E3779B97F4A7C15ULL | 1;
    if (opt_real) { SB = 3; SBM = 7; }
    mksyms();
    if (nu * SB > 64 || (nu > 16)) { fprintf(stderr, "nu too large for packing (nu*SB > 64 or nu > 16)\n"); return 3; }
    nJ = 0; Jl[nJ++] = 0; for (int j = 1; j < nu; j++) if (j < neg_[j]) Jl[nJ++] = j;
    { u64 st = 88172645463325252ULL; for (int jj = 0; jj < nJ; jj++) for (int c = 0; c < 2; c++) { cf1[jj][c] = sm64(&st) | 1; cfI[jj][c] = sm64(&st) | 1; cfQ[jj][c] = sm64(&st) | 1; } }
    for (int a = 0; a < nu; a++) for (int b = 0; b < nu; b++) { u64 x = 0, y = 0, z = 0, w = 0; for (int jj = 0; jj < nJ; jj++) { int n = Nt[Jl[jj]][a][b]; x += (u64)n * cfI[jj][0]; y += (u64)n * cfI[jj][1]; z += (u64)n * cfQ[jj][0]; w += (u64)n * cfQ[jj][1]; } al[a][b] = x; be[a][b] = y; ga[a][b] = z; de_[a][b] = w; }
    T1 = cf1[0][0] * (u64)(4LL * Ksz * q); T2 = cfI[0][0] * (u64)(16LL * Ksz * q);
    gwcap = opt_wcap >= 0 ? opt_wcap : 4 * Ksz;
    clock_t t0 = clock();
    int wr[MAXNU], wi[MAXNU];
    { int npair = 0; for (int a = 1; a < nu; a++) if (a < neg_[a]) npair++; rthin_skew = opt_sample[0] < 1.0 ? pow(opt_sample[0], 1.0 / (npair > 0 ? npair : 1)) : 1.0; rthin_free = opt_sample[1] < 1.0 ? pow(opt_sample[1], 1.0 / (nu - 1)) : 1.0; }
    rec_skew(&LA, 0, 0, 0, wr, wi, opt_sample[0]);
    rec_free(&LL, 0, 0, 0, wr, wi, opt_sample[1]);
    if (opt_bE) { int gw0 = gwcap; double rt0 = rthin_free; gwcap = 0; rthin_free = 1.0; rec_free(&LB_, 0, 0, 0, wr, wi, 1.0); gwcap = gw0; rthin_free = rt0; qsort(LB_.v, LB_.n, sizeof(El), cmp_cls); }
    qsort(LA.v, LA.n, sizeof(El), cmp_cls); qsort(LL.v, LL.n, sizeof(El), cmp_cls);
    fprintf(stderr, "lists built: |A|=%ld |L|=%ld  %.1fs\n", LA.n, LL.n, (double)(clock() - t0) / CLOCKS_PER_SEC);
    if (opt_listonly) { long nA_, nL_; Cls *c1 = mkcls(&LA, &nA_), *c2 = mkcls(&LL, &nL_); printf("LISTS |A|=%ld classesA=%ld |L|=%ld classesL=%ld\n", LA.n, nA_, LL.n, nL_); return 0; }
    join1();
    printf("STATS l2_hit_instances=%lld t_l2=%.1f l1_instances=%lld l1_cand=%lld left_pairs=%lld probes=%lld l2_instances=%lld l2_solutions=%lld l2_toobig=%lld l2_work=%lld dcfail=%lld l2cprobe=%lld time=%.1f\n", n_l2hit, t_l2, n_l1, n_cand, n_pairs_left, n_probe, n_l2inst, n_l2sol, n_l2big, n_l2work, n_dcfail, n_l2cprobe, (double)(clock() - t0) / CLOCKS_PER_SEC);
    printf("COUNT %lld (normalised: B,C,D x0(0)=1; ordered = 64x)\n", n_l2sol);
    return 0;
}
