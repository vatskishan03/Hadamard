// pl2f_fast.h : AVX2 fibre builder for blocks whose variable orbits are all singletons and whose lag lists are 1..nI / 1..nQ (cyclic K, H = 1: the C, D blocks of the plug-in).
// Same Gray-code enumeration, same filters (allowT, PSD <= T + 0.01), same float quantisation and same hash as build_fibre() (so the fibres are identical), but the per-flip updates are vectorised:
//   P (PAF of t at lags 1..nI) and Qd (Q at lags 1..nQ) as int16 vectors using doubled/reflected copies of t (no per-lag wraparound),  ar/ai (PSD channel sums) as double vectors.
static double CGRZ[MAXN][PS] __attribute__((aligned(32))), CGIZ[MAXN][PS] __attribute__((aligned(32)));   // [z][c] = cgr[c][z], cgi[c][z] (zero padded)
static void fast_tables(void) {
    for (int z = 0; z < n; z++) for (int c = 0; c < PS; c++) { CGRZ[z][c] = c < nch ? cgr[c][z] : 0.0; CGIZ[z][c] = c < nch ? cgi[c][z] : 0.0; }
}
static int fast_min_vars = 12;
static int fast_ok(int b) {
    int kd = b == 0 ? 0 : 1;
    if (nI > 16 || nQ > 16 || nI < 1 || nQ < 1) return 0;
    for (int a = 0; a < nI; a++) if (kI[a] != a + 1) return 0;
    for (int j = 0; j < nQ; j++) if (kQ[j] != j + 1) return 0;
    for (int j = 0; j < nvar[b]; j++) if (omm[kd][vorb[b][j]] != 1) return 0;
    if (nvar[b] < fast_min_vars) return 0;
    return 1;
}
static void build_fibre_fast(int b, Fib *F, int fixfirst) {
    int kd = b == 0 ? 0 : 1; int nvb = nvar[b];
    int nbits = (fixfirst && nvb > 0) ? nvb - 1 : nvb;
    long N = 1L << nbits;
    F->n = 0; for (int t = 0; t < PS; t++) F->fmin[t] = 1e30f;
    fib_reserve(F, 1024);
    enum { LEN = 3 * MAXN + 48 };
    int16_t t2[LEN] __attribute__((aligned(32))), tm2[LEN] __attribute__((aligned(32)));
    int t[MAXN]; memset(t, 0, sizeof(t));
    int vz[MAXV]; int vsg[MAXV];
    for (int j = 0; j < nvb; j++) { int o = vorb[b][j]; vz[j] = opos[kd][o][0]; vsg[j] = osg[kd][o][0]; t[vz[j]] = vsg[j]; }
    int *s = S[b];
    int P[MAXK], Qd[MAXK], pafS[MAXK]; int tau = 0;
    for (int i = 0; i < nI; i++) {
        int k = kI[i], pe = 0, ps_ = 0;
        for (int z = 0; z < n; z++) { int z2 = (z + k) % n; pe += t[z] * t[z2]; ps_ += s[z] * s[z2]; }
        P[i] = pe; pafS[i] = ps_;
    }
    for (int j = 0; j < nQ; j++) { int a = 0; for (int z = 0; z < n; z++) a += t[z] * wQ[b][z][j]; Qd[j] = a; }
    double ar[PS] __attribute__((aligned(32))), ai[PS] __attribute__((aligned(32)));
    for (int c = 0; c < PS; c++) { ar[c] = 0; ai[c] = 0; }
    for (int c = 0; c < nch; c++) {
        int psi = chpsi[c]; double re = 0, im = 0;
        for (int z = 0; z < n; z++) { re += s[z] * sre[psi][z] + t[z] * cgr[c][z]; im += s[z] * sim_[psi][z] + t[z] * cgi[c][z]; }
        ar[c] = re; ai[c] = im;
    }
    for (int z = 0; z < n; z++) tau += t[z];
    // vector state
    const int lenu = n + 40;
    for (int i = 0; i < lenu; i++) { t2[i] = (int16_t)t[i % n]; tm2[i] = (int16_t)t[((-i) % n + n) % n]; }
    int16_t Pa[16] __attribute__((aligned(32))), Qa[16] __attribute__((aligned(32))), WQ2[MAXN][16] __attribute__((aligned(32)));
    for (int a = 0; a < 16; a++) { Pa[a] = a < nI ? (int16_t)P[a] : 0; Qa[a] = a < nQ ? (int16_t)Qd[a] : 0; }
    for (int z = 0; z < n; z++) for (int j = 0; j < 16; j++) WQ2[z][j] = j < nQ ? (int16_t)(2 * wQ[b][z][j]) : 0;
    __m256i Pv = _mm256_load_si256((const __m256i *)Pa), Qv = _mm256_load_si256((const __m256i *)Qa);
    const int ncv = (nch + 3) >> 2;
    __m256d arv[PS / 4], aiv[PS / 4];
    for (int v = 0; v < ncv; v++) { arv[v] = _mm256_load_pd(ar + 4 * v); aiv[v] = _mm256_load_pd(ai + 4 * v); }
    const __m256d Tl = _mm256_set1_pd(T + 0.01);
    const __m256 DVv = _mm256_set1_ps((float)DV);
    __m256 fminv[PS / 8 + 1]; for (int v = 0; v < PS / 8 + 1; v++) fminv[v] = _mm256_set1_ps(1e30f);
    uint32_t mask = 0;
    const int *allow = allowT[b];
    for (long i = 0; i < N; i++) {
        if (i > 0) {
            int bit = __builtin_ctzl(i);
            mask ^= (1u << bit);
            int z = vz[bit]; int tv = t2[z];
            __m256i tvv = _mm256_set1_epi16((int16_t)tv);
            __m256i plus = _mm256_loadu_si256((const __m256i *)(t2 + z + 1));
            __m256i minus = _mm256_loadu_si256((const __m256i *)(tm2 + n - z + 1));
            __m256i d = _mm256_slli_epi16(_mm256_add_epi16(plus, minus), 1);
            Pv = _mm256_sub_epi16(Pv, _mm256_sign_epi16(d, tvv));
            Qv = _mm256_sub_epi16(Qv, _mm256_sign_epi16(_mm256_load_si256((const __m256i *)WQ2[z]), tvv));
            __m256d c2 = _mm256_set1_pd(2.0 * tv);
            const double *gr = CGRZ[z], *gi = CGIZ[z];
            for (int v = 0; v < ncv; v++) {
                arv[v] = _mm256_sub_pd(arv[v], _mm256_mul_pd(c2, _mm256_load_pd(gr + 4 * v)));
                aiv[v] = _mm256_sub_pd(aiv[v], _mm256_mul_pd(c2, _mm256_load_pd(gi + 4 * v)));
            }
            tau -= 2 * tv;
            int16_t nv = (int16_t)(-tv);
            for (int k = z; k < lenu; k += n) t2[k] = nv;
            int zz = z ? n - z : 0;
            for (int k = zz; k < lenu; k += n) tm2[k] = nv;
        }
        if (!allow[tau + MAXN]) continue;
        __m256d vv[PS / 4]; int bad = 0;
        for (int v = 0; v < ncv; v++) {
            vv[v] = _mm256_add_pd(_mm256_mul_pd(arv[v], arv[v]), _mm256_mul_pd(aiv[v], aiv[v]));
            bad |= _mm256_movemask_pd(_mm256_cmp_pd(vv[v], Tl, _CMP_GT_OQ));
        }
        if (bad) continue;
        if (F->n >= F->cap) fib_reserve(F, F->n + 1);
        {
            uint8_t *ps = F->psd + (size_t)F->n * PS;
            int32_t qi[PS] __attribute__((aligned(32)));
            for (int v = 0; v < ncv; v += 2) {
                __m128 f0 = _mm256_cvtpd_ps(vv[v]); __m128 f1 = (v + 1 < ncv) ? _mm256_cvtpd_ps(vv[v + 1]) : _mm_setzero_ps();
                __m256 fl = _mm256_insertf128_ps(_mm256_castps128_ps256(f0), f1, 1);
                fminv[v >> 1] = _mm256_min_ps(fminv[v >> 1], fl);
                _mm256_store_si256((__m256i *)(qi + 4 * v), _mm256_cvttps_epi32(_mm256_div_ps(fl, DVv)));
            }
            for (int c = 0; c < nch; c++) ps[c] = (uint8_t)qi[c];
            for (int c = nch; c < PS; c++) ps[c] = 0;
        }
        _mm256_store_si256((__m256i *)Pa, Pv); _mm256_store_si256((__m256i *)Qa, Qv);
        uint64_t hs = 0;
        for (int a = 0; a < nI; a++) hs += rc[a] * (uint64_t)(int64_t)(int16_t)(pafS[a] + q * (int)Pa[a]);
        for (int j = 0; j < nQ; j++) hs += rc[nI + j] * (uint64_t)(int64_t)(int16_t)(int)Qa[j];
        F->hs[F->n] = hs; F->mask[F->n] = mask; F->n++;
    }
    { float fm[PS + 8] __attribute__((aligned(32))); for (int v = 0; v < PS / 8 + 1; v++) _mm256_store_ps(fm + 8 * v, fminv[v]); for (int c = 0; c < nch; c++) F->fmin[c] = fm[c]; }
}
