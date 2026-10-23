// achievable tau sums of a block: DP over variables with weight w = sum of sg over the orbit
static void tau_set(int b, char *set /* size 2*MAXN+1 */) {
    int kd = b == 0 ? 0 : 1; static char cur[2 * MAXN * 4 + 1]; int OFF = 2 * MAXN;
    memset(cur, 0, sizeof(cur)); cur[OFF] = 1;
    for (int j = 0; j < nvar[b]; j++) {
        int o = vorb[b][j]; int w = 0; for (int m = 0; m < omm[kd][o]; m++) w += osg[kd][o][m];
        static char nxt[2 * MAXN * 4 + 1]; memset(nxt, 0, sizeof(nxt));
        for (int x = 0; x <= 2 * OFF; x++) if (cur[x]) { if (x + w >= 0 && x + w <= 2 * OFF) nxt[x + w] = 1; if (x - w >= 0 && x - w <= 2 * OFF) nxt[x - w] = 1; }
        memcpy(cur, nxt, sizeof(cur));
    }
    memset(set, 0, 2 * MAXN + 1);
    for (int x = 0; x <= 2 * OFF; x++) if (cur[x]) { int tv = x - OFF; if (tv >= -MAXN && tv <= MAXN) set[tv + MAXN] = 1; }
}

int main(int argc, char **argv) {
    mallopt(M_MMAP_THRESHOLD, 1 << 30); mallopt(M_TRIM_THRESHOLD, 1 << 30); mallopt(M_TOP_PAD, 16 << 20);
    if (argc < 2) { fprintf(stderr, "usage: pl2 file [maxsol_print] [tlimit] [seed]\n"); return 1; }
    FILE *f = fopen(argv[1], "r"); if (!f) return 2;
    maxsol_print = argc > 2 ? atoi(argv[2]) : 5;
    double tlimit = argc > 3 ? atof(argv[3]) : 1e30;
    if (argc > 4) rs ^= (uint64_t)atol(argv[4]) * 0xD1B54A32D192ED03ULL;
    if (fscanf(f, "%d %d %d %d", &n, &q, &eps, &T) != 4) return 2;
    if (fscanf(f, "%d", &nI) != 1) return 2; for (int i = 0; i < nI; i++) if (fscanf(f, "%d", &kI[i]) != 1) return 2;
    if (fscanf(f, "%d", &nQ) != 1) return 2; for (int i = 0; i < nQ; i++) if (fscanf(f, "%d", &kQ[i]) != 1) return 2;
    for (int kd = 0; kd < 2; kd++) {
        if (fscanf(f, "%d", &nu[kd]) != 1) return 2;
        for (int o = 0; o < nu[kd]; o++) { if (fscanf(f, "%d", &omm[kd][o]) != 1) return 2; for (int m = 0; m < omm[kd][o]; m++) if (fscanf(f, "%d %d", &opos[kd][o][m], &osg[kd][o][m]) != 2) return 2; }
    }
    int nchl, chl[MAXN]; if (fscanf(f, "%d", &nchl) != 1) return 2; for (int i = 0; i < nchl; i++) if (fscanf(f, "%d", &chl[i]) != 1) return 2;
    int ninst; if (fscanf(f, "%d", &ninst) != 1) return 2;
    nk = nI + nQ; sqq = sqrt((double)q);
    if (eps != -1) { fprintf(stderr, "skf2: eps = -1 only\n"); return 4; }
    if (nchl > PS) { fprintf(stderr, "PS too small\n"); return 4; }
    // channels: the psi listed in the header (orbit representatives of psi under H_K; all psi = 1..n-1 for H_K = 1)
    nch = 0;
    for (int psi = 1; psi < n; psi++) for (int z = 0; z < n; z++) { double a = 2 * M_PI * ((psi * z) % n) / n; sre[psi][z] = cos(a); sim_[psi][z] = sin(a); }
    if (eps == -1) {
        for (int ii = 0; ii < nchl; ii++) { int psi = chl[ii]; chpsi[nch] = psi; for (int z = 0; z < n; z++) { // gamma = i sqrt(q):  gamma*(c + i s) = -sqrt(q) s + i sqrt(q) c
                cgr[nch][z] = -sqq * sim_[psi][z]; cgi[nch][z] = sqq * sre[psi][z]; } nch++; }
    } else {
        for (int psi = 1; psi <= (n - 1) / 2; psi++) for (int sg = 1; sg >= -1; sg -= 2) { chpsi[nch] = psi; for (int z = 0; z < n; z++) { cgr[nch][z] = sg * sqq * sre[psi][z]; cgi[nch][z] = sg * sqq * sim_[psi][z]; } nch++; }
    }
    DV = (T + 249) / 250; if (DV < 1) DV = 1;
    int nojoin = getenv("PL2F_NOJOIN") ? 1 : 0;
    fast_tables(); int use_fast = !(getenv("PL2F_NOFAST") && atoi(getenv("PL2F_NOFAST"))); if (getenv("PL2F_FASTMIN")) fast_min_vars = atoi(getenv("PL2F_FASTMIN"));
    double spc = getenv("PL2F_SPC") ? atof(getenv("PL2F_SPC")) : 1500.0;
    for (int j = 0; j < 2 * MAXK; j++) rc[j] = rnd() | 1;
    for (int inst = 0; inst < ninst; inst++) {
        for (int b = 0; b < 4; b++) for (int z = 0; z < n; z++) if (fscanf(f, "%d", &S[b][z]) != 1) return 2;
        for (int b = 0; b < 4; b++) { int kd = b == 0 ? 0 : 1; for (int o = 0; o < nu[kd]; o++) if (fscanf(f, "%d", &Dfl[b][o]) != 1) return 2; }
        double t0 = now();
        for (int b = 0; b < 4; b++) {
            int kd = b == 0 ? 0 : 1; nvar[b] = 0;
            for (int o = 0; o < nu[kd]; o++) if (!Dfl[b][o]) { if (nvar[b] >= MAXV) { fprintf(stderr, "too many vars\n"); return 5; } vorb[b][nvar[b]++] = o; }
            for (int z = 0; z < n; z++) for (int j = 0; j < nQ; j++) { int zp = (z + kQ[j]) % n, zm = ((z - kQ[j]) % n + n) % n; wQ[b][z][j] = S[b][zp] + eps * S[b][zm]; }
        }
        // order the variables so that the fixed one (last) is in the block with the most variables (C or D preferred)
        int fixk = 2; for (int b = 0; b < 4; b++) if (nvar[b] > nvar[fixk]) fixk = b;
        if (nvar[fixk] == 0) fixk = -1;
        // row-sum / tau filters
        int sg2 = 0; int sgm[4]; for (int b = 0; b < 4; b++) { int sr = 0; for (int z = 0; z < n; z++) sr += S[b][z]; sgm[b] = sr; sg2 += sr * sr; }
        static char tset[4][2 * MAXN + 1];
        for (int b = 0; b < 4; b++) tau_set(b, tset[b]);
        memset(allowT, 0, sizeof(allowT));
        {
            // left (A,B) sums, right (C,D)
            int rem0 = T - sg2;
            if (rem0 >= 0) {
                for (int ta = -MAXN; ta <= MAXN; ta++) if (tset[0][ta + MAXN])
                for (int tb = -MAXN; tb <= MAXN; tb++) if (tset[1][tb + MAXN]) {
                    int r1 = rem0 - q * (ta * ta + tb * tb); if (r1 < 0) continue;
                    int cr1 = sgm[0] * ta + sgm[1] * tb;
                    for (int tc = -MAXN; tc <= MAXN; tc++) if (tset[2][tc + MAXN]) {
                        int r2 = r1 - q * tc * tc; if (r2 < 0) continue;
                        for (int td = -MAXN; td <= MAXN; td++) if (tset[3][td + MAXN]) {
                            if (q * td * td != r2) continue;
                            if (eps == 1 && cr1 + sgm[2] * tc + sgm[3] * td != 0) continue;
                            allowT[0][ta + MAXN] = allowT[1][tb + MAXN] = allowT[2][tc + MAXN] = allowT[3][td + MAXN] = 1;
                        }
                    }
                }
            }
        }
        int anyok = 1; for (int b = 0; b < 4; b++) { int a = 0; for (int x = 0; x <= 2 * MAXN; x++) a |= allowT[b][x]; if (!a) anyok = 0; }
        Fib F[4]; memset(F, 0, sizeof(F));
        long long compl = 0, matches = 0, surv1 = 0, surv2 = 0; int timed_out = 0; double npairs1 = 0, npairs2 = 0;
        int big = 0; for (int b = 0; b < 4; b++) if (nvar[b] - (b == fixk ? 1 : 0) > 24) big = 1;
        if (!anyok) { printf("INST %d rowfilter\n", inst); fflush(stdout); continue; }
        if (big) { printf("INST %d TOOBIG nv=%d,%d,%d,%d\n", inst, nvar[0], nvar[1], nvar[2], nvar[3]); fflush(stdout); continue; }
        // move the fixed variable to the end of its block's list
        for (int b = 0; b < 4; b++) { if (use_fast && fast_ok(b)) build_fibre_fast(b, &F[b], b == fixk); else build_fibre(b, &F[b], b == fixk); }
        double tfib = now() - t0;
        int sx, sy, sz, sw;      // skew plug-in: pair the smallest fibre with the largest and the two middle ones (as gl2.c); minimises |X||Y| + |Z||W|
        { int ord[4] = {0, 1, 2, 3};
          for (int i = 0; i < 4; i++) for (int j = i + 1; j < 4; j++) if (F[ord[j]].n < F[ord[i]].n) { int t_ = ord[i]; ord[i] = ord[j]; ord[j] = t_; }
          sx = ord[0]; sy = ord[3]; sz = ord[1]; sw = ord[2]; }
        npairs1 = (double)F[sx].n * F[sy].n; npairs2 = (double)F[sz].n * F[sw].n;
        int emptyf = 0; for (int b = 0; b < 4; b++) if (F[b].n == 0) emptyf = 1;
        compl_g = 0; matches_g = 0; cur_inst = inst;
        if (!emptyf && !nojoin) {
            // per-channel limits (bytes): same quantisation as pl2.c
            uint8_t L1[PS], L2[PS];
            for (int c = 0; c < PS; c++) {
                if (c < nch) {
                    float l1 = T + 0.02f - F[sz].fmin[c] - F[sw].fmin[c], l2 = T + 0.02f - F[sx].fmin[c] - F[sy].fmin[c];
                    L1[c] = (uint8_t)(l1 < 0 ? 0 : (l1 / DV > 255 ? 255 : (int)(l1 / DV)));
                    L2[c] = (uint8_t)(l2 < 0 ? 0 : (l2 / DV > 255 ? 255 : (int)(l2 / DV)));
                } else { L1[c] = L2[c] = 255; }
            }
            double bigp = npairs1 > npairs2 ? npairs1 : npairs2;
            int lgM = 0; while ((bigp * 0.1 / (double)(1L << lgM)) > spc && lgM < 22) lgM++;
            if (getenv("PL2F_LGM")) lgM = atoi(getenv("PL2F_LGM"));
            uint64_t M = 1ull << lgM, mm = M - 1;
            Side SL, SR; memset(&SL, 0, sizeof(SL)); memset(&SR, 0, sizeof(SR));
            side_build(&SL, sx, sy, F, L1); side_build(&SR, sz, sw, F, L2);
            CB *cbL = calloc(M, sizeof(CB)), *cbR = calloc(M, sizeof(CB));
            surv1 = scan_side(&SL, mm, 0, cbL, t0, tlimit, &timed_out);
            if (!timed_out) surv2 = scan_side(&SR, mm, 1, cbR, t0, tlimit, &timed_out);
            if (!timed_out) join_classes(&SL, &SR, F, mm, M, cbL, cbR, t0, tlimit, &timed_out);
            cbs_free(cbL, M); cbs_free(cbR, M); side_free(&SL); side_free(&SR);
        }
        compl = compl_g; matches = matches_g;
        long long mult = (fixk >= 0) ? 2 : 1;
        total_compl += compl * mult;
        printf("INST %d nv=%d,%d,%d,%d fib=%ld,%ld,%ld,%ld pairs=%.3g,%.3g surv=%lld,%lld matches=%lld completions=%lld%s t=%.3f tfib=%.3f\n",
               inst, nvar[0], nvar[1], nvar[2], nvar[3], F[0].n, F[1].n, F[2].n, F[3].n, npairs1, npairs2, surv1, surv2, matches,
               compl * mult, timed_out ? " TIMEOUT" : "", now() - t0, tfib);
        fflush(stdout);
        for (int b = 0; b < 4; b++) fib_free(&F[b]);
    }
    printf("TOTAL completions %lld\n", total_compl);
    return 0;
}
