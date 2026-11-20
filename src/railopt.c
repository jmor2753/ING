/* ============================================================================
 * src/railopt.c  --  C CODE: exact forward solvers (compiled once, called from R)
 * ----------------------------------------------------------------------------
 *  rail_opt_c    optimum total intake on fixed-composition rails (no-choice
 *                diets): enumerates the at most three intervals between target
 *                crossings and takes the best clipped stationary point
 *                (Methods: "Solution under asymmetric regulation"; Appendix S2).
 *  Both accept an intake capacity lam[2] (total volume, uL/day; INFINITY = none):
 *  on a rail it truncates t at vmax/g. In the choice arena the capacity is
 *  shared (VY + VS <= vmax, with the shared volume cost) or per food
 *  (VY, VS <= vmax, with the food-specific cost); its boundary lines join the
 *  enumeration and infeasible candidates are rejected.
 *  choice_opt_c  optimum yeast and sucrose volumes in the two-food choice
 *                arena: enumerates stationary points of the four quadratic
 *                pieces, of lines along the axes and kink lines, and their
 *                intersections (Methods: "The choice experiment").
 * Build:  R CMD SHLIB src/railopt.c     (run_all.R does this automatically)
 * Called by R/ing_model.R via .C(); identical R implementations (rail_opt_R,
 * choice_opt_R) are used if the library is not compiled.
 * ========================================================================== */
#include <math.h>
/* Exact optimum on fixed-composition rails for the ING objective
   (piecewise-quadratic regulatory loss in proportional coordinates
    + lambda1 V + lambda2 V^2, V = (w1 c1 + w2 c2) t / F). */
void rail_opt_c(int *n, double *c1v, double *c2v, double *Fv, double *Qz,
                double *mu, double *kap, double *lam, double *w, double *out) {
  double q11 = Qz[0], q12 = Qz[1], q22 = Qz[3];
  double m1 = mu[0], m2 = mu[1], k1 = kap[0], k2 = kap[1], l1 = lam[0], l2 = lam[1];
  double vmax = lam[2];                      /* intake capacity (uL/day); INFINITY = none */
  for (int i = 0; i < *n; i++) {
    double c1 = c1v[i], c2 = c2v[i], g = (w[0]*c1 + w[1]*c2) / Fv[i];
    double tcap = isfinite(vmax) ? vmax / g : INFINITY;   /* V <= vmax  <=>  t <= vmax/g */
    double b[3]; int nb = 0; b[nb++] = 0.0;
    if (c1 > 0) b[nb++] = m1 / c1;
    if (c2 > 0) b[nb++] = m2 / c2;
    for (int a = 1; a < nb; a++) for (int c = a; c > 0 && b[c] < b[c-1]; c--) { double z=b[c]; b[c]=b[c-1]; b[c-1]=z; }
    double best_t = 0, best_J = INFINITY;
    for (int k = 0; k < nb; k++) {
      double lo = b[k], hi = (k + 1 < nb) ? b[k+1] : INFINITY;
      if (lo >= tcap) continue;
      if (hi > tcap) hi = tcap;
      if (hi <= lo) continue;
      double mid = isfinite(hi) ? 0.5*(lo+hi) : lo + 1.0;
      double s1 = (mid*c1 > m1) ? k1 : 1.0, s2 = (mid*c2 > m2) ? k2 : 1.0;
      double a1 = s1*c1/m1, a2 = s2*c2/m2;
      double A2 = q11*a1*a1 + 2*q12*a1*a2 + q22*a2*a2 + l2*g*g;
      double A1 = -2*(q11*a1*s1 + q12*(a1*s2 + a2*s1) + q22*a2*s2) + l1*g;
      double A0 = q11*s1*s1 + 2*q12*s1*s2 + q22*s2*s2;
      double tt = lo;
      if (A2 > 0) { tt = -A1/(2*A2); if (tt < lo) tt = lo; if (tt > hi) tt = hi; }
      if (!isfinite(tt)) continue;
      double J = A2*tt*tt + A1*tt + A0;
      if (J < best_J) { best_J = J; best_t = tt; }
    }
    out[i] = best_t;
  }
}

/* ---- exact two-food choice solver (yeast + sucrose) --------------------- */
static int SEP = 0;   /* 1: separable per-food volume cost l2 (VY^2 + VS^2) */
static double VMAX = INFINITY;   /* intake capacity on total volume VY + VS */
static double Jval(double vy, double vs, double *aY, double FS, double *Qz,
                   double *mu, double *kap, double l1, double l2) {
  double x1 = vy*aY[0], x2 = vy*aY[1] + vs*FS;
  double d1 = (x1 - mu[0]) / mu[0], d2 = (x2 - mu[1]) / mu[1];
  if (d1 > 0) d1 *= kap[0];
  if (d2 > 0) d2 *= kap[1];
  double V = vy + vs, V2 = SEP ? (vy*vy + vs*vs) : V*V;
  return Qz[0]*d1*d1 + 2*Qz[1]*d1*d2 + Qz[3]*d2*d2 + l1*V + l2*V2;
}
static void consider(double vy, double vs, double *best, double *bv,
                     double *aY, double FS, double *Qz, double *mu, double *kap, double l1, double l2) {
  if (!isfinite(vy) || !isfinite(vs)) return;
  if (vy < 0) vy = 0;
  if (vs < 0) vs = 0;
  if (SEP ? (vy > VMAX * (1 + 1e-10) || vs > VMAX * (1 + 1e-10))   /* per-food capacity */
          : (vy + vs > VMAX * (1 + 1e-10))) return;                  /* shared capacity */
  double J = Jval(vy, vs, aY, FS, Qz, mu, kap, l1, l2);
  if (J < *best) { *best = J; bv[0] = vy; bv[1] = vs; }
}
void choice_opt_c(int *n, double *FYv, double *FSv, double *Qz, double *mu,
                  double *kap, double *lam, double *yeast, int *sep, double *out) {
  double l1 = lam[0], l2 = lam[1]; SEP = *sep; VMAX = lam[2];
  for (int i = 0; i < *n; i++) {
    double FY = FYv[i], FS = FSv[i];
    double aY[2] = {FY*yeast[0], FY*yeast[1]};
    double best = INFINITY, bv[2] = {0, 0};
    double vyP = mu[0] / aY[0];
    #define LC(vy) ((mu[1] - aY[1]*(vy)) / FS)
    consider(0, 0, &best, bv, aY, FS, Qz, mu, kap, l1, l2);
    consider(vyP, 0, &best, bv, aY, FS, Qz, mu, kap, l1, l2);
    consider(0, LC(0), &best, bv, aY, FS, Qz, mu, kap, l1, l2);
    consider(vyP, LC(vyP), &best, bv, aY, FS, Qz, mu, kap, l1, l2);
    consider(mu[1]/aY[1], 0, &best, bv, aY, FS, Qz, mu, kap, l1, l2);
    if (isfinite(VMAX) && !SEP) {   /* shared capacity line VY + VS = VMAX: ends and crossings */
      consider(0, VMAX, &best, bv, aY, FS, Qz, mu, kap, l1, l2);
      consider(VMAX, 0, &best, bv, aY, FS, Qz, mu, kap, l1, l2);
      consider(vyP, VMAX - vyP, &best, bv, aY, FS, Qz, mu, kap, l1, l2);
      double den = 1 - aY[1]/FS;
      if (fabs(den) > 1e-12) { double vy = (VMAX - mu[1]/FS)/den; consider(vy, VMAX - vy, &best, bv, aY, FS, Qz, mu, kap, l1, l2); }
    }
    if (isfinite(VMAX) && SEP) {    /* per-food capacity box VY, VS <= VMAX: corners and crossings */
      consider(VMAX, 0, &best, bv, aY, FS, Qz, mu, kap, l1, l2);
      consider(0, VMAX, &best, bv, aY, FS, Qz, mu, kap, l1, l2);
      consider(VMAX, VMAX, &best, bv, aY, FS, Qz, mu, kap, l1, l2);
      consider(VMAX, LC(VMAX), &best, bv, aY, FS, Qz, mu, kap, l1, l2);
      consider(vyP, VMAX, &best, bv, aY, FS, Qz, mu, kap, l1, l2);
      consider((mu[1] - FS*VMAX)/aY[1], VMAX, &best, bv, aY, FS, Qz, mu, kap, l1, l2);
    }
    for (int r1 = 0; r1 < 2; r1++) for (int r2 = 0; r2 < 2; r2++) {
      double s1 = r1 ? kap[0] : 1.0, s2 = r2 ? kap[1] : 1.0;
      /* M = diag(s/mu) A,  A = [[aY0, 0],[aY1, FS]] */
      double M11 = s1*aY[0]/mu[0], M12 = 0, M21 = s2*aY[1]/mu[1], M22 = s2*FS/mu[1];
      /* QM = Qz M */
      double QM11 = Qz[0]*M11 + Qz[2]*M21, QM12 = Qz[0]*M12 + Qz[2]*M22;
      double QM21 = Qz[1]*M11 + Qz[3]*M21, QM22 = Qz[1]*M12 + Qz[3]*M22;
      double H11 = M11*QM11 + M21*QM21 + l2, H12 = M11*QM12 + M21*QM22 + (SEP ? 0.0 : l2);
      double H22 = M12*QM12 + M22*QM22 + l2;
      /* g = -2 M' Qz s + l1 */
      double Qs1 = Qz[0]*s1 + Qz[2]*s2, Qs2 = Qz[1]*s1 + Qz[3]*s2;
      double g1 = -2*(M11*Qs1 + M21*Qs2) + l1, g2 = -2*(M12*Qs1 + M22*Qs2) + l1;
      double det = H11*H22 - H12*H12;
      if (det > 0) consider(-(H22*g1 - H12*g2)/(2*det), -(-H12*g1 + H11*g2)/(2*det), &best, bv, aY, FS, Qz, mu, kap, l1, l2);
      /* 1-D minima along lines v0 + s dv */
      double lines[6][4] = {{0,0,1,0},{0,0,0,1},{vyP,0,0,1},{0,LC(0),1,-aY[1]/FS},{0,VMAX,1,-1},{0,0,0,0}};
      int nL = 4;
      if (isfinite(VMAX) && !SEP) nL = 5;                              /* VY + VS = VMAX */
      if (isfinite(VMAX) && SEP) {                                     /* VY = VMAX and VS = VMAX */
        lines[4][0] = VMAX; lines[4][1] = 0; lines[4][2] = 0; lines[4][3] = 1;
        lines[5][0] = 0; lines[5][1] = VMAX; lines[5][2] = 1; lines[5][3] = 0; nL = 6; }
      for (int L = 0; L < nL; L++) {
        double v0a = lines[L][0], v0b = lines[L][1], da = lines[L][2], db = lines[L][3];
        double a = H11*da*da + 2*H12*da*db + H22*db*db;
        double b = 2*(H11*da*v0a + H12*(da*v0b + db*v0a) + H22*db*v0b) + g1*da + g2*db;
        if (a > 0) { double s = -b/(2*a); consider(v0a + s*da, v0b + s*db, &best, bv, aY, FS, Qz, mu, kap, l1, l2); }
      }
    }
    out[4*i] = bv[0]*aY[0]; out[4*i+1] = bv[0]*aY[1] + bv[1]*FS; out[4*i+2] = bv[0]; out[4*i+3] = bv[1];
  }
}

/* ---------------------------------------------------------------------------
 * rail_opt_ext_c: forward solver for the EXTENDED models (5_supplementary/
 * analysis/S07). Regulatory deviations on a Box-Cox scale,
 *   z_j = ((x_j/mu_j)^beta - 1)/beta   (beta = 1: proportional, beta -> 0: log),
 * excess rescaled by kap_j, R = delta' Qz delta, plus a diet-specific linear
 * volume term l1e_i (= lambda1 - gamma h_i, sucrose palatability) and lambda2.
 * J(t) is no longer piecewise quadratic, so the minimum over t is found on a
 * dense grid (including the breakpoints t = mu_j/c_j) followed by golden-section
 * refinement of the best bracket. Checked against the exact solver at beta = 1.
 * l2 = (lambda2, Vmax): the search is truncated at the intake capacity.
 * ------------------------------------------------------------------------- */
static double Jext(double t, double c1, double c2, double g, double l1e, double *Qz,
                   double *mu, double *kap, double l2, double beta) {
  double r1 = t*c1/mu[0], r2 = t*c2/mu[1], z1, z2;
  if (fabs(beta - 1.0) < 1e-12) { z1 = r1 - 1; z2 = r2 - 1; }
  else { z1 = (pow(r1, beta) - 1)/beta; z2 = (pow(r2, beta) - 1)/beta; }
  if (z1 > 0) z1 *= kap[0];
  if (z2 > 0) z2 *= kap[1];
  return Qz[0]*z1*z1 + 2*Qz[2]*z1*z2 + Qz[3]*z2*z2 + l1e*g*t + l2*g*g*t*t;  /* l2 = lambda2 */
}
void rail_opt_ext_c(int *n, double *c1v, double *c2v, double *gv, double *l1ev,
                    double *Qz, double *mu, double *kap, double *l2, double *beta,
                    int *ngrid, double *out) {
  const double gr = 0.6180339887498949;
  for (int i = 0; i < *n; i++) {
    double c1 = c1v[i], c2 = c2v[i], g = gv[i], l1e = l1ev[i];
    double tmax = 0;
    if (c1 > 0 && mu[0]/c1 > tmax) tmax = mu[0]/c1;
    if (c2 > 0 && mu[1]/c2 > tmax) tmax = mu[1]/c2;
    tmax *= 6.0;
    if (isfinite(l2[1]) && l2[1] / g < tmax) tmax = l2[1] / g;     /* intake capacity: V <= Vmax */
    int N = *ngrid; double h = tmax/(N - 1), bestJ = INFINITY, bestt = 0; int bk = 0;
    for (int k = 0; k < N; k++) {
      double J = Jext(k*h, c1, c2, g, l1e, Qz, mu, kap, l2[0], *beta);
      if (J < bestJ) { bestJ = J; bestt = k*h; bk = k; }
    }
    double lo = (bk > 0) ? (bk - 1)*h : 0, hi = (bk < N - 1) ? (bk + 1)*h : tmax;
    double a = lo, b = hi, x1 = b - gr*(b - a), x2 = a + gr*(b - a);
    double f1 = Jext(x1, c1, c2, g, l1e, Qz, mu, kap, l2[0], *beta), f2 = Jext(x2, c1, c2, g, l1e, Qz, mu, kap, l2[0], *beta);
    for (int it = 0; it < 100; it++) {
      if (f1 < f2) { b = x2; x2 = x1; f2 = f1; x1 = b - gr*(b - a); f1 = Jext(x1, c1, c2, g, l1e, Qz, mu, kap, l2[0], *beta); }
      else { a = x1; x1 = x2; f1 = f2; x2 = a + gr*(b - a); f2 = Jext(x2, c1, c2, g, l1e, Qz, mu, kap, l2[0], *beta); }
    }
    double tg = 0.5*(a + b), Jg = Jext(tg, c1, c2, g, l1e, Qz, mu, kap, l2[0], *beta);
    if (Jg < bestJ) { bestJ = Jg; bestt = tg; }
    /* exact corners: breakpoints */
    double bks[2] = {c1 > 0 ? mu[0]/c1 : -1, c2 > 0 ? mu[1]/c2 : -1};
    for (int q = 0; q < 2; q++) if (bks[q] > 0 && bks[q] <= tmax) {
      double Jb = Jext(bks[q], c1, c2, g, l1e, Qz, mu, kap, l2[0], *beta);
      if (Jb < bestJ) { bestJ = Jb; bestt = bks[q]; }
    }
    out[i] = bestt;
  }
}
