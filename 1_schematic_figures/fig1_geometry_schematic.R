# Fig. 1: geometry of ING (schematic)
source("R/ing_model.R"); source("R/fig_theme.R")

MU   <- c(1, 2)
QZ_A <- Qz_uv(-0.212, -0.367)
QZ_I <- Qz_uv(0, -0.8)
NOCOST <- 1e9
RR <- c(1/16, 1/8, 1/4, 1/2, 1, 2, 4)
COMP_R <- cbind(RR, 1) / (1 + RR)
LIM <- c(0, 2.6)

mk <- function(Qz, eta = c(P = 0, C = 0), l2 = 0) list(mu = MU, Qz = Qz, eta = eta, l1 = 0, l2 = l2)
dense_all <- function(n = 3000) { th <- seq(0.002, pi / 2 - 0.002, length.out = n)
  d <- cbind(cos(th), sin(th)); d / rowSums(d) }
curve_of <- function(p, Fc = NOCOST, comp = dense_all()) rail_opt(p, comp, rep(Fc, nrow(comp))) * comp
frame <- function(lab, title) {
  blank(LIM, LIM, "Protein (P)", "Carbohydrate (C)"); panel_lab(lab); rails_draw(RR)
  mtext(title, 3, .3, cex = .72, col = MUTED) }
target <- function() { points(MU[1], MU[2], pch = 4, lwd = 2, col = INK)
  text(MU[1], MU[2], expression(bold(mu)), pos = 4, col = INK) }
curv_centre <- function(X, k) { d1 <- (X[k + 1, ] - X[k - 1, ]) / 2; d2 <- X[k + 1, ] - 2 * X[k, ] + X[k - 1, ]
  kap <- (d1[1] * d2[2] - d1[2] * d2[1]) / sum(d1^2)^1.5; N <- c(-d1[2], d1[1]) / sqrt(sum(d1^2)); X[k, ] + N / kap }
isoloss <- function(p, x0, col = REF) { g <- seq(-0.2, 3, length.out = 300)
  z <- outer(g, g, function(a, b) loss_R(cbind(a, b), p))
  contour(g, g, z, levels = loss_R(matrix(x0, 1), p), add = TRUE, drawlabels = FALSE, col = col, lwd = 1, lty = 2) }
inset_J <- function(p, comp1, box = c(1.45, 2.5, 0.12, 0.95), tmax = 4) {
  tt <- seq(0, tmax, length.out = 800); J <- loss_R(outer(tt, comp1), p)
  bk <- MU / comp1; bk <- bk[is.finite(bk) & bk < tmax]
  sx <- function(t) box[1] + (box[2] - box[1]) * t / tmax
  sy <- function(j) box[3] + (box[4] - box[3]) * (j - min(J)) / (quantile(J, .8) - min(J))
  rect(box[1], box[3], box[2], box[4], col = "white", border = GRID)
  keep <- sy(J) <= box[4]; lines(sx(tt)[keep], sy(J)[keep], col = ACC, lwd = 1.5)
  segments(sx(bk), box[3], sx(bk), box[4], col = ACC2, lty = 3)
  ts <- rail_opt(p, matrix(comp1, 1), NOCOST); points(sx(ts), sy(loss_R(matrix(ts * comp1, 1), p)), pch = 19, col = INK, cex = .7)
  text(box[1], box[4], "J(t)", adj = c(-0.1, 1.3), cex = .6, col = MUTED)
  text(box[2], box[3], "t", adj = c(1.4, -0.4), cex = .6, col = MUTED) }

fig("Fig1_geometry_schematic", 7.2, 10, function() {
  par(mfrow = c(3, 2), mar = c(3.6, 3.6, 2, 1), oma = c(0, 1.2, 0, 0))

  frame("a", "Asymmetry without interaction: curvature changes only")
  lines(curve_of(mk(diag(.5, 2))), col = MUTED, lwd = 1.3, lty = 2)
  lines(brk(curve_of(mk(diag(.5, 2), eta = c(P = log(.5), C = log(3)))), .05), col = ACC, lwd = 2)
  abline(v = MU[1], h = MU[2], col = GRID, lty = 3); target()
  legend("topright", bty = "n", cex = .65, lty = c(2, 1), lwd = c(1.3, 2), col = c(MUTED, ACC), text.col = INK,
         legend = c("symmetric", "protein excess cheap,\ncarbohydrate excess costly"))

  frame("b", "Asymmetry + interaction: corner solution")
  pe <- mk(QZ_I, eta = c(P = 0, C = log(3))); Xd <- curve_of(mk(QZ_I)); lines(Xd, col = MUTED, lwd = 1.3, lty = 2)
  Xe <- curve_of(pe); lines(brk(Xe, .05), col = ACC, lwd = 2)
  onC <- abs(Xe[, 2] - MU[2]) < 1e-9; lines(Xe[onC, ], col = ACC2, lwd = 4)
  abline(h = MU[2], col = GRID, lty = 3)
  cr <- which(onC)[round(sum(onC) / 2)]; c1 <- dense_all()[cr, ]
  segments(0, 0, 4 * c1[1], 4 * c1[2], col = INK, lwd = 1.2); points(Xe[cr, 1], Xe[cr, 2], pch = 19, col = INK, cex = .8)
  inset_J(pe, c1, tmax = 2 * MU[2] / c1[2]); target()

  frame("c", "Asymmetry + interaction: abrupt switch")
  pf <- mk(QZ_I, eta = c(P = log(.5), C = 0)); lines(Xd, col = MUTED, lwd = 1.3, lty = 2)
  Xf <- curve_of(pf); lines(brk(Xf, .05), col = ACC, lwd = 2)
  st <- sqrt(rowSums(diff(Xf)^2)); k <- which.max(st)
  segments(Xf[k, 1], Xf[k, 2], Xf[k + 1, 1], Xf[k + 1, 2], col = ACC2, lty = 2, lwd = 1.2)
  points(rbind(Xf[k, ], Xf[k + 1, ]), pch = 21, bg = ACC2, col = "white", cex = 1.1)
  c1 <- dense_all()[k, ]; segments(0, 0, 4 * c1[1], 4 * c1[2], col = INK, lwd = 1.2)
  abline(v = MU[1], col = GRID, lty = 3)
  inset_J(pf, c1, box = c(1.55, 2.55, 0.1, 0.9), tmax = 2.5 * MU[1] / c1[1]); target()
  legend(1.45, 2.62, bty = "n", cex = .65, lty = c(2, 1), lwd = c(1.3, 2), col = c(MUTED, ACC), text.col = INK,
         legend = c("symmetric, same Q", "asymmetric"))

  p <- mk(diag(.5, 2)); frame("d", "Isotropic: circle; evolute is a single point")
  ec <- effective_conic(p, NOCOST); lines(clip_pos(ellipse_pts(ec$Qe, ec$centre)), col = ACC, lwd = 2)
  X <- rail_opt(p, COMP_R, rep(NOCOST, 7)) * COMP_R
  segments(X[, 1], X[, 2], ec$centre[1], ec$centre[2], col = ACC2, lwd = .7)
  hl <- 3; segments(0, 0, 4 * COMP_R[hl, 1], 4 * COMP_R[hl, 2], col = MUTED, lwd = 1.4)
  isoloss(p, X[hl, ]); segments(X[hl, 1], X[hl, 2], MU[1], MU[2], col = INK, lty = 3)
  points(X, pch = 21, bg = "white", col = ACC, lwd = 1.4)
  points(ec$centre[1], ec$centre[2], pch = 19, col = ACC2, cex = 1.2)
  text(ec$centre[1], ec$centre[2], expression(E == mu / 2), pos = 2, col = ACC2, cex = .8); target()

  p <- mk(QZ_A); frame("e", expression("Anisotropic: ellipse centred at " * mu / 2 * "; extended evolute"))
  ec <- effective_conic(p, NOCOST); E <- ellipse_pts(ec$Qe, ec$centre, 6000)
  lines(clip_pos(E), col = ACC, lwd = 2); lines(ellipse_evolute(ec$Qe, ec$centre), col = ACC2, lwd = 1.2)
  for (k in round(seq(300, 5700, length.out = 9))) { if (any(E[k, ] < 0)) next
    cc <- curv_centre(E, k); segments(E[k, 1], E[k, 2], cc[1], cc[2], col = ACC2, lwd = .5) }
  X <- rail_opt(p, COMP_R, rep(NOCOST, 7)) * COMP_R
  hl <- 3; segments(0, 0, 4 * COMP_R[hl, 1], 4 * COMP_R[hl, 2], col = MUTED, lwd = 1.4)
  isoloss(p, X[hl, ]); segments(X[hl, 1], X[hl, 2], MU[1], MU[2], col = INK, lty = 3)
  points(X, pch = 21, bg = "white", col = ACC, lwd = 1.4)
  points(ec$centre[1], ec$centre[2], pch = 19, col = ACC2, cex = .9); target()

  p <- mk(QZ_A, l2 = .25); frame("f", "Volume cost: the conic deforms; apparent target moves")
  Fs <- c(.5, 1, 2, 8); traj <- NULL
  for (k in seq_along(Fs)) { ec <- effective_conic(p, Fs[k]); lines(clip_pos(ellipse_pts(ec$Qe, ec$centre)), col = SEQ[k], lwd = 1.8)
    traj <- rbind(traj, ec$mu_eff) }
  traj <- rbind(traj, MU); lines(traj, col = INK, lty = 2)
  points(traj[1:4, ], pch = 23, bg = SEQ, col = "white", cex = 1.2); target()
  legend("topright", bty = "n", cex = .65, lwd = 1.8, col = SEQ, text.col = INK,
         legend = paste0("concentration ", Fs, c(" (most dilute)", "", "", "")))

})
cat("Fig. 1 written to figures/main\n")
