# Fig. S14: construction of the evolute (schematic)
INK <- "#1f2937"; MUTED <- "#6b7280"; GRID <- "#d1d5db"
ACC <- "#2563eb"; ACC2 <- "#c2410c"; REF <- "#9ca3af"
LEE_RATIOS <- c(0, 1/16, 1/8, 1/4, 1/2, 1, 1.875)

regulation <- function(mu, Q, w_excess = c(1, 1)) list(mu = mu, Q = Q, w = w_excess)
metric_of <- function(reg, x) {
  s <- ifelse(x > reg$mu, reg$w, 1); S <- diag(s, 2); S %*% reg$Q %*% S }

rail_optimum <- function(reg, d) {
  bk <- sort(unique(c(0, (reg$mu / d)[d > 0], 1e6)))
  best <- list(f = Inf)
  for (k in seq_len(length(bk) - 1)) {
    lo <- bk[k]; hi <- bk[k + 1]; M <- metric_of(reg, ((lo + hi) / 2) * d)
    a <- drop(t(d) %*% M %*% d); b <- drop(t(d) %*% M %*% reg$mu)
    t0 <- min(max(b / a, lo), hi)
    f <- drop(t(t0 * d - reg$mu) %*% M %*% (t0 * d - reg$mu))
    if (f < best$f - 1e-12) best <- list(f = f, t = t0, corner = (b / a < lo || b / a > hi))
  }
  list(x = best$t * d, corner = best$corner)
}
rail_dir <- function(r) c(r, 1) / (1 + r)
rule_of_compromise <- function(reg, n = 1500, ratios = NULL) {
  if (is.null(ratios)) { th <- seq(1e-4, pi / 2 - 1e-4, length.out = n); ratios <- cos(th) / sin(th) }
  out <- t(sapply(ratios, function(r) { o <- rail_optimum(reg, rail_dir(r)); c(o$x, o$corner) }))
  colnames(out) <- c("P", "C", "corner"); out
}

fit_conic <- function(X) {
  Z <- cbind(X[, 1]^2, X[, 1] * X[, 2], X[, 2]^2, X[, 1], X[, 2])
  cf <- svd(Z)$v[, 5]
  M <- matrix(c(cf[1], cf[2] / 2, cf[2] / 2, cf[3]), 2); g <- cf[4:5]
  if (sum(diag(M)) < 0) { M <- -M; g <- -g }
  centre <- -solve(M, g) / 2
  list(M = M, centre = centre, r2 = drop(t(centre) %*% M %*% centre))
}
conic_points <- function(cn, s) {
  L <- t(chol(cn$M / cn$r2))
  U <- rbind(cos(s), sin(s)); t(solve(t(L), U) + cn$centre)
}
conic_geom <- function(x, M, centre) {
  g <- drop(2 * M %*% (x - centre)); H <- 2 * M; ng <- sqrt(sum(g^2))
  kap <- (g[2]^2 * H[1, 1] - 2 * g[1] * g[2] * H[1, 2] + g[1]^2 * H[2, 2]) / ng^3
  N <- -g / ng; Tg <- c(-N[2], N[1])
  list(x = x, T = Tg, N = N, kappa = kap, rho = 1 / kap, E = x + N / kap)
}
geom_rows <- function(X, M_fun, centre) lapply(seq_len(nrow(X)), function(i)
  conic_geom(X[i, 1:2], M_fun(X[i, 1:2]), centre))
evolute_numeric <- function(X) {
  n <- nrow(X); i <- 2:(n - 1)
  d1 <- (X[i + 1, ] - X[i - 1, ]) / 2; d2 <- X[i + 1, ] - 2 * X[i, ] + X[i - 1, ]
  sp <- sqrt(rowSums(d1^2)); kap <- (d1[, 1] * d2[, 2] - d1[, 2] * d2[, 1]) / sp^3
  N <- cbind(-d1[, 2], d1[, 1]) / sp
  X[i, ] + N / kap
}

MU1 <- c(60, 100)
reg1 <- regulation(MU1, diag(2))
arr1 <- rule_of_compromise(reg1, ratios = LEE_RATIOS)[, 1:2]
cn1  <- fit_conic(arr1)
geo1 <- geom_rows(arr1, function(x) cn1$M, cn1$centre)
s_all <- seq(0, 2 * pi, length.out = 721)
curve1 <- conic_points(cn1, s_all)
E1 <- t(sapply(geom_rows(curve1, function(x) cn1$M, cn1$centre), `[[`, "E"))
MU1_hat <- 2 * cn1$centre

MU2 <- c(82.67, 129.47)
QZ2 <- matrix(c(0.7710, -0.0455, -0.0455, 0.2290), 2)
W2  <- exp(c(-1.2581, 0.6506))
Dinv <- diag(1 / MU2)
reg2  <- regulation(MU2, Dinv %*% QZ2 %*% Dinv, W2)
reg2s <- regulation(MU2, Dinv %*% QZ2 %*% Dinv)
arr2  <- rule_of_compromise(reg2, ratios = LEE_RATIOS)
rc2   <- rule_of_compromise(reg2)
smooth2 <- rc2[rc2[, "corner"] == 0, 1:2, drop = FALSE]
E2 <- t(sapply(geom_rows(smooth2, function(x) metric_of(reg2, x), MU2 / 2), `[[`, "E"))
branch <- apply(smooth2, 1, function(x) paste(x > MU2, collapse = ""))
E2b <- do.call(rbind, lapply(split(seq_len(nrow(E2)), cumsum(c(1, branch[-1] != branch[-length(branch)]))),
  function(ii) rbind(E2[ii, , drop = FALSE], NA)))
sym_conic <- list(M = reg2s$Q, centre = MU2 / 2, r2 = drop(t(MU2 / 2) %*% reg2s$Q %*% (MU2 / 2)))

cat("Circle example: fitted centre =", round(cn1$centre, 3), "; 2 x centre =", round(MU1_hat, 3),
    "; true mu =", MU1, "\n")
cat("  spread of the evolute (should be 0):", signif(max(abs(sweep(E1, 2, cn1$centre))), 3), "\n")
cat("  radius of curvature at the array points:", round(sapply(geo1, `[[`, "rho"), 3), "\n")
ell <- conic_points(sym_conic, seq(0.3, 1.2, length.out = 4001))
Enum <- evolute_numeric(ell)
Eana <- t(sapply(geom_rows(ell[2:4000, ], function(x) reg2s$Q, MU2 / 2), `[[`, "E"))
cat("Ellipse: max |numeric - analytic| centre of curvature:", signif(max(abs(Enum - Eana)), 3), "ug\n")
cat("Paper-like example: corner solutions on", sum(arr2[, "corner"]), "of 7 Lee rails;",
    "corners over", round(100 * mean(rc2[, "corner"]), 1), "% of all rails\n")

draw_rails <- function(ratios = LEE_RATIOS, col = GRID)
  for (r in ratios) segments(0, 0, 1e4 * r / (1 + r), 1e4 / (1 + r), col = col, lwd = .8)
frame <- function(lab, title, lim) {
  plot(NA, xlim = lim, ylim = lim, asp = 1, xlab = "Protein (P)", ylab = "Carbohydrate (C)",
       bty = "l", xaxs = "i", yaxs = "i")
  mtext(lab, 3, .6, adj = -0.2, font = 2, cex = .95); mtext(title, 3, .35, cex = .7, col = INK)
  draw_rails()
}
pts <- function(X, ...) points(X[, 1], X[, 2], pch = 21, bg = "white", col = INK, cex = 1.05, lwd = 1.3, ...)
mark_mu <- function(mu, half = TRUE, lab = TRUE) {
  points(mu[1], mu[2], pch = 4, lwd = 2.2, col = ACC2, cex = 1.1)
  if (lab) text(mu[1], mu[2], expression(bold(mu)), pos = 4, col = ACC2, cex = .9)
  if (half) { points(mu[1] / 2, mu[2] / 2, pch = 19, col = ACC2, cex = .9)
    if (lab) text(mu[1] / 2, mu[2] / 2, expression(mu/2), pos = 4, col = ACC2, cex = .8) }
}
in_range <- function(X) {
  a <- atan2(X[, 1], X[, 2]); a >= 0 & a <= atan(max(LEE_RATIOS)) + 1e-9 }
arc <- function(X) { X <- X[in_range(X), , drop = FALSE]
  X[order(atan2(X[, 1], X[, 2])), , drop = FALSE] }

draw <- function() {
  par(mfrow = c(2, 3), mar = c(3.4, 3.4, 2.4, .8), oma = c(3.0, 0, 0, 0), las = 1, mgp = c(2.1, .55, 0), tcl = -.3,
      cex.axis = .75, cex.lab = .85, col.axis = INK, col.lab = INK, fg = INK)
  L1 <- c(0, 125)

  frame("a", "The array: optimal intake on each rail", L1)
  pts(arr1)
  text(1e4 * 0, 120, "0:1", pos = 4, cex = .6, col = MUTED)
  text(115 * 1.875 / 2.875 * 1.55, 115 / 2.875 * 1.55, "1.9:1", cex = .6, col = MUTED)

  frame("b", expression("The curve "*gamma*"(s): conic through the origin"), L1)
  lines(curve1, col = ACC, lwd = 1, lty = 2); lines(arc(curve1), col = ACC, lwd = 2.2)
  pts(arr1); points(0, 0, pch = 19, cex = .7, col = INK)
  text(5, 118, expression(AP^2+BPC+CC^2+DP+EC==0), adj = 0, cex = .62, col = ACC)

  frame("c", expression("Tangents "*bold(T)*" and normals "*bold(N)), L1)
  lines(arc(curve1), col = ACC, lwd = 1.6)
  for (g in geo1) {
    segments(g$x[1] - 14 * g$T[1], g$x[2] - 14 * g$T[2], g$x[1] + 14 * g$T[1], g$x[2] + 14 * g$T[2],
             col = MUTED, lwd = 1.1)
    arrows(g$x[1], g$x[2], g$x[1] + 15 * g$N[1], g$x[2] + 15 * g$N[2], length = .05, col = INK, lwd = 1.1)
    segments(g$x[1] + 15 * g$N[1], g$x[2] + 15 * g$N[2], g$E[1], g$E[2], col = INK, lty = 3, lwd = .8)
  }
  pts(arr1)

  frame("d", expression("Centres of curvature "*bold(E)==gamma+kappa^-1*bold(N)), L1)
  lines(arc(curve1), col = ACC, lwd = 1.6)
  for (k in c(2, 6)) { g <- geo1[[k]]; th <- seq(0, 2 * pi, length.out = 300)
    lines(g$E[1] + g$rho * cos(th), g$E[2] + g$rho * sin(th), col = ACC2, lwd = .8, lty = 2) }
  for (g in geo1) arrows(g$x[1], g$x[2], g$E[1], g$E[2], length = .05, col = INK, lwd = .9)
  pts(arr1); points(cn1$centre[1], cn1$centre[2], pch = 19, col = ACC2)
  text(cn1$centre[1], cn1$centre[2] - 6, expression(paste("radius ", kappa^-1)), cex = .62, col = INK, adj = c(.2, 1))

  frame("e", expression("Evolute = one point, so "*mu==2%*%E), L1)
  lines(arc(curve1), col = ACC, lwd = 1.6); lines(curve1, col = ACC, lty = 2)
  points(E1[, 1], E1[, 2], pch = 19, col = ACC2, cex = 1.1)
  arrows(0, 0, MU1_hat[1], MU1_hat[2], length = .07, col = ACC2, lwd = 1.3)
  pts(arr1); mark_mu(MU1)

  L2 <- c(0, 175)
  frame("f", "As fitted: anisotropic, asymmetric", L2)
  abline(v = MU2[1], h = MU2[2], col = GRID, lty = 3)
  text(MU2[1], L2[2] * .985, "P at target", cex = .55, col = MUTED, adj = c(-.05, 1))
  text(L2[2] * .985, MU2[2], "C at target", cex = .55, col = MUTED, adj = c(1, -.4))
  lines(rc2[in_range(rc2), 1:2], col = ACC, lwd = 2.2)
  lines(rc2[!in_range(rc2), 1:2], col = ACC, lwd = 1, lty = 2)
  lines(E2b, col = ACC2, lwd = 2)
  sel <- c(3, 5)
  for (k in sel) { x <- arr2[k, 1:2]; if (arr2[k, "corner"] == 1) next
    g <- conic_geom(x, metric_of(reg2, x), MU2 / 2); th <- seq(0, 2 * pi, length.out = 400)
    lines(g$E[1] + g$rho * cos(th), g$E[2] + g$rho * sin(th), col = ACC2, lwd = .7, lty = 2)
    arrows(x[1], x[2], g$E[1], g$E[2], length = .05, col = INK, lwd = .9) }
  pts(arr2[, 1:2]); mark_mu(MU2)
  op <- par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE); on.exit(par(op))
  plot.new()
  legend("bottom", bty = "n", xpd = NA, ncol = 4, cex = .72, text.col = INK,
         legend = c("rails (fixed P:C)", "intakes on the rails", "curve (fitted conic / rule of compromise)",
                    "tangent", "normal (extended)", "osculating circle", "evolute (centres of curvature)", "intake target"),
         lty = c(1, NA, 1, 1, 3, 2, 1, NA), lwd = c(1, NA, 2, 1.1, .8, .8, 2, NA), pch = c(NA, 21, NA, NA, NA, NA, NA, 4),
         col = c(GRID, INK, ACC, MUTED, INK, ACC2, ACC2, ACC2), pt.bg = "white", pt.lwd = c(1, 1.3, 1, 1, 1, 1, 1, 2.2))
}

OUT <- if (dir.exists("figures/main")) "figures/main/Fig2_evolute_construction" else "Fig2_evolute_construction"
for (dev in c("pdf", "png")) {
  if (dev == "pdf") pdf(paste0(OUT, ".pdf"), 7.4, 5.9, useDingbats = FALSE)
  else png(paste0(OUT, ".png"), 7.4, 5.9, units = "in", res = 300)
  draw(); dev.off()
}
cat("Written:", OUT, "(.pdf, .png)\n")
