# Fig. 5: specialist-generalist caterpillars
source("R/ing_model.R"); source("R/fig_theme.R")
L <- readRDS("out/larvae_fits.rds"); res <- L$res; SPP <- L$SPP
SPP <- SPP[SPP$diet_breadth != "hybrid", ]
prof_dros <- readRDS("out/profile_ratio.rds")
comp_of <- function(lbl) { v <- as.numeric(regmatches(lbl, gregexpr("[0-9]+", lbl))[[1]])[1:2]; v / sum(v) }
COLSP <- c(generalist = ACC, specialist = ACC2)
XLAB <- "Protein eaten (mg)"; YLAB <- "Carbohydrate eaten (mg)"
lin_t <- function(psi, lk, d) exp(lk) / drop(d %*% c(cos(psi), sin(psi)))

fig("Fig8_specialist_generalist", 7.2, 8.8, function() {
  par(mfrow = c(3, 2), mar = c(3.6, 4.2, 2, 1), oma = c(3.6, 0, 0, 0)); par(cex = .9)
  al <- seq(atan2(1, 5), atan2(5, 1), length.out = 300); dense <- cbind(sin(al), cos(al)); dense <- dense / rowSums(dense)
  for (k in seq_len(nrow(SPP))) { sp <- SPP$species[k]; r <- res[[sp]]; col <- COLSP[SPP$diet_breadth[k]]
    nc <- r$nc; ch <- r$ch; lim <- max(c(nc$P, nc$C, ch$P, ch$C)) * 1.05
    blank(c(0, lim), c(0, lim), XLAB, YLAB); panel_lab(letters[k])
    mtext(bquote(italic(.(SPP$short[k])) ~ "(" * .(SPP$diet_breadth[k]) * ")"), 3, .3, cex = .78)
    for (dl in unique(nc$diet)) { d <- comp_of(dl); segments(0, 0, 3 * lim * d[1], 3 * lim * d[2], col = GRID, lwd = .8) }
    for (f in strsplit(r$pair, " vs ")[[1]]) { d <- comp_of(f); segments(0, 0, 3 * lim * d[1], 3 * lim * d[2], col = MUTED, lty = 2, lwd = 1) }
    points(nc$P, nc$C, pch = 16, cex = .45, col = adjustcolor(INK, .25))
    TB <- apply(r$boot[, c("psi", "lk")], 1, function(b) lin_t(b[1], b[2], dense))
    qb <- apply(TB, 1, quantile, c(.025, .975))
    polygon(c(qb[1, ] * dense[, 1], rev(qb[2, ] * dense[, 1])), c(qb[1, ] * dense[, 2], rev(qb[2, ] * dense[, 2])),
            col = adjustcolor(col, .2), border = NA)
    a <- c(cos(r$lin$par[1]), sin(r$lin$par[1])); kk <- exp(r$lin$par[2])
    segments(0, kk / a[2], kk / a[1], 0, col = adjustcolor(col, .6), lty = 3)
    Xl <- lin_t(r$lin$par[1], r$lin$par[2], dense) * dense; lines(Xl, col = col, lwd = 2)
    fd <- t(sapply(strsplit(r$pair, " vs ")[[1]], comp_of)); tt <- lin_t(r$lin$par[1], r$lin$par[2], fd)
    segments(tt[1] * fd[1, 1], tt[1] * fd[1, 2], tt[2] * fd[2, 1], tt[2] * fd[2, 2], col = adjustcolor(INK, .7), lwd = 3.5, lend = 1)
    for (i in seq_along(r$dat$diets)) { y <- r$dat$lt[r$dat$trt == i]; d <- r$dat$comp[i, ]
      m <- mean(y); se <- sd(y) / sqrt(length(y))
      segments(exp(m - 1.96 * se) * d[1], exp(m - 1.96 * se) * d[2], exp(m + 1.96 * se) * d[1], exp(m + 1.96 * se) * d[2], col = INK, lwd = 1.2)
      points(exp(m) * d[1], exp(m) * d[2], pch = 21, bg = "white", col = INK, cex = 1.1, lwd = 1.3) }
    lP <- log(ch$P); lC <- log(ch$C); mP <- mean(lP); mC <- mean(lC); sP <- sd(lP) / sqrt(nrow(ch)); sC <- sd(lC) / sqrt(nrow(ch))
    segments(exp(mP - 1.96 * sP), exp(mC), exp(mP + 1.96 * sP), exp(mC), col = INK, lwd = 1.2)
    segments(exp(mP), exp(mC - 1.96 * sC), exp(mP), exp(mC + 1.96 * sC), col = INK, lwd = 1.2)
    points(exp(mP), exp(mC), pch = 23, bg = col, col = INK, cex = 1.5, lwd = 1.2)
    text(lim * .98, lim * .97, sprintf("w = %.2f", r$lin$w), adj = c(1, 1), cex = .75, col = col)
  }
  par(mar = c(3.6, 4.2, 2, 1))
  plot(NA, xlim = c(10, 80), ylim = c(0, 25), xlab = expression("Direction of the target, " * theta * " (degrees)"),
       ylab = expression(Delta * " negative log-likelihood"), bty = "l")
  abline(h = qchisq(.95, 1) / 2, col = GRID, lty = 2)
  s <- prof_dros[prof_dros$model == ING_BEST & !prof_dros$joint, ]; s <- s[order(s$r), ]
  lines(atan(s$r) * 180 / pi, s$dnll, col = INK, lwd = 2)
  for (k in seq_len(nrow(SPP))) { p <- res[[SPP$species[k]]]$profile
    lines(p$angle, p$dnll + (k - 1) * .25, col = COLSP[SPP$diet_breadth[k]], lwd = 1.6, lty = c(1, 1, 2, 2)[k]) }
  text(66, 14, "Drosophila\n(curved array)", font = 3, cex = .7, col = INK, adj = 1)
  text(12, 2.6, "caterpillars (straight arrays)", cex = .7, col = MUTED, adj = 0)
  panel_lab("e"); mtext("Is the target located?", 3, .3, cex = .78)
  w <- sapply(SPP$species, function(sp) res[[sp]]$lin$w)
  ci <- sapply(SPP$species, function(sp) quantile(res[[sp]]$boot[, "w"], c(.025, .975)))
  yy <- c(4, 3, 1.6, .6); par(mar = c(3.6, 6.4, 2, 1))
  plot(NA, xlim = c(.4, 3), ylim = c(0, 4.8), log = "x", yaxt = "n", ylab = "", bty = "l",
       xlab = expression("Relative weight of protein, " * w == a[P] / a[C]))
  abline(v = 1, col = GRID, lty = 2)
  axis(2, at = yy, labels = SPP$short, las = 1, font = 3, cex.axis = .78, tick = FALSE)
  segments(ci[1, ], yy, ci[2, ], yy, col = COLSP[SPP$diet_breadth], lwd = 2)
  points(w, yy, pch = 21, bg = COLSP[SPP$diet_breadth], col = "white", cex = 1.6)
  text(.42, c(4.6, 2.2), c("Heliothis", "Spodoptera"), adj = 0, font = 3, cex = .75, col = MUTED)
  text(1, 0.05, "equal weight\n(total P + C)", cex = .6, col = MUTED, adj = c(-.05, 0))
  panel_lab("f"); mtext("What is regulated", 3, .3, cex = .78)
  legend_below(ncol = 3, cex = .7,
    legend = c("larvae", "diet mean (95% CI)", "linear rule, 95% band", "foods of the choice test",
               "predicted choice (line within the foods)", "observed choice (95% CI)", "generalist", "specialist"),
    pch = c(16, 21, NA, NA, NA, 23, 15, 15), lty = c(NA, NA, 1, 2, 1, NA, NA, NA), lwd = c(NA, 1.3, 2, 1, 3.5, NA, NA, NA),
    col = c(adjustcolor(INK, .4), INK, MUTED, MUTED, INK, INK, ACC, ACC2), pt.bg = c(NA, "white", NA, NA, NA, "white", NA, NA),
    pt.cex = c(.8, 1, 1, 1, 1, 1.2, 1.4, 1.4))
})
cat("Figure 8 written\n")
