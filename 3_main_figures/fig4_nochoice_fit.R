# Fig. 2: no-choice arrays and the best model
source("R/ing_model.R"); source("R/fig_theme.R")
F0 <- readRDS("out/fits_nochoice.rds"); dat <- F0$dat; fits <- F0$fits
best <- fits[[ING_BEST]]; symm <- fits[[ING_SYMM]]
CB <- readRDS("out/cv_boot.rds"); BP <- boot_params(CB, best)
lee_ratios <- LEE_RATIOS

fig("Fig4_nochoice_fit", 7.2, 7.6, function() {
  par(mfrow = c(2, 2), mar = c(3.6, 3.8, 2, 1), oma = c(3.2, 0, 0, 0))
  al <- seq(0, atan(1.95), length.out = 600); dense <- cbind(sin(al), cos(al)); dense <- dense / rowSums(dense)
  for (k in seq_along(c(45, 90, 180, 360))) { F <- c(45, 90, 180, 360)[k]; kt <- which(dat$F == F)
    blank(c(0, 145), c(0, 215), expression("Protein eaten ("*mu*"g day"^-1*")"), expression("Carbohydrate eaten ("*mu*"g day"^-1*")"))
    panel_lab(letters[k]); rails_draw(lee_ratios, 400)
    segments(0, 0, 400 * .2, 400 * .8, col = "#9ca3af", lty = 2)
    obs <- exp(dat$ybar[kt]) * dat$comp[kt, ]
    se <- sqrt(dat$s2n[kt])
    lo <- exp(dat$ybar[kt] - 1.96 * se) * dat$comp[kt, ]; hi <- exp(dat$ybar[kt] + 1.96 * se) * dat$comp[kt, ]
    segments(lo[, 1], lo[, 2], hi[, 1], hi[, 2], col = INK, lwd = 1)
    TB <- sapply(BP, function(pb) rail_opt(pb, dense, rep(F, nrow(dense))))
    qb <- apply(TB, 1, quantile, c(.025, .975)); lo_b <- qb[1, ] * dense; hi_b <- qb[2, ] * dense
    polygon(c(lo_b[, 1], rev(hi_b[, 1])), c(lo_b[, 2], rev(hi_b[, 2])), col = adjustcolor(ACC, .18), border = NA)
    Xs <- rail_opt(symm$p, dense, rep(F, nrow(dense))) * dense; lines(Xs, col = MUTED, lty = 2, lwd = 1.3)
    Xb <- rail_opt(best$p, dense, rep(F, nrow(dense))) * dense; lines(brk(Xb), col = ACC, lwd = 2)
    pr <- best$pred[kt] * dat$comp[kt, ]; points(pr, pch = 21, bg = "white", col = ACC, lwd = 1.4, cex = 1.05)
    points(obs, pch = 19, col = INK, cex = .8)
    points(best$p$mu[1], best$p$mu[2], pch = 4, lwd = 2, col = INK)
    mtext(paste0("Food concentration ", F, " g/L"), 3, .3, cex = .78, col = INK)
  }
  legend_below(ncol = 4,
      legend = c("Observed mean (95% CI)", "ING prediction", "Best model", "95% bootstrap band", "Same model, symmetric",
                 "Latent intake target", "P:C 1:4 (self-selected, choice)"),
      pch = c(19, 21, NA, 15, NA, 4, NA), lty = c(NA, NA, 1, NA, 2, NA, 2), lwd = c(NA, 1.4, 2, NA, 1.3, 2, 1), pt.cex = c(1, 1, 1, 1.8, 1, 1, 1),
      col = c(INK, ACC, ACC, adjustcolor(ACC, .18), MUTED, INK, "#9ca3af"), pt.bg = c(NA, "white", NA, NA, NA, NA, NA))
})
cat("Figure 4 written\n")
