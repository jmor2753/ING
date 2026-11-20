# Fig. S16: performance surfaces and the ING target
source("R/ing_model.R"); source("R/fig_theme.R")
F0 <- readRDS("out/fits_nochoice.rds"); dat <- F0$dat; fits <- F0$fits
best <- fits[[ING_BEST]]; symm <- fits[[ING_SYMM]]
ch <- load_choice(6); CP180 <- mean(ch$C_per_P[ch$sucrose_gL == 180])
TN <- readRDS("out/traits_nutrigon.rds")

suppressMessages(library(nutrigon))
lab <- c(lifespan = "Lifespan", lifetimeegg = "Lifetime egg production", dailyeggs = "Egg production rate")
fig("Fig7_trait_surfaces", 7.2, 7.3, function() {
  par(mfrow = c(2, 2), mar = c(3.6, 3.8, 2, 1), oma = c(1.6, 0, 0, 0))
  pal <- SURF
  for (k in seq_along(lab)) { y <- names(lab)[k]; f <- TN$fits[[paste(y, "lm")]]
    g <- nutrigon_grid(f, n = 150, quiet = TRUE)
    xs <- sort(unique(g$protein_eaten)); ys <- sort(unique(g$carb_eaten))
    z <- matrix(NA, length(xs), length(ys)); z[cbind(match(g$protein_eaten, xs), match(g$carb_eaten, ys))] <- g$predicted
    image(xs, ys, z, col = pal, xlim = c(0, 185), ylim = c(0, 305), useRaster = TRUE, bty = "l",
          xlab = expression("Protein eaten ("*mu*"g day"^-1*")"), ylab = expression("Carbohydrate eaten ("*mu*"g day"^-1*")"))
    contour(xs, ys, z, add = TRUE, nlevels = 6, col = "#ffffff", lwd = .6, labcex = .5)
    lines(nutrigon_hull(f), col = MUTED, lwd = .8)
    rp <- nutrigon_region_polygon(f, quiet = TRUE); polygon(rp[[1]], rp[[2]], border = ACC2, lwd = 1.6, lty = 1)
    o <- f$optimum; points(o$protein_eaten, o$carb_eaten, pch = 24, bg = ACC2, col = "white", cex = 1.4)
    segments(0, 0, 400 * cos(o$theta * pi / 180), 400 * sin(o$theta * pi / 180), col = ACC2, lty = 3)
    points(best$p$mu[1], best$p$mu[2], pch = 4, lwd = 2, col = INK)
    panel_lab(letters[k]); mtext(lab[k], 3, .3, cex = .78) }
  B <- TN$boot; est <- TN$est
  nm <- c(ing = "ING latent target", lifespan = "Lifespan", lifetimeegg = "Lifetime eggs", dailyeggs = "Egg rate")
  yy <- 4:1; ci <- t(apply(B[, names(nm)], 2, quantile, c(.025, .975)))
  par(mar = c(3.6, 3.8, 3.6, 1)); plot(NA, xlim = c(40, 95), ylim = c(.4, 5.2), yaxt = "n", ylab = "", bty = "l",
       xlab = expression("Angle from protein axis, " * theta * " (degrees)"))
  axis(2, at = yy, labels = nm, las = 1, cex.axis = .72, tick = FALSE)
  rat <- c(1, 2, 4, 8, 16); axis(3, at = atan(rat) * 180 / pi, labels = paste0("1:", rat), cex.axis = .65, col.axis = MUTED, tcl = -.2, mgp = c(2, .3, 0))
  abline(v = atan(CP180) * 180 / pi, col = ACC2, lty = 2); text(atan(CP180) * 180 / pi, 5.05, "choice (sucrose 180)", col = ACC2, cex = .62, pos = 2)
  set.seed(1)
  for (k in seq_along(nm)) { v <- B[, names(nm)[k]]
    points(v, yy[k] + runif(length(v), -.22, .22), pch = 16, cex = .45, col = adjustcolor(if (k == 1) INK else ACC, .35))
    segments(ci[k, 1], yy[k] - .32, ci[k, 2], yy[k] - .32, col = if (k == 1) INK else ACC, lwd = 1.6) }
  points(est[names(nm)], yy, pch = c(4, 24, 24, 24), bg = c(NA, ACC2, ACC2, ACC2), col = c(INK, "white", "white", "white"), cex = 1.4, lwd = 2)
  mtext("d", 3, 2.0, adj = -0.18, font = 2, cex = .95); mtext("Optima as angles (bootstrap draws, 95% range)", 3, 2.1, cex = .78)
  legend_below(ncol = 4, pch = c(24, NA, 4, NA), lty = c(NA, 1, NA, 3), col = c("white", ACC2, INK, ACC2), pt.bg = c(ACC2, NA, NA, NA),
               lwd = c(NA, 1.6, 2, 1), pt.cex = 1.2,
               legend = c("trait optimum (nutrigon)", "peak region (top 5%)", "ING latent target", "ray through the optimum"))
})

cat("Figure 7 written\n")
