# Supplementary figures
source("R/ing_model.R"); source("R/fig_theme.R"); suppressMessages(library(nutrigon)); suppressMessages(library(mgcv))
DIR <- "figures/supplementary"
F0 <- readRDS("out/fits_nochoice.rds"); dat <- F0$dat; fits <- F0$fits
best <- fits[[ING_BEST]]; symm <- fits[[ING_SYMM]]
cv <- read.csv("out/table_cv.csv")
short <- function(n) sub("\\+NoAcq", " (no cost)", gsub("\\+", " + ", n))

fig("FigS3_model_comparison", 7.2, 4.6, dir = DIR, function() {
  tab <- data.frame(model = names(fits), AICc = sapply(fits, `[[`, "AICc"))
  tab <- merge(tab, cv[, setdiff(names(cv), "AICc")], by = "model"); tab <- tab[order(tab$AICc, decreasing = TRUE), ]
  tab <- tab[!grepl("NoAcq", tab$model), ]
  n <- nrow(tab); yy <- seq_len(n)
  par(mfrow = c(1, 2), mar = c(3.6, 9.5, 2, 1), oma = c(1.4, 0, 0, 0))
  plot(tab$AICc - min(tab$AICc), yy, pch = 19, col = ACC, yaxt = "n", ylab = "", xlab = expression(Delta * "AICc"), bty = "l", ylim = c(.5, n + 2))
  axis(2, yy, short(tab$model), las = 1, cex.axis = .68, tick = FALSE); abline(v = 2, col = GRID, lty = 2)
  mtext("a", side = 3, adj = -0.9, line = 0.6, font = 2, cex = .95); panel_title("In-sample support")
  par(mar = c(3.6, 1, 2, 1))
  plot(NA, xlim = range(c(tab$pct_rail, tab$pct_conc)) + c(-1, 1), ylim = c(.5, n + 2), yaxt = "n", ylab = "",
       xlab = "Prediction error (%)", bty = "l")
  segments(tab$pct_rail, yy, tab$pct_conc, yy, col = GRID)
  points(tab$pct_rail, yy, pch = 19, col = ACC); points(tab$pct_conc, yy, pch = 21, bg = "white", col = INK)
  mtext("b", side = 3, adj = 0, line = 0.6, font = 2, cex = .95); panel_title("Out-of-sample prediction")
  legend_below(ncol = 2, pch = c(19, 21), col = c(ACC, INK), pt.bg = c(NA, "white"),
               legend = c("(b) leave one rail out", "(b) leave one concentration out"))
})

fig("FigS4_residuals", 7.2, 3.6, dir = DIR, function() {
  res <- dat$ybar - log(best$pred); se <- sqrt(dat$s2n)
  rails <- c("(0:1)", "(1:16)", "(1:8)", "(1:4)", "(1:2)", "(1:1)", "(2:1)")
  rl <- c("0:1", "1:16", "1:8", "1:4", "1:2", "1:1", "1.9:1")
  Fs <- c(45, 90, 180, 360); off <- c(-.24, -.08, .08, .24)
  par(mar = c(3.6, 4, 2, 1))
  plot(NA, xlim = c(.5, 7.5), ylim = range(c(res - 2 * se, res + 2 * se)), xaxt = "n", bty = "l",
       xlab = "Nutritional rail (P:C)", ylab = "Observed - predicted log intake")
  axis(1, 1:7, rl); abline(h = 0, col = GRID, lty = 2)
  for (k in 1:4) { ii <- which(dat$F == Fs[k]); x <- match(dat$tr$Ratio[ii], rails) + off[k]
    segments(x, res[ii] - 1.96 * se[ii], x, res[ii] + 1.96 * se[ii], col = SEQ[k])
    points(x, res[ii], pch = 19, col = SEQ[k], cex = .9) }
  legend("topright", bty = "n", cex = .68, pch = 19, col = SEQ, legend = paste(Fs, "g/L"), horiz = TRUE, text.col = INK)
  panel_title(paste0("Best model (", ING_BEST, "); bars: 95% sampling interval"))
})

CB <- readRDS("out/cv_boot.rds"); bb <- CB$boots[[ING_BEST]]
fig("FigS5_bootstrap_parameters", 7.2, 4.8, dir = DIR, function() {
  par(mfrow = c(2, 3), mar = c(3.6, 3.4, 2, 1))
  items <- list(list("C_per_P", "Target C:P ratio", identity), list("omega", expression("Anisotropy " * omega), identity),
                list("phi_deg", expression("Principal axis " * phi * " (deg)"), wrap180),
                list("excess_deficit_P", "Protein excess:deficit (log10)", log10),
                list("excess_deficit_C", "Carbohydrate excess:deficit (log10)", log10),
                list("lambda2", expression("Volume cost " * lambda[2]), identity))
  for (k in seq_along(items)) { it <- items[[k]]; v <- it[[3]](bb[, it[[1]]]); e <- it[[3]](best$desc[[it[[1]]]])
    hist(v, breaks = 25, col = "#dbeafe", border = "white", main = "", xlab = it[[2]], ylab = if (k %in% c(1, 4)) "Bootstrap replicates" else "")
    abline(v = e, col = ACC2, lwd = 2); abline(v = quantile(v, c(.025, .975)), col = INK, lty = 2)
    if (it[[1]] %in% c("excess_deficit_P", "excess_deficit_C")) abline(v = 0, col = MUTED, lty = 3)
    panel_lab(letters[k]) }
})

fig("FigS1_effective_conics", 7.2, 3.8, dir = DIR, function() {
  par(mfrow = c(1, 2), mar = c(3.6, 3.8, 2, 1)); Fs <- c(45, 90, 180, 360)
  blank(c(0, 145), c(0, 215)); rails_draw(LEE_RATIOS, 400); traj <- NULL
  sy <- fits[["Aniso+L2"]]
  for (k in 1:4) { ec <- effective_conic(sy$p, Fs[k]); lines(clip_pos(ellipse_pts(ec$Qe, ec$centre)), col = SEQ[k], lwd = 1.6)
    kt <- which(dat$F == Fs[k]); points(exp(dat$ybar[kt]) * dat$comp[kt, ], pch = 19, col = SEQ[k], cex = .8)
    traj <- rbind(traj, ec$mu_eff) }
  points(traj, pch = 23, bg = SEQ, col = INK, cex = 1.2); points(sy$p$mu[1], sy$p$mu[2], pch = 4, lwd = 2, col = INK)
  legend("topright", bty = "n", cex = .65, lwd = 1.6, col = SEQ, legend = paste(Fs, "g/L"), text.col = INK)
  panel_lab("a"); panel_title("Exact conic at each concentration")
  invF <- 1 / Fs; ang <- atan2(traj[, 2], traj[, 1]) * 180 / pi; mag <- sqrt(rowSums(traj^2))
  plot(invF, mag, type = "b", pch = 23, bg = SEQ, col = INK, xlab = "1 / food concentration (L/g)",
       ylab = expression("|" * mu[eff] * "| (" * mu * "g day"^-1 * ")"), bty = "l", ylim = c(0, sqrt(sum(sy$p$mu^2)) * 1.15))
  abline(h = sqrt(sum(sy$p$mu^2)), col = ACC2, lty = 2)
  text(mean(range(invF)), sqrt(sum(sy$p$mu^2)), expression("latent target |" * mu * "|"), col = ACC2, cex = .65, pos = 3, offset = .3)
  panel_lab("b"); panel_title("Apparent target shrinks with dilution")
})

J1 <- readRDS("out/fits_joint.rds"); J2 <- readRDS("out/fits_joint_sep.rds")
ch <- load_choice(6)
fig("FigS6_joint_choice", 7.2, 3.6, dir = DIR, function() {
  par(mfrow = c(1, 2), mar = c(3.6, 3.8, 2, 1))
  for (k in 1:2) { pr <- if (k == 1) J1$jf[[ING_BEST]]$pred_ch else J2$js[[ING_BEST]]$pred_ch
    blank(c(0, 80), c(0, 95)); segments(0, 0, 23.75, 95, col = REF, lty = 2)
    segments(ch$P_ug_day, ch$C_ug_day, pr[, 1], pr[, 2], col = GRID)
    points(pr[, 1], pr[, 2], pch = 21, bg = "white", col = ACC, lwd = 1.4)
    points(ch$P_ug_day, ch$C_ug_day, pch = 19, col = SEQ[c(4, 3, 2)][match(ch$sucrose_gL, c(180, 90, 45))])
    panel_lab(letters[k]); panel_title(if (k == 1) "Shared volume cost (gut capacity)" else "Food-specific volume cost")
    if (k == 1) legend("bottomright", bty = "n", cex = .65, pch = c(19, 19, 19, 21), col = c(SEQ[4:2], ACC), pt.bg = c(NA, NA, NA, "white"),
      legend = c("observed, sucrose 180", "observed, sucrose 90", "observed, sucrose 45", "joint ING fit"), text.col = INK) }
})

dm <- read.csv("out/table_trait_domain.csv"); doms <- readRDS("out/trait_domains.rds")
lab <- c(lifespan = "Lifespan", lifetimeegg = "Lifetime egg production", dailyeggs = "Egg production rate")
L <- read_lee("data/lee2008_nochoice.csv"); L <- L[complete.cases(L[, c("protein_eaten", "carb_eaten", names(lab))]), ]
fig("FigS7_thinplate_surfaces", 7.2, 2.9, dir = DIR, function() {
  par(mfrow = c(1, 3), mar = c(3.6, 3.6, 2, .8))
  gx <- seq(0, 185, length.out = 150); gy <- seq(0, 305, length.out = 150); G <- expand.grid(P = gx, C = gy)
  hall <- doms[["individuals (all)"]]
  inh <- function(px, py, h) { n <- nrow(h); ins <- rep(FALSE, length(px)); j <- n
    for (i in 1:n) { cr <- ((h[i, 2] > py) != (h[j, 2] > py)) & (px < (h[j, 1] - h[i, 1]) * (py - h[i, 2]) / (h[j, 2] - h[i, 2]) + h[i, 1]); ins <- xor(ins, cr); j <- i }; ins }
  inside <- inh(G$P, G$C, hall); pch_dom <- c(`treatment means` = 22, `individuals (95% peeled)` = 23, `individuals (all)` = 24)
  for (k in seq_along(lab)) { y <- names(lab)[k]; g <- readRDS(paste0("out/gam_", y, ".rds"))
    z <- rep(NA, nrow(G)); z[inside] <- predict(g, newdata = G[inside, ]); z <- matrix(z, length(gx))
    image(gx, gy, z, col = SURF, xlab = XLAB_P, ylab = YLAB_C, useRaster = TRUE, bty = "l")
    contour(gx, gy, z, add = TRUE, nlevels = 6, col = "white", lwd = .6, labcex = .5)
    for (nm in names(doms)) polygon(doms[[nm]], border = MUTED, lty = match(nm, names(doms)), lwd = .7)
    s <- dm[dm$trait == y & dm$surface == "thinplate", ]
    points(s$P, s$C, pch = pch_dom[s$domain], bg = ACC2, col = "white", cex = 1.4)
    panel_lab(letters[k]); panel_title(lab[k])
    if (k == 1) legend("bottomright", bty = "n", cex = .68, pch = pch_dom, pt.bg = ACC2, col = "white", pt.cex = 1.4, lty = 1:3,
                       legend = paste("optimum,", names(pch_dom)), text.col = INK) }
})

AT <- readRDS("out/adaptive_test.rds"); EH <- readRDS("out/eggrate_hessian.rds")
fig("FigS9_adaptive_profiles", 5, 3.8, dir = DIR, function() {
  par(mar = c(3.6, 3.8, 1, 1))
  plot(NA, xlim = c(-2, 6), ylim = c(0, 30), bty = "l", xlab = expression(log[10] * lambda[2] * "  (weight of volume cost relative to trait)"),
       ylab = expression(Delta * " negative log-likelihood"))
  cl <- c(lifespan = MUTED, lifetimeegg = INK, dailyeggs = ACC); lt <- c(lifespan = 3, lifetimeegg = 2, dailyeggs = 1)
  for (y in names(AT)) lines(AT[[y]]$profile$log10_l2, AT[[y]]$profile$dnll, col = cl[y], lty = lt[y], lwd = 1.8)
  ep <- EH$eta_profile; ll <- seq(-10, 2, by = .5) / log(10); lines(ll, ep - min(ep), col = ACC2, lwd = 1.8)
  legend("topright", bty = "n", cex = .7, lwd = 1.8, col = c(cl, ACC2), lty = c(lt, 1), text.col = INK,
         legend = c("lifespan", "lifetime eggs", "egg rate", "log egg rate (interior peak)"))
  mtext("left: trait dominates feeding; right: constant-volume feeding", 1, 2.6, cex = .6, col = MUTED)
})

fig("FigS10_curvature_matching", 7.2, 4.0, dir = DIR, function() {
  par(mfrow = c(1, 2), mar = c(3.6, 3.8, 2, 1), oma = c(1.4, 0, 0, 0))
  to_Qz <- function(Qraw, mu) { D <- diag(mu); M <- D %*% Qraw %*% D; M / sum(diag(M)) }
  unit_ell <- function(M) { e <- eigen(M, symmetric = TRUE); th <- seq(0, 2 * pi, length.out = 300)
    t(e$vectors %*% rbind(cos(th) / sqrt(e$values[1]), sin(th) / sqrt(e$values[2]))) }
  plot(NA, xlim = c(-2.6, 2.6), ylim = c(-2.4, 3.4), asp = 1, bty = "l",
       xlab = "Proportional protein deviation", ylab = "Proportional carbohydrate deviation")
  abline(h = 0, v = 0, col = GRID)
  lines(unit_ell(best$p$Qz), col = ACC, lwd = 2); lines(unit_ell(symm$p$Qz), col = ACC, lwd = 1.4, lty = 2)
  lines(unit_ell(to_Qz(-EH$H, best$p$mu)), col = ACC2, lwd = 2)
  panel_lab("a"); panel_title("Equal-cost contours around the target")
  v <- EH$boot[, paste0(ING_BEST, ":axis_diff")]
  hist(v, breaks = seq(0, 90, by = 5), col = "#dbeafe", border = "white", main = "", xlab = "Angle between principal axes (deg)",
       ylab = "Bootstrap replicates"); abline(v = EH$point["axis_diff", ING_BEST], col = ACC2, lwd = 2)
  panel_lab("b"); panel_title("Joint bootstrap (replicates with an interior peak)")
  legend_below(ncol = 4, lwd = c(2, 1.4, 2, 2), lty = c(1, 2, 1, 1), col = c(ACC, ACC, ACC2, ACC2),
               legend = c("(a) ING metric, best model", "(a) ING metric, symmetric", "(a) -Hessian of log egg rate", "(b) point estimate"))
})
cat("supplementary figures written to figures/supplementary\n")
