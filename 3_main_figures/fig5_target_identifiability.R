# Fig. 3: profile likelihood of the target ratio
source("R/ing_model.R"); source("R/fig_theme.R")
suppressMessages(library(nutrigon))
ch <- load_choice(6); CP180 <- mean(ch$C_per_P[ch$sucrose_gL == 180])
TN <- readRDS("out/traits_nutrigon.rds")

prof <- readRDS("out/profile_ratio.rds")
ratio_of <- function(y) { r <- TN$tab_opt[TN$tab_opt$trait == y & TN$tab_opt$model == "lm", ]; r$C_opt / r$P_opt }
fig("Fig5_target_identifiability", 7.2, 4.0, function() {
  par(mfrow = c(1, 2), mar = c(3.6, 3.8, 2, 1), oma = c(1.6, 0, 0, 0))
  for (jt in c(FALSE, TRUE)) {
    plot(NA, xlim = c(1, 6), ylim = c(0, 40), log = "", xlab = "Target C:P ratio", ylab = expression(Delta * " negative log-likelihood"), bty = "l")
    panel_lab(if (jt) "b" else "a")
    abline(h = qchisq(.95, 1) / 2, col = GRID, lty = 2)
    tv <- c(lifespan = ratio_of("lifespan"), lifetimeegg = ratio_of("lifetimeegg"), dailyeggs = ratio_of("dailyeggs"))
    abline(v = CP180, col = ACC2, lwd = 1.5); text(CP180, 39, "choice (sucrose 180)", col = ACC2, cex = .7, pos = 2)
    abline(v = tv[2:3], col = "#9ca3af", lty = 3); text(tv[2], 33, "eggs", col = MUTED, cex = .65, pos = 4); text(tv[3], 36, "egg rate", col = MUTED, cex = .65, pos = 2)
    mods <- ING_REP; lt <- c(3, 2, 1); cl <- c(MUTED, INK, ACC)
    for (m in seq_along(mods)) { s <- prof[prof$model == mods[m] & prof$joint == jt, ]; s <- s[order(s$r), ]
      lines(s$r, s$dnll, lty = lt[m], col = cl[m], lwd = 1.8) }
    mtext(if (jt) "No-choice + choice (shared regulation)" else "No-choice array only", 3, .3, cex = .78)
  }
  legend_below(ncol = 3, lty = c(3, 2, 1), col = c(MUTED, INK, ACC), lwd = 1.8,
               legend = c("isotropic, linear cost", "anisotropic, quadratic cost", "best model (asymmetric, capacity)"))
})
cat("Figure 5 written\n")
