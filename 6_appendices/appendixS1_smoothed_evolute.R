# Appendix S1: evolutes from smoothed arrays
source("R/ing_model.R"); source("R/fig_theme.R"); source("R/table_utils.R")
F0  <- readRDS("out/fits_nochoice.rds"); dat <- F0$dat
fit <- F0$fits[["Aniso+L2"]]; p <- fit$p
NREP   <- 500
ratios <- c(0, 1/16, 1/8, 1/4, 1/2, 1, 1.875)
th_obs <- atan2(1, ratios)
th_ev  <- seq(min(th_obs), max(th_obs), length.out = 61)
SE     <- median(sqrt(dat$s2n))
NOISE  <- c(none = 0, sampling = SE, total = sqrt(fit$tau^2 + SE^2))
DFS    <- 3:6

r_true <- function(th, F) { comp <- cbind(cos(th), sin(th)); comp <- comp / rowSums(comp)
  rail_opt(p, comp, rep(F, length(th))) * sqrt(rowSums(comp^2)) }

centre <- function(th, r, r1, r2) {
  k  <- (r^2 + 2 * r1^2 - r * r2) / (r^2 + r1^2)^1.5
  dx <- r1 * cos(th) - r * sin(th); dy <- r1 * sin(th) + r * cos(th)
  nm <- sqrt(dx^2 + dy^2)
  cbind(x = r * cos(th) - dy / nm / k, y = r * sin(th) + dx / nm / k, R = 1 / k) }

res <- list()
for (F in c(45, 90, 180, 360)) {
  h <- 1e-4; r0 <- r_true(th_ev, F)
  r1 <- (r_true(th_ev + h, F) - r_true(th_ev - h, F)) / (2 * h)
  r2 <- (r_true(th_ev + h, F) - 2 * r0 + r_true(th_ev - h, F)) / h^2
  C_true <- centre(th_ev, r0, r1, r2)
  scale  <- sqrt(sum(effective_conic(p, F)$mu_eff^2))
  r_obs  <- r_true(th_obs, F)
  for (nz in names(NOISE)) for (df in DFS) {
    set.seed(1e4 * F + 100 * match(nz, names(NOISE)) + df)
    m <- t(replicate(NREP, {
      s  <- smooth.spline(th_obs, log(r_obs) + rnorm(7, 0, NOISE[[nz]]), df = df)
      g  <- predict(s, th_ev)$y; g1 <- predict(s, th_ev, deriv = 1)$y; g2 <- predict(s, th_ev, deriv = 2)$y
      r  <- exp(g); C <- centre(th_ev, r, r * g1, r * (g2 + g1^2))
      err <- sqrt(rowSums((C[, 1:2] - C_true[, 1:2])^2)) / scale
      c(median(err), mean(sign(C[, "R"]) != sign(C_true[, "R"])), mean(abs(C[, "R"]) > 3 * abs(C_true[, "R"]))) }))
    res[[length(res) + 1]] <- data.frame(F_gL = F, noise = nz, noise_sd_log = round(NOISE[[nz]], 3), spline_df = df,
      median_error_pct_mu = round(100 * mean(m[, 1]), 1), sign_flip_pct = round(100 * mean(m[, 2]), 1),
      radius_3x_pct = round(100 * mean(m[, 3]), 1))
  }
}
tab <- do.call(rbind, res)
dir.create("tables/appendix", recursive = TRUE, showWarnings = FALSE)
write_table(tab, "AppendixS1_Table1_smoothed_evolute",
  paste("Evolutes estimated by smoothing the no-choice array (cubic spline in the angle of each rail) recover",
        "the true evolute poorly. Error: median distance from the true centre of curvature, in percent of",
        "|mu|; sign flips: percentage of points whose radius of curvature has the wrong sign; radius 3x:",
        "percentage of points with radius more than three times the true value. Noise: none, within-treatment",
        "sampling error, or total (sampling plus treatment-level) error of the Lee et al. data."),
  "tab:S1.1", "tables/appendix", digits = 1)

s <- aggregate(cbind(median_error_pct_mu, sign_flip_pct, radius_3x_pct) ~ noise + spline_df, tab, mean)
print(s[order(s$noise, s$spline_df), ], row.names = FALSE)

F <- 180; h <- 1e-4; r0 <- r_true(th_ev, F)
C_true <- centre(th_ev, r0, (r_true(th_ev + h, F) - r_true(th_ev - h, F)) / (2 * h),
                 (r_true(th_ev + h, F) - 2 * r0 + r_true(th_ev - h, F)) / h^2)
set.seed(7); y <- log(r_true(th_obs, F)) + rnorm(7, 0, NOISE[["total"]])
fig("AppendixS1_Fig1_smoothed_evolute", 7.2, 4.2, dir = "figures/appendix", function() {
  par(mfrow = c(1, 2), mar = c(3.6, 3.8, 2, 1), oma = c(2.2, 0, 0, 0))
  blank(c(-60, 160), c(-40, 180)); rails_draw(LEE_RATIOS, 400)
  lines(r0 * cos(th_ev), r0 * sin(th_ev), col = INK, lwd = 2); lines(C_true[, 1:2], col = INK, lwd = 1.2, lty = 2)
  cols <- c(`4` = ACC, `6` = ACC2)
  for (df in c(4, 6)) { sp <- smooth.spline(th_obs, y, df = df)
    g <- predict(sp, th_ev)$y; g1 <- predict(sp, th_ev, deriv = 1)$y; g2 <- predict(sp, th_ev, deriv = 2)$y
    r <- exp(g); C <- centre(th_ev, r, r * g1, r * (g2 + g1^2))
    lines(r * cos(th_ev), r * sin(th_ev), col = cols[as.character(df)], lwd = 1.4)
    C[abs(C[, "R"]) > 5 * max(abs(C_true[, "R"])), 1:2] <- NA
    lines(C[, 1:2], col = cols[as.character(df)], lty = 2, lwd = 1) }
  points(exp(y) * cos(th_obs), exp(y) * sin(th_obs), pch = 19, col = INK, cex = .8)
  panel_lab("a"); panel_title("One noisy realisation, 180 g/L")
  s2 <- aggregate(median_error_pct_mu ~ noise + spline_df, tab, mean)
  par(mar = c(3.6, 4.4, 2, 1))
  plot(NA, xlim = c(3, 6), ylim = c(0, 100), xaxt = "n", bty = "l", xlab = "Spline degrees of freedom",
       ylab = "Curvature-centre error (% of |target|)"); axis(1, 3:6)
  nz <- c(none = MUTED, sampling = ACC, total = ACC2)
  for (n in names(nz)) { z <- s2[s2$noise == n, ]; lines(z$spline_df, z$median_error_pct_mu, type = "b", pch = 19, col = nz[n], lwd = 1.6) }
  legend("topleft", bty = "n", cex = .65, col = nz, lwd = 1.6, pch = 19, text.col = INK,
         legend = c("no noise (smoothing bias)", "sampling error", "Lee-level deviation"))
  panel_lab("b"); panel_title("Mean over concentrations")
  legend_below(ncol = 2, lty = c(1, 2, 1, 1), lwd = c(2, 1.2, 1.4, 1.4), col = c(INK, INK, ACC, ACC2),
               legend = c("(a) true array", "(a) true evolute", "(a) smoothed, df = 4 (dashed: its evolute)", "(a) smoothed, df = 6 (dashed: its evolute)"))
})
cat("Appendix S1 table and figure written\n")
