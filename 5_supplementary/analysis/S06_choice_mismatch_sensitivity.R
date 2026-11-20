# Choice mismatch: yeast composition and sucrose palatability
source("R/ing_model.R"); source("R/fig_theme.R"); source("R/table_utils.R"); library(parallel)
if (!file.exists("out/fits_nochoice.rds")) stop("Run 2_main_analysis/01_fit_nochoice_models.R first")
F0 <- readRDS("out/fits_nochoice.rds"); dat <- F0$dat; fits <- F0$fits
fb <- fits[[ING_BEST]]; p0 <- fb$p; s0 <- fb$spec

ch <- load_choice(6)
FY <- ch$yeast_gL; FS <- ch$sucrose_gL
VY <- ch$V_yeast_uL_6d / 6; VS <- ch$V_sucrose_uL_6d / 6
YEAST0 <- YEAST
obs_intake <- function(comp) cbind(P = ch$P_ug_day * comp[["P"]] / YEAST0[["P"]],
                                   C = ch$C_ug_day + VY * FY * (comp[["C"]] - YEAST0[["C"]]))
OBS0 <- obs_intake(YEAST0)

score <- function(predfun, obs) {
  lP <- log(obs[, "P"]); lC <- log(obs[, "C"])
  ss <- function(lr) { pr <- predfun(exp(lr)); if (any(pr[, 1:2] <= 1e-8)) return(1e10)
    sum((lP - log(pr[, "P"]))^2 + (lC - log(pr[, "C"]))^2) }
  o <- optimize(ss, c(log(.1), log(3))); pr <- predfun(exp(o$minimum))
  lr <- log(pr[, "C"] / pr[, "P"]); lro <- lC - lP
  list(rho = exp(o$minimum), pred = pr, lr = lr, SS = o$objective,
       ratio_RMSE = sqrt(mean((lr - lro)^2)), ratio_cor = cor(lr, lro),
       CP_S180_pred = mean(exp(lr[FS == 180])), CP_S180_obs = mean(exp(lro[FS == 180])))
}
pred_exact <- function(p, comp = YEAST0) function(rho) { pc <- p; pc$mu <- rho * p$mu
  assign("YEAST", comp, envir = globalenv()); on.exit(assign("YEAST", YEAST0, envir = globalenv()))
  choice_opt_vec(pc, FY, FS) }
aicc_ch <- function(SS, k, n = 18) { nll <- n / 2 * (log(2 * pi * SS / n) + 1); 2 * nll + 2 * k + 2 * k * (k + 1) / (n - k - 1) }

choice_num <- function(p, extra = NULL, comp = YEAST0)
  choice_num_gen(p, FY, FS, extra = if (is.null(extra)) NULL else function(G, i) extra(G), comp = comp)
chk <- max(abs(choice_num(p0)[, 1:2] / choice_opt_vec(p0, FY, FS)[, 1:2] - 1))
cat(sprintf("numerical vs exact choice solver, max relative difference: %.1e\n", chk))
stopifnot(chk < 1e-3)

BASE <- score(pred_exact(p0), OBS0)
cat(sprintf("Baseline: ratio RMSE %.3f, r = %.3f, C:P at S180 predicted %.2f vs observed %.2f\n",
            BASE$ratio_RMSE, BASE$ratio_cor, BASE$CP_S180_pred, BASE$CP_S180_obs))

fPs <- round(seq(.10, .60, by = .025), 3); fCs <- sort(unique(round(c(seq(.10, .40, by = .025), YEAST0[["C"]]), 3)))
grid <- expand.grid(fP = fPs, fC = fCs)
resA <- do.call(rbind, lapply(seq_len(nrow(grid)), function(k) {
  comp <- c(P = grid$fP[k], C = grid$fC[k]); sc <- score(pred_exact(p0, comp), obs_intake(comp))
  data.frame(fP = comp[["P"]], fC = comp[["C"]], rho = sc$rho, ratio_RMSE = sc$ratio_RMSE,
             ratio_cor = sc$ratio_cor, CP_S180_pred = sc$CP_S180_pred, CP_S180_obs = sc$CP_S180_obs) }))
resA$log_gap <- log(resA$CP_S180_obs / resA$CP_S180_pred)
write.csv(resA, "out/table_choice_yeast_grid.csv", row.names = FALSE)
A_best <- resA[which.min(resA$ratio_RMSE), ]
A_avail <- resA[abs(resA$fC - YEAST0[["C"]]) < 1e-9, ]
A_gap0 <- A_avail[which.min(abs(A_avail$log_gap)), ]
if (min(abs(A_avail$log_gap)) > .05) cat("(A) note: no protein fraction in the grid closes the S180 gap\n")
cat("(A) best composition:", A_best$fP, A_best$fC, " ratio RMSE", round(A_best$ratio_RMSE, 3), "\n")
cat("(A) protein fraction closing the S180 gap (fC = 0.24):", A_gap0$fP, "\n")

pred_pal <- function(gam) function(rho) { pc <- p0; pc$mu <- rho * p0$mu
  choice_num(pc, extra = function(G) -gam * G[, 2]) }
prof_g <- data.frame(gamma = seq(-0.3, 0.6, by = .025))
prof_g$SS <- sapply(prof_g$gamma, function(g) score(pred_pal(g), OBS0)$SS)
g_hat <- optimize(function(g) score(pred_pal(g), OBS0)$SS,
                  range(prof_g$gamma[order(prof_g$SS)[1:3]]) + c(-.025, .025))$minimum
C_ <- score(pred_pal(g_hat), OBS0)
cat(sprintf("(C) gamma = %.3f; ratio RMSE %.3f; C:P at S180 %.2f\n", g_hat, C_$ratio_RMSE, C_$CP_S180_pred))

row <- function(scen, extra, aicc_nc, sc, aicc_c) data.frame(Scenario = scen, Extra_parameters = extra,
  AICc_nochoice = aicc_nc, rho = sc$rho, Ratio_RMSE = sc$ratio_RMSE, Ratio_r = sc$ratio_cor,
  CP_S180_predicted = sc$CP_S180_pred, CP_S180_observed = sc$CP_S180_obs, AICc_choice = aicc_c)
scA <- score(pred_exact(p0, c(P = A_best$fP, C = A_best$fC)), obs_intake(c(P = A_best$fP, C = A_best$fC)))
scA0 <- score(pred_exact(p0, c(P = A_gap0$fP, C = A_gap0$fC)), obs_intake(c(P = A_gap0$fP, C = A_gap0$fC)))
tab <- rbind(
  row(sprintf("Baseline: %s (yeast 45%% P, 24%% C)", ING_BEST), "--", fb$AICc, BASE, aicc_ch(BASE$SS, 2)),
  row(sprintf("A: yeast %.1f%% P, %.1f%% C (best ratio fit)", 100 * A_best$fP, 100 * A_best$fC), "(composition)", fb$AICc, scA, NA),
  row(sprintf(if (abs(A_gap0$log_gap) > .05) "A: available protein %.1f%% (largest tested; S180 gap not closed)" else "A: available protein %.1f%% (closes S180 gap)", 100 * A_gap0$fP), "(composition)", fb$AICc, scA0, NA),
  row(sprintf("C: sucrose palatability (gamma = %.3f)", g_hat), "gamma (from choice)", fb$AICc, C_, aicc_ch(C_$SS, 3)))
print(tab, row.names = FALSE, digits = 3)
write_table(tab, "TableS3_choice_mismatch",
  paste("Sensitivity of the choice-experiment prediction to yeast composition (A)",
        "and to a sucrose palatability term estimated from the choice data (C). AICc\\_nochoice: fit to the 28",
        "no-choice treatments (best model, ING\\_BEST). AICc\\_choice: Gaussian AICc of the choice data (18 log intakes) given the",
        "no-choice parameters, counting rho, sigma and gamma; not comparable across compositions (A),",
        "because the observed intakes change with composition. C:P at S180: mean over the three yeast",
        "concentrations at sucrose 180 g/L."),
  "tab:S3", "tables/supplementary", digits = 3)
saveRDS(list(base = BASE, gridA = resA, A_best = A_best, A_gap0 = A_gap0,
             gamma = g_hat, prof_gamma = prof_g, C = C_, table = tab), "out/choice_mismatch.rds")

lro <- log(OBS0[, "C"] / OBS0[, "P"])
fig("FigS2_choice_mismatch", 7.2, 6.6, dir = "figures/supplementary", function() {
  layout(matrix(c(1, 2, 3, 4), 2, byrow = TRUE)); par(mar = c(3.8, 3.8, 2.2, 1))
  Z <- matrix(exp(resA$log_gap), length(fPs), length(fCs))
  image(fPs * 100, fCs * 100, Z, col = colorRampPalette(c("white", "#fdba74", "#c2410c"))(40),
        zlim = c(1, max(Z, 1.01)), xlab = "Yeast protein (% mass)", ylab = "Yeast carbohydrate (% mass)")
  contour(fPs * 100, fCs * 100, Z, levels = c(1, 1.5, 2, 3, 4, 6, 8), add = TRUE, labcex = .65, col = INK, lwd = .8)
  points(45, 24, pch = 4, lwd = 2, col = INK); text(45, 24, "assumed", pos = 4, cex = .7, col = INK)
  panel_lab("a"); panel_title("A: observed / predicted C:P at sucrose 180")
  plot(A_avail$fP * 100, A_avail$CP_S180_obs, type = "l", lwd = 2, col = INK, bty = "l",
       ylim = range(c(A_avail$CP_S180_obs, A_avail$CP_S180_pred, 0)), xlab = "Available yeast protein (% mass)",
       ylab = "C:P at sucrose 180")
  lines(A_avail$fP * 100, A_avail$CP_S180_pred, lwd = 2, col = ACC); abline(v = 45, col = GRID, lty = 3)
  legend("topright", bty = "n", cex = .7, lwd = 2, col = c(INK, ACC), legend = c("observed", "predicted"))
  panel_lab("b"); panel_title("A: protein availability (yeast C = 24%)")
  cols <- c(MUTED, ACC2); sc <- list(BASE, C_)
  lim <- range(c(lro, unlist(lapply(sc, `[[`, "lr"))))
  plot(NA, xlim = lim, ylim = lim, bty = "l", xlab = "Predicted log(C:P)", ylab = "Observed log(C:P)")
  abline(0, 1, col = GRID, lty = 2)
  for (k in 1:2) points(sc[[k]]$lr, lro, pch = c(1, 16)[k], col = cols[k], cex = 1.1)
  legend("bottomright", bty = "n", cex = .7, pch = c(1, 16), col = cols,
         legend = c(sprintf("best model (RMSE %.2f)", BASE$ratio_RMSE), sprintf("+ palatability (RMSE %.2f)", C_$ratio_RMSE)))
  panel_lab("c"); panel_title("Choice ratios: observed vs predicted")
  plot(prof_g$gamma, prof_g$SS - min(prof_g$SS), type = "l", lwd = 2, col = ACC, bty = "l",
       xlab = expression("Sucrose palatability " * gamma), ylab = expression(Delta * " sum of squares (log intakes)"))
  abline(v = c(0, g_hat), col = c(MUTED, ACC2), lty = c(3, 2))
  panel_lab("d"); panel_title("C: palatability estimated from choice data")
})
cat("S06 done: Table S3, Fig. S2\n")
