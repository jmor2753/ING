# Leave-one-rail-out and leave-one-concentration-out CV; parametric bootstrap
source("R/ing_model.R"); library(parallel)
F0 <- readRDS("out/fits_nochoice.rds"); dat <- F0$dat; fits <- F0$fits
cv_one <- function(n, scheme) {
  f <- fits[[n]]; grp <- if (scheme == "rail") dat$tr$Ratio else dat$F; err <- c()
  for (g in unique(grp)) { tin <- which(grp != g); tout <- which(grp == g)
    fi <- fit_ing(subset_trt(dat, tin), f$spec, n_start = 20, seed = 3, start = f$par)
    pr <- rail_opt(fi$p, dat$comp[tout, , drop = FALSE], dat$F[tout])
    err <- c(err, dat$ybar[tout] - log(pmax(pr, 1e-6))) }
  sqrt(mean(err^2))
}
jobs <- expand.grid(model = names(fits), scheme = c("rail", "conc"), stringsAsFactors = FALSE)
jobs$RMSE <- unlist(mclapply(seq_len(nrow(jobs)), function(i) cv_one(jobs$model[i], jobs$scheme[i]), mc.cores = ING_CORES))
cvt <- reshape(jobs, idvar = "model", timevar = "scheme", direction = "wide")
cvt$AICc <- sapply(cvt$model, function(n) round(fits[[n]]$AICc, 2))
cvt$pct_rail <- round(100 * (exp(cvt$RMSE.rail) - 1), 1); cvt$pct_conc <- round(100 * (exp(cvt$RMSE.conc) - 1), 1)
print(cvt[order(cvt$RMSE.rail), ], row.names = FALSE)
write.csv(cvt, "out/table_cv.csv", row.names = FALSE)

boot_model <- function(n, B = 200) {
  f <- fits[[n]]
  rows <- mclapply(seq_len(B), function(b) { set.seed(1000 + b)
    db <- simulate_dat(dat, f$p, f$tau)
    fb <- fit_ing(db, f$spec, n_start = 12, seed = b, start = f$par)
    c(fb$desc, tau = fb$tau, setNames(fb$par, paste0(".par", seq_along(fb$par)))) }, mc.cores = ING_CORES)
  do.call(rbind, rows)
}
boots <- lapply(setNames(nm = unique(c("Iso+L1", "Aniso+L2", "Aniso+AsymPC+L2", ING_BEST))), boot_model)
ci <- do.call(rbind, lapply(names(boots), function(n) { b <- boots[[n]]; b <- b[, !grepl("^\\.par", colnames(b))]; est <- c(fits[[n]]$desc, tau = fits[[n]]$tau)
  data.frame(model = n, parameter = colnames(b), estimate = signif(est[colnames(b)], 4),
             lo = signif(apply(b, 2, quantile, .025, na.rm = TRUE), 4), hi = signif(apply(b, 2, quantile, .975, na.rm = TRUE), 4)) }))
print(ci, row.names = FALSE)
saveRDS(list(cv = cvt, boots = boots, ci = ci), "out/cv_boot.rds"); write.csv(ci, "out/table_param_ci.csv", row.names = FALSE)
