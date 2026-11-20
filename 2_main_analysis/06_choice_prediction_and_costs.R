# Prospective and calibrated prediction of the choice experiment
source("R/ing_model.R"); source("R/ing_joint.R"); library(parallel)
F0 <- readRDS("out/fits_nochoice.rds"); dat <- F0$dat; fits <- F0$fits
J0 <- readRDS("out/fits_joint.rds"); jf <- J0$jf; CH <- J0$CH; ch <- J0$ch
lr_obs <- log(ch$C_ug_day / ch$P_ug_day)
oos <- do.call(rbind, lapply(names(fits)[!grepl("NoAcq", names(fits))], function(n) { p <- fits[[n]]$p
  do.call(rbind, lapply(c(FALSE, TRUE), function(sp) {
    rr <- optimize(function(lr) { pc <- p; pc$mu <- exp(lr) * p$mu; pr <- choice_opt_vec(pc, CH$FY, CH$FS, sp)
      sum((CH$lP - log(pr[, 1]))^2 + (CH$lC - log(pr[, 2]))^2) }, c(log(.1), log(3)))
    pr1 <- choice_opt_vec(p, CH$FY, CH$FS, sp); lr1 <- log(pr1[, 2] / pr1[, 1])
    pc <- p; pc$mu <- exp(rr$minimum) * p$mu; pr <- choice_opt_vec(pc, CH$FY, CH$FS, sp)
    lr <- log(pr[, 2] / pr[, 1])
    data.frame(model = n, cost = if (sp) "separable" else "shared", rho = round(exp(rr$minimum), 3),
               ratio_RMSE = round(sqrt(mean((lr - lr_obs)^2)), 3), ratio_cor = round(cor(lr, lr_obs), 3),
               amount_RMSE = round(sqrt(rr$objective / 18), 3),
               CP_S180 = round(mean(exp(lr[ch$sucrose_gL == 180])), 2),
               p1_ratio_RMSE = round(sqrt(mean((lr1 - lr_obs)^2)), 3), p1_ratio_cor = round(cor(lr1, lr_obs), 3),
               p1_amount_RMSE = round(sqrt(mean(c(CH$lP - log(pr1[, 1]), CH$lC - log(pr1[, 2]))^2)), 3),
               p1_CP_S180 = round(mean(exp(lr1[ch$sucrose_gL == 180])), 2)) })) }))
print(oos[order(oos$ratio_RMSE), ], row.names = FALSE)
write.csv(oos, "out/table_choice_oos.csv", row.names = FALSE)
fit_joint_sep <- function(s, start) { nb <- n_par(s); best <- NULL; set.seed(4)
  for (k in 1:40) { p0 <- start + rnorm(length(start), 0, .3) * (k > 1)
    if (k > 20) p0[1:2] <- log(c(30, 90)) + rnorm(2, 0, .6)
    f <- try(optim(p0, negll_joint, dat = dat, s = s, CH = CH, sep = TRUE, control = list(maxit = 8000, reltol = 1e-11)), silent = TRUE)
    if (inherits(f, "try-error") || f$value >= 1e10) next
    f2 <- try(optim(f$par, negll_joint, dat = dat, s = s, CH = CH, sep = TRUE, method = "BFGS"), silent = TRUE)
    if (!inherits(f2, "try-error") && f2$value <= f$value) f <- f2
    if (is.null(best) || f$value < best$value) best <- f }
  kk <- length(best$par); p <- unpack(best$par[1:nb], s); rho <- exp(best$par[nb + 2])
  pc <- p; pc$mu <- rho * p$mu
  list(spec = s, par = best$par, nll = best$value, k = kk, AICc = 2 * best$value + 2 * kk + 2 * kk * (kk + 1) / (46 - kk - 1),
       p = p, rho = rho, tau = exp(best$par[nb + 1]), sig_c = exp(best$par[kk]), desc = describe(p),
       pred_nc = rail_opt(p, dat$comp, dat$F), pred_ch = choice_opt_vec(pc, CH$FY, CH$FS, TRUE), sep = TRUE) }
js <- mclapply(names(jf), function(n) fit_joint_sep(jf[[n]]$spec, jf[[n]]$par), mc.cores = ING_CORES); names(js) <- names(jf)
tab <- do.call(rbind, lapply(names(js), function(n) { f <- js[[n]]; g <- jf[[n]]
  data.frame(model = n, AICc_shared = round(g$AICc, 2), AICc_separable = round(f$AICc, 2),
    sig_c_shared = round(g$sig_c, 3), sig_c_sep = round(f$sig_c, 3), C_per_P = round(f$desc[["C_per_P"]], 2),
    omega = round(f$desc[["omega"]], 3), phi = round(f$desc[["phi_deg"]]), exP = round(f$desc[["excess_deficit_P"]], 3),
    exC = round(f$desc[["excess_deficit_C"]], 3), l1 = round(f$p$l1, 3), l2 = round(f$p$l2, 4), rho = round(f$rho, 3)) }))
tab$dAICc <- round(pmin(tab$AICc_shared, tab$AICc_separable) - min(c(tab$AICc_shared, tab$AICc_separable)), 2)
print(tab[order(tab$AICc_separable), ], row.names = FALSE)
saveRDS(list(js = js, tab = tab, oos = oos), "out/fits_joint_sep.rds")
write.csv(tab, "out/table_joint_shared_vs_separable.csv", row.names = FALSE)
