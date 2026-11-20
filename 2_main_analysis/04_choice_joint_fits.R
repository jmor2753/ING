# Joint fits of the no-choice and choice data
source("R/ing_model.R"); library(parallel)
F0 <- readRDS("out/fits_nochoice.rds"); dat <- F0$dat; fits <- F0$fits
ch <- load_choice(6)
CH <- choice_CH(ch)

source("R/ing_joint.R")

specs <- list(spec(FALSE, l1 = TRUE, l2 = FALSE), spec(FALSE), spec(TRUE),
              spec(TRUE, l1 = TRUE, l2 = FALSE),
              spec(FALSE, asym = "P"), spec(FALSE, asym = "C"), spec(FALSE, asym = c("P", "C")),
              spec(TRUE, asym = "P"), spec(TRUE, asym = "C"), spec(TRUE, asym = c("P", "C")),
              fits[[ING_BEST]]$spec,
              fits[["Aniso+AsymP+L1+Cap"]]$spec, fits[["Aniso+AsymPC+L1+Cap"]]$spec,
              fits[["Iso+AsymP+L1+Cap"]]$spec, fits[["Iso+AsymP+L1"]]$spec)
jf <- mclapply(specs, function(s) {
  st <- if (!is.null(fits[[s$name]])) fits[[s$name]]$par else NULL
  fit_joint(dat, s, CH, start_nc = st, n_start = 40, seed = 7)
}, mc.cores = ING_CORES)
names(jf) <- sapply(specs, `[[`, "name")
tab <- do.call(rbind, lapply(names(jf), function(n) { f <- jf[[n]]
  nb <- n_par(f$spec); nc <- negll(f$par[1:(nb + 1)], dat, f$spec)
  data.frame(model = n, k = f$k, nll = round(f$nll, 2), nll_nochoice = round(nc, 2),
             nll_choice = round(f$nll - nc, 2), AICc = round(f$AICc, 2),
             C_per_P = round(f$desc[["C_per_P"]], 2), muP = round(f$p$mu[1], 1), muC = round(f$p$mu[2], 1),
             rho = round(f$rho, 3), omega = round(f$desc[["omega"]], 3), phi = round(f$desc[["phi_deg"]], 0),
             exP = round(f$desc[["excess_deficit_P"]], 3), exC = round(f$desc[["excess_deficit_C"]], 3),
             l1 = round(f$p$l1, 3), l2 = round(f$p$l2, 4), tau = round(f$tau, 3), sig_c = round(f$sig_c, 3)) }))
tab$dAICc <- round(tab$AICc - min(tab$AICc), 2)
print(tab[order(tab$AICc), ], row.names = FALSE)
saveRDS(list(jf = jf, tab = tab, CH = CH, ch = ch), "out/fits_joint.rds")
write.csv(tab[order(tab$AICc), ], "out/table_joint_models.csv", row.names = FALSE)
