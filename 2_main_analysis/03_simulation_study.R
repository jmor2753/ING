# Simulation study on the Lee design (7 scenarios x 100 data sets)
source("R/ing_model.R"); library(parallel)
F0 <- readRDS("out/fits_nochoice.rds"); dat <- F0$dat; fits <- F0$fits
args <- commandArgs(TRUE); R_REP <- if (length(args)) as.integer(args[1]) else 100
cands <- list(spec(FALSE), spec(TRUE), spec(FALSE, asym = c("P", "C")), spec(TRUE, asym = c("P", "C")),
              spec(FALSE, l1 = TRUE, l2 = FALSE), fits[[ING_BEST]]$spec)
names(cands) <- sapply(cands, `[[`, "name")
gen <- list(S1 = "Iso+L2", S2 = "Aniso+L2", S3 = "Iso+AsymP+L2", S4 = ING_BEST,
            S5 = "Aniso+L2", S6 = "Aniso+L2", S7 = "Aniso+L2")
rails4 <- which(dat$tr$Ratio %in% c("(1:16)", "(1:4)", "(1:1)", "(2:1)"))
one_rep <- function(sc, r) {
  set.seed(10000 * match(sc, names(gen)) + r)
  g <- fits[[gen[[sc]]]]; dd <- dat
  if (sc == "S6") dd <- subset_trt(dat, rails4)
  db <- switch(sc,
    S5 = simulate_dat(dd, g$p, g$tau, het_cv = 0.25),
    S7 = simulate_dat(dd, g$p, 2 * g$tau, s_i = 2 * dd$s_i),
    simulate_dat(dd, g$p, g$tau))
  fs <- lapply(cands, function(s) {
    st <- fits[[s$name]]$par
    fit_ing(db, s, n_start = 8, seed = r, start = st) })
  a <- sapply(fs, `[[`, "AICc")
  gs <- if (gen[[sc]] %in% names(fs)) fs[[gen[[sc]]]] else
        fit_ing(db, g$spec, n_start = 8, seed = r, start = g$par)
  c(sel = unname(which.min(a)), gs$desc, sapply(fs, function(f) f$desc[["C_per_P"]]))
}
out <- list()
for (sc in names(gen)) {
  rows <- mclapply(seq_len(R_REP), function(r) tryCatch(one_rep(sc, r), error = function(e) NULL), mc.cores = ING_CORES)
  out[[sc]] <- do.call(rbind, rows); cat(sc, "done", nrow(out[[sc]]), "\n")
  saveRDS(list(out = out, gen = gen, cands = names(cands)), "out/simulations.rds")
}
truth <- lapply(gen, function(n) fits[[n]]$desc)
summ <- do.call(rbind, lapply(names(out), function(sc) { o <- out[[sc]]; tr <- truth[[sc]]
  sel <- table(factor(names(cands)[o[, "sel"]], levels = names(cands))) / nrow(o)
  q <- function(v) sprintf("%.3g [%.3g, %.3g]", median(v), quantile(v, .025), quantile(v, .975))
  data.frame(scenario = sc, generating = gen[[sc]], n = nrow(o),
    true_CP = round(tr[["C_per_P"]], 2), est_CP = q(o[, "C_per_P"]),
    true_omega = round(tr[["omega"]], 3), est_omega = q(o[, "omega"]),
    true_phi = round(tr[["phi_deg"]]), est_phi = q(o[, "phi_deg"]),
    true_exP = round(tr[["excess_deficit_P"]], 2), est_exP = q(o[, "excess_deficit_P"]),
    true_l2 = signif(tr[["lambda2"]], 3), est_l2 = q(o[, "lambda2"]),
    setNames(as.list(round(100 * as.numeric(sel))), paste0("pct_sel_", names(cands))),
    check.names = FALSE) }))
print(summ, row.names = FALSE)
write.csv(summ, "out/table_simulations.csv", row.names = FALSE)
