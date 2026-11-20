# Fit the 64 ING models to the no-choice array
source("R/ing_model.R")
d   <- read_lee("data/lee2008_nochoice.csv"); dat <- prep_intake(d)
grid <- expand.grid(cost = c("none", "L1", "L2", "L12"), asym = c("", "P", "C", "PC"), aniso = c(FALSE, TRUE),
                    cap = c(FALSE, TRUE), stringsAsFactors = FALSE)
specs <- lapply(seq_len(nrow(grid)), function(i) { g <- grid[i, ]
  spec(g$aniso, asym = switch(g$asym, P = "P", C = "C", PC = c("P", "C"), character(0)),
       l1 = g$cost %in% c("L1", "L12"), l2 = g$cost %in% c("L2", "L12"), cap = g$cap) })
library(parallel); fits <- list()
nocap <- Filter(function(s) !s$cap, specs); withcap <- Filter(function(s) s$cap, specs)
fl <- mclapply(nocap, function(s) fit_ing(dat, s, n_start = 80, seed = 11), mc.cores = ING_CORES); names(fl) <- sapply(nocap, `[[`, "name")
fc <- mclapply(withcap, function(s) {
  s0 <- spec(s$aniso, s$asym, s$l1, s$l2); f0 <- fl[[s0$name]]; nb <- n_par(s0); best <- fit_ing(dat, s, n_start = 80, seed = 11)
  for (vm in log(c(1.5, 1.9, 2.5, 4))) { st <- c(f0$par[1:nb], vm, f0$par[nb + 1])
    f <- fit_ing(dat, s, n_start = 20, seed = 12, start = st); if (f$nll < best$nll) best <- f }
  best }, mc.cores = ING_CORES); names(fc) <- sapply(withcap, `[[`, "name")
fl <- c(fl, fc)
for (s in specs) { fits[[s$name]] <- fl[[s$name]]
  f <- fits[[s$name]]; cat(sprintf("%-24s k=%2d nll=%8.3f AICc=%8.3f tau=%.3f  C:P=%.2f omega=%.3f phi=%.0f exP=%.2f exC=%.2f l1=%.3f l2=%.3f Vmax=%.2f\n",
    s$name, f$k, f$nll, f$AICc, f$tau, f$desc["C_per_P"], f$desc["omega"], f$desc["phi_deg"],
    f$desc["excess_deficit_P"], f$desc["excess_deficit_C"], f$desc["lambda1"], f$desc["lambda2"], vmax_of(f$p))) }
best <- names(fits)[which.min(sapply(fits, `[[`, "AICc"))]
cat("Lowest AICc:", best, "\n")
if (best != ING_BEST) warning("ING_BEST in R/ing_model.R (", ING_BEST, ") is not the lowest-AICc model (", best, ")")
saveRDS(list(dat = dat, fits = fits), "out/fits_nochoice.rds")
