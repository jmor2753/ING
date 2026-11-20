# Appendix S3: Jensen et al. (2015) no-choice array
source("R/ing_model.R"); source("R/table_utils.R"); library(parallel)
W_FOOD <- c(1, 1)
J <- read.csv("data/jensen2015_nochoice_means.csv"); J <- J[J$F_gL > 0, ]

prep_means <- function(d) {
  t <- d$P_ug_day + d$C_ug_day
  se_t <- ifelse(d$P_ug_day > 0, d$P_se * t / d$P_ug_day, d$C_se * t / d$C_ug_day)
  r <- d$P_gL / d$C_gL
  list(comp = cbind(d$P_gL, d$C_gL) / d$F_gL, F = d$F_gL, ybar = log(t), s2n = (se_t / t)^2,
       tr = data.frame(diet = d$diet, ratio = d$ratio, F = d$F_gL), t = t)
}
grid <- expand.grid(cost = c("none", "L1", "L2", "L12"), asym = c("", "P", "C", "PC"), aniso = c(FALSE, TRUE),
                    cap = c(FALSE, TRUE), stringsAsFactors = FALSE)
specs <- lapply(seq_len(nrow(grid)), function(i) { g <- grid[i, ]
  spec(g$aniso, asym = switch(g$asym, P = "P", C = "C", PC = c("P", "C"), character(0)),
       l1 = g$cost %in% c("L1", "L12"), l2 = g$cost %in% c("L2", "L12"), cap = g$cap) })
fit_all <- function(dat) {
  nocap <- Filter(function(s) !s$cap, specs); withcap <- Filter(function(s) s$cap, specs)
  fl <- mclapply(nocap, function(s) fit_ing(dat, s, n_start = 40, seed = 21), mc.cores = ING_CORES)
  names(fl) <- sapply(nocap, `[[`, "name")
  fc <- mclapply(withcap, function(s) {
    s0 <- spec(s$aniso, s$asym, s$l1, s$l2); f0 <- fl[[s0$name]]; nb <- n_par(s0); best <- fit_ing(dat, s, n_start = 40, seed = 21)
    for (vm in log(c(0.6, 0.9, 1.2, 2))) { st <- c(f0$par[1:nb], vm, f0$par[nb + 1])
      f <- fit_ing(dat, s, n_start = 12, seed = 22, start = st); if (f$nll < best$nll) best <- f }
    best }, mc.cores = ING_CORES)
  names(fc) <- sapply(withcap, `[[`, "name"); c(fl, fc)
}
profile_ratio <- function(dat, fit, rs = c(1, 1.25, 1.5, 1.75, 2, 2.5, 3, 3.5, 4, 5, 6), n_start = 20) {
  s <- fit$spec; np <- length(fit$par)
  sapply(rs, function(r) {
    f <- function(q) { m <- exp(q[1]); par <- c(log(m / (1 + r)), log(m * r / (1 + r)), q[-1]); negll(par, dat, s) }
    q0 <- c(log(sum(exp(fit$par[1:2]))), fit$par[3:np]); set.seed(31)
    best <- Inf
    for (k in seq_len(n_start)) { st <- if (k == 1) q0 else q0 + rnorm(length(q0), 0, .5)
      o <- try(optim(st, f, control = list(maxit = 4000)), silent = TRUE)
      if (!inherits(o, "try-error")) best <- min(best, o$value) }
    best - fit$nll }) }

res <- list()
for (sx in c("F", "M")) {
  dat <- prep_means(J[J$sex == sx, ]); fits <- fit_all(dat)
  tab <- data.frame(model = names(fits), k = sapply(fits, `[[`, "k"), AICc = sapply(fits, `[[`, "AICc"),
                    tau = sapply(fits, `[[`, "tau"), CP = sapply(fits, function(f) f$desc[["C_per_P"]]),
                    omega = sapply(fits, function(f) f$desc[["omega"]]), Vmax = sapply(fits, function(f) f$desc[["Vmax"]]))
  tab$dAICc <- tab$AICc - min(tab$AICc); tab <- tab[order(tab$AICc), ]
  best <- fits[[tab$model[1]]]
  prof <- data.frame(r = c(1, 1.25, 1.5, 1.75, 2, 2.5, 3, 3.5, 4, 5, 6))
  prof$best <- profile_ratio(dat, best, prof$r)
  if ("Aniso+AsymPC+L2+Cap" %in% names(fits)) prof$lee_best <- profile_ratio(dat, fits[["Aniso+AsymPC+L2+Cap"]], prof$r)
  res[[sx]] <- list(dat = dat, fits = fits, table = tab, best = best, profile = prof)
  cat("\n==== Jensen et al. 2015,", sx, "====\n"); print(head(tab, 8), digits = 3)
  cat("Best:", tab$model[1], " target C:P =", round(best$desc[["C_per_P"]], 2), "\n"); print(round(prof, 2))
}
saveRDS(res, "out/jensen2015_fits.rds")
cat("S10 done\n")
