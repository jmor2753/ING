# Profile likelihood of the target ratio
source("R/ing_model.R"); source("R/ing_joint.R"); library(parallel)
F0 <- readRDS("out/fits_nochoice.rds"); dat <- F0$dat; fits <- F0$fits
J0 <- readRDS("out/fits_joint.rds"); jf <- J0$jf; CH <- J0$CH
ratios <- c(1, 1.25, 1.5, 1.75, 2, 2.5, 3, 3.5, 4, 5, 6)
prof_one <- function(s, r, joint, start) {
  nb <- n_par(s)
  build <- function(q) { lm <- q[1]; c(lm - log(1 + r) + log(1), lm - log(1 + r) + log(r), q[-1]) }
  f_obj <- function(q) { par <- build(q); if (joint) negll_joint(par, dat, s, CH) else negll(par, dat, s) }
  q0 <- c(log(sum(exp(start[1:2]))), start[-(1:2)]); best <- NULL; set.seed(1)
  for (k in 1:25) { qq <- q0 + c(rnorm(1, 0, .3), rnorm(length(q0) - 1, 0, .5)) * (k > 1)
    f <- try(optim(qq, f_obj, control = list(maxit = 6000, reltol = 1e-11)), silent = TRUE)
    if (inherits(f, "try-error")) next
    if (is.null(best) || f$value < best$value) best <- f }
  best$value
}
jobs <- expand.grid(model = ING_REP, joint = c(FALSE, TRUE), r = ratios,
                    stringsAsFactors = FALSE)
spec_of <- function(n) fits[[n]]$spec
vals <- unlist(mclapply(seq_len(nrow(jobs)), function(i) {
  n <- jobs$model[i]; s <- spec_of(n)
  st <- if (jobs$joint[i]) jf[[n]]$par else fits[[n]]$par
  prof_one(s, jobs$r[i], jobs$joint[i], st) }, mc.cores = ING_CORES))
jobs$nll <- vals
jobs$dnll <- ave(jobs$nll, jobs$model, jobs$joint, FUN = function(v) v - min(v))
print(reshape(jobs[, c("model", "joint", "r", "dnll")], idvar = c("model", "joint"), timevar = "r", direction = "wide"), digits = 3)
saveRDS(jobs, "out/profile_ratio.rds")

negll_joint_free <- function(par, dat, s, CH) {
  nb <- n_par(s); a <- negll(par[1:(nb + 1)], dat, s); if (a >= 1e10) return(1e10)
  p <- unpack(par[1:nb], s); pc <- p; pc$mu <- exp(par[nb + 2:3])
  pr <- choice_opt_vec(pc, CH$FY, CH$FS); if (any(pr[, 1:2] <= 1e-6)) return(1e10)
  sig <- exp(par[nb + 4]); r <- c(CH$lP - log(pr[, "P"]), CH$lC - log(pr[, "C"]))
  v <- sig^2 + c(CH$se_lP, CH$se_lC)^2
  a + 0.5 * sum(log(2 * pi * v) + r^2 / v)
}
free <- mclapply(ING_REP, function(n) {
  s <- spec_of(n); st <- jf[[n]]$par; nb <- n_par(s); best <- NULL; set.seed(2)
  for (k in 1:30) {
    p0 <- c(st[1:(nb + 1)] + rnorm(nb + 1, 0, .2) * (k > 1), log(c(20, 75)) + rnorm(2, 0, .3), log(.2))
    f <- try(optim(p0, negll_joint_free, dat = dat, s = s, CH = CH, control = list(maxit = 8000, reltol = 1e-11)), silent = TRUE)
    if (inherits(f, "try-error")) next
    if (is.null(best) || f$value < best$value) best <- f }
  kk <- length(best$par); N <- 28 + 18
  list(model = n, nll = best$value, k = kk, AICc = 2 * best$value + 2 * kk + 2 * kk * (kk + 1) / (N - kk - 1),
       mu_nc = exp(best$par[1:2]), mu_ch = exp(best$par[nb + 2:3]), sig = exp(best$par[kk]), par = best$par)
}, mc.cores = ING_CORES)
for (f in free) cat(sprintf("%-16s free choice target: AICc=%.2f (shared-geometry AICc=%.2f)  no-choice C:P=%.2f  choice C:P=%.2f  sig_c=%.3f\n",
  f$model, f$AICc, jf[[f$model]]$AICc, f$mu_nc[2] / f$mu_nc[1], f$mu_ch[2] / f$mu_ch[1], f$sig))
saveRDS(free, "out/fits_joint_freetarget.rds")
