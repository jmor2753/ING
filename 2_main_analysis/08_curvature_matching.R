# Curvature of regulation vs curvature of log egg rate
source("R/ing_model.R"); suppressMessages(library(nutrigon)); library(parallel)
F0 <- readRDS("out/fits_nochoice.rds"); dat <- F0$dat; fits_ing <- F0$fits
XC <- "protein_eaten"; YC <- "carb_eaten"; data(lee_drosophila)
fit_glm <- function(dd) { set.seed(1); nutrigon_fit(dd, XC, YC, "dailyeggs", model_type = "glm",
  family = quasipoisson(), n_cloud = 1000, n_reps = 3) }
hess_of <- function(f) { cf <- nutrigon_coefs(f)
  matrix(c(2 * cf[["I(.x^2)"]], cf[[".x:.y"]], cf[[".x:.y"]], 2 * cf[["I(.y^2)"]]), 2) }
to_Qz <- function(Qraw, mu) { D <- diag(mu); M <- D %*% Qraw %*% D; M / sum(diag(M)) }
axis_deg <- function(M) { v <- eigen(M, symmetric = TRUE)$vectors[, 1]; a <- atan2(v[2], v[1]) * 180 / pi
  a <- a %% 180; if (a > 90) a - 180 else a }
compare <- function(Qz_ing, negH, mu) { W <- to_Qz(negH, mu)
  a1 <- axis_deg(Qz_ing); a2 <- axis_deg(W); dd <- abs(a1 - a2); dd <- min(dd, 180 - dd)
  c(axis_ING = a1, axis_W = a2, axis_diff = dd,
    omega_ING = eigen(Qz_ing, symmetric = TRUE)$values[1], omega_W = eigen(W, symmetric = TRUE)$values[1]) }

g0 <- fit_glm(lee_drosophila); H0 <- hess_of(g0)
cat("log egg-rate surface: closed-form interior peak =", !g0$optimum$used_fallback,
    " eigenvalues of H:", signif(eigen(H0)$values, 3), "\n")
mods <- c("Aniso+L2", ING_BEST)
pt <- sapply(mods, function(n) compare(fits_ing[[n]]$p$Qz, -H0, fits_ing[[n]]$p$mu))

cf <- nutrigon_coefs(g0); T_MAX <- 2 * max(exp(dat$ybar))
rail_eta <- function(l1, l2) sapply(seq_len(nrow(dat$comp)), function(i) {
  c1 <- dat$comp[i, 1]; c2 <- dat$comp[i, 2]; g <- sum(W_FOOD * dat$comp[i, ]) / dat$F[i]
  a2 <- -(cf[["I(.x^2)"]] * c1^2 + cf[["I(.y^2)"]] * c2^2 + cf[[".x:.y"]] * c1 * c2) + l2 * g^2
  a1 <- -(cf[[".x"]] * c1 + cf[[".y"]] * c2) + l1 * g
  cand <- c(1e-3, T_MAX); if (a2 > 0) cand <- c(cand, min(max(-a1 / (2 * a2), 1e-3), T_MAX))
  cand[which.min(a2 * cand^2 + a1 * cand)] })
nll_eta <- function(q) { ts <- rail_eta(q[1], exp(q[2])); v <- exp(2 * q[3]) + dat$s2n
  0.5 * sum(log(2 * pi * v) + (dat$ybar - log(ts))^2 / v) }
best <- NULL; set.seed(3)
for (k in 1:60) { f <- optim(c(rnorm(1, 0, .05), rnorm(1, -6, 3), log(.2)), nll_eta, control = list(maxit = 5000))
  if (is.null(best) || f$value < best$value) best <- f }
prof <- sapply(seq(-10, 2, by = .5), function(ll2) { b <- NULL
  for (s in c(-.1, -.01, 0, .01, .1)) { f <- optim(c(s, log(.2)), function(q) nll_eta(c(q[1], ll2, q[2])))
    if (is.null(b) || f$value < b$value) b <- f }; b$value })
eta_tab <- data.frame(objective = "maximise log egg rate (interior peak)", k = 3,
  AICc = round(2 * best$value + 6 + 24 / 24, 2), log_lambda2 = round(best$par[2], 2),
  dnll_trait_dominated = round(prof[1] - min(prof), 2))
print(eta_tab, row.names = FALSE)

draw <- read_lee("data/lee2008_nochoice.csv"); idx_by <- split(seq_len(nrow(draw)), paste(draw$Ratio, draw$Food))
bt <- do.call(rbind, mclapply(1:200, function(b) { set.seed(7000 + b)
  ii <- unlist(lapply(idx_by, function(ix) ix[sample.int(length(ix), replace = TRUE)])); db <- draw[ii, ]
  dt <- db[complete.cases(db[, c(XC, YC, "dailyeggs")]), ]
  g <- tryCatch(fit_glm(dt), error = function(e) NULL); if (is.null(g)) return(NULL)
  H <- hess_of(g); nd <- all(eigen(-H)$values > 0) && !g$optimum$used_fallback
  out <- c(interior = nd, theta_W = g$optimum$theta)
  for (n in mods) { fi <- fit_ing(prep_intake(db), fits_ing[[n]]$spec, n_start = 6, seed = b, start = fits_ing[[n]]$par)
    cm <- if (all(eigen(-H)$values > 0)) compare(fi$p$Qz, -H, fi$p$mu) else rep(NA, 5)
    out <- c(out, setNames(cm, paste0(n, ":", c("axis_ING", "axis_W", "axis_diff", "omega_ING", "omega_W")))) }
  out }, mc.cores = ING_CORES))
cat("\nBootstrap: interior log-scale peak in", round(100 * mean(bt[, "interior"]), 1), "% of replicates\n")
ci <- function(v) sprintf("%.3g [%.3g, %.3g]", median(v, na.rm = TRUE), quantile(v, .025, na.rm = TRUE), quantile(v, .975, na.rm = TRUE))
tab <- do.call(rbind, lapply(mods, function(n) data.frame(ING_model = n,
  axis_ING_deg = sprintf("%.1f", pt["axis_ING", n]), axis_eggrate_deg = sprintf("%.1f", pt["axis_W", n]),
  axis_difference_deg = paste0(sprintf("%.1f", pt["axis_diff", n]), "; boot ", ci(bt[, paste0(n, ":axis_diff")])),
  omega_ING = paste0(sprintf("%.3f", pt["omega_ING", n]), "; boot ", ci(bt[, paste0(n, ":omega_ING")])),
  omega_eggrate = paste0(sprintf("%.3f", pt["omega_W", n]), "; boot ", ci(bt[, paste0(n, ":omega_W")])))))
print(tab, row.names = FALSE)
write.csv(tab, "out/table_eggrate_hessian.csv", row.names = FALSE); write.csv(eta_tab, "out/table_eggrate_adaptive.csv", row.names = FALSE)
saveRDS(list(glm = g0, H = H0, point = pt, boot = bt, eta = eta_tab, eta_profile = prof), "out/eggrate_hessian.rds")
