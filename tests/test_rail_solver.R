# Test: rail solver vs brute force
source("R/ing_model.R"); set.seed(1)
comp <- cbind(c(0, 1/16, 1/8, 1/4, 1/2, 1, 1.875), 1); comp <- comp / rowSums(comp)
p <- list(mu = c(60, 150), Qz = Qz_uv(0.5, -0.4), eta = c(P = 0.6, C = -0.3), l1 = 0.05, l2 = 0.4)
worst_brute <- 0; worst_conic <- 0; worst_CR <- 0
for (F in c(45, 90, 180, 360)) {
  ts <- rail_opt(p, comp, rep(F, 7))
  br <- sapply(1:7, function(i) optimize(function(t) { x <- matrix(t * comp[i, ], 1); V <- sum(W_FOOD * t * comp[i, ]) / F
    loss_R(x, p) + p$l1 * V + p$l2 * V^2 }, c(0, 3000), tol = 1e-10)$minimum)
  worst_brute <- max(worst_brute, abs(ts - br))
  worst_CR <- max(worst_CR, abs(ts - rail_opt_R(p, comp, rep(F, 7))))
  ps <- p; ps$eta <- c(P = 0, C = 0)
  X <- rail_opt(ps, comp, rep(F, 7)) * comp; ec <- effective_conic(ps, F)
  r <- apply(X, 1, function(x) drop(t(x - ec$centre) %*% ec$Qe %*% (x - ec$centre))) / drop(t(ec$centre) %*% ec$Qe %*% ec$centre) - 1
  worst_conic <- max(worst_conic, abs(r))
}
cat(sprintf("rail solver vs brute force: max |dt| = %.2e\nsymmetric optima vs conic: max relative residual = %.2e\nC vs R implementation: max |dt| = %.2e (C used: %s)\n",
            worst_brute, worst_conic, worst_CR, ING_USE_C))
stopifnot(worst_brute < 1e-4, worst_conic < 1e-10, worst_CR < 1e-8)
cat("PASS\n")

set.seed(21); worst_cap <- -Inf; worst_capCR <- 0
comp7 <- cbind(c(0, 1/16, 1/8, 1/4, 1/2, 1, 1.875), 1); comp7 <- comp7 / rowSums(comp7)
for (r in 1:200) { s <- spec(TRUE, c("P", "C"), TRUE, TRUE, cap = TRUE)
  p <- unpack(c(log(c(80, 128)) + rnorm(2, 0, .3), rnorm(2, 0, .7), rnorm(2, 0, .8), rnorm(1, 0, .2),
                log(runif(1, .01, .2)), log(runif(1, .3, 2.5))), s)
  for (F in c(45, 90, 180, 360)) { g <- as.vector(comp7 %*% W_FOOD) / F
    tc <- rail_opt(p, comp7, rep(F, 7)); worst_capCR <- max(worst_capCR, abs(tc - rail_opt_R(p, comp7, rep(F, 7))))
    for (i in 1:7) { tt <- seq(0, p$vmax / g[i], length.out = 4001)
      J <- function(t) loss_R(outer(t, comp7[i, ]), p) + p$l1 * g[i] * t + p$l2 * (g[i] * t)^2
      stopifnot(tc[i] <= p$vmax / g[i] * (1 + 1e-12))
      worst_cap <- max(worst_cap, J(tc[i]) - min(J(tt))) } } }
cat(sprintf("with capacity: exact minus brute-force objective (max) = %.2e; C vs R max |dt| = %.2e\n", worst_cap, worst_capCR))
stopifnot(worst_cap <= 1e-12, worst_capCR < 1e-8)
cat("PASS (capacity)\n")
