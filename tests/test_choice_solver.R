# Test: choice-arena solver vs numerical optimisation
source("R/ing_model.R"); set.seed(5)
J_at <- function(p, FY, FS, r, sep) { aY <- FY * YEAST; x <- matrix(r[3] * aY + c(0, FS) * r[4], 1)
  V <- r[3] + r[4]; loss_R(x, p) + p$l1 * V + p$l2 * (if (sep) r[3]^2 + r[4]^2 else V^2) }
worse <- 0; dCR <- 0
for (k in 1:300) {
  p <- list(mu = exp(rnorm(2, c(3.5, 4.5), .6)), Qz = Qz_uv(rnorm(1), rnorm(1)),
            eta = c(P = rnorm(1, 0, .8), C = rnorm(1, 0, .8)), l1 = rnorm(1, 0, .1), l2 = exp(rnorm(1, -2.5, 1)))
  FY <- sample(c(45, 90, 180), 1); FS <- sample(c(45, 90, 180), 1)
  for (sep in c(FALSE, TRUE)) {
    e <- choice_opt(p, FY, FS, sep); n <- choice_opt_num(p, FY, FS, sep)
    worse <- max(worse, J_at(p, FY, FS, e, sep) - n[["J"]])
    dCR <- max(dCR, abs(e - choice_opt_R(p, FY, FS, sep)))
  }
}
cat(sprintf("exact minus numerical objective (max over 600 problems): %.2e  (<= ~0: exact never worse)\nC vs R implementation: max difference %.2e (C used: %s)\n",
            worse, dCR, ING_USE_C))
stopifnot(worse < 1e-8, dCR < 1e-6)
cat("PASS\n")

set.seed(22); worst_cap <- -Inf; capCR <- 0
FY <- rep(c(180, 90, 45), 3); FS <- rep(c(180, 90, 45), each = 3)
for (r in 1:150) { s <- spec(TRUE, c("P", "C"), TRUE, TRUE, cap = TRUE)
  p <- unpack(c(log(c(80, 128)) + rnorm(2, 0, .3), rnorm(2, 0, .7), rnorm(2, 0, .8), rnorm(1, 0, .2),
                log(runif(1, .01, .2)), log(runif(1, .2, 2))), s)
  for (sep in c(FALSE, TRUE)) {
    cc <- choice_opt_vec(p, FY, FS, sep); capCR <- max(capCR, max(abs(cc - t(sapply(1:9, function(i) choice_opt_R(p, FY[i], FS[i], sep))))))
    for (i in c(1, 5, 9)) { A <- cbind(FY[i] * YEAST, c(0, FS[i])); VM <- p$vmax
      G <- as.matrix(expand.grid(seq(0, VM, length.out = 250), seq(0, VM, length.out = 250))); if (!sep) G <- G[rowSums(G) <= VM, ]
      Jm <- function(G) { V <- rowSums(G); loss_R(G %*% t(A), p) + p$l1 * V + p$l2 * (if (sep) rowSums(G^2) else V^2) }
      worst_cap <- max(worst_cap, Jm(matrix(cc[i, 3:4], 1)) - min(Jm(G))) } } }
cat(sprintf("with capacity: exact minus grid objective (max) = %.2e; C vs R max difference = %.2e\n", worst_cap, capCR))
stopifnot(worst_cap <= 1e-10, capCR < 1e-6)
cat("PASS (capacity)\n")
