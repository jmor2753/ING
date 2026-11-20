# Appendix S5: profile likelihoods and bootstrap correlations (Table S16, Fig. S13)
source("R/ing_model.R"); source("R/table_utils.R"); library(parallel)
F0 <- readRDS("out/fits_nochoice.rds"); dat <- F0$dat
best <- F0$fits[[ING_BEST]]; s <- best$spec; par0 <- best$par
CB <- readRDS("out/cv_boot.rds"); BT <- CB$boots[[ING_BEST]]
PN <- c("log muP", "log muC", "a", "b", "etaP", "etaC", "log lambda2", "log Vmax", "log tau")
stopifnot(length(par0) == length(PN))

prof_at <- function(j, v, starts) {
  f <- function(q) { p <- par0; p[-j] <- q; p[j] <- v; negll(p, dat, s) }
  bv <- Inf; bq <- starts[[1]]
  for (st in starts) {
    o <- try(optim(st, f, control = list(maxit = 4000, reltol = 1e-11)), silent = TRUE)
    if (inherits(o, "try-error")) next
    o2 <- try(optim(o$par, f, method = "BFGS", control = list(maxit = 2000)), silent = TRUE)
    if (!inherits(o2, "try-error") && o2$value < o$value) o <- o2
    if (o$value < bv) { bv <- o$value; bq <- o$par } }
  list(nll = bv, q = bq) }
GRID <- list(etaP = seq(-4, 2, by = .25), etaC = seq(-2.5, 4, by = .25),
             a = seq(0, 1.4, by = .05), b = seq(-0.5, 0.3, by = .025),
             `log lambda2` = seq(-7, -1, by = .25), `log Vmax` = seq(0.2, 1.2, by = .05))
GRID <- lapply(setNames(names(GRID), names(GRID)), function(n) sort(unique(c(GRID[[n]], par0[match(n, PN)]))))
profile_one <- function(n) {
  j <- match(n, PN); g <- GRID[[n]]; k0 <- which.min(abs(g - par0[j])); out <- rep(Inf, length(g))
  set.seed(17 + j)
  for (dir in list(k0:length(g), k0:1)) { q <- par0[-j]
    for (k in dir) { st <- list(q, par0[-j], par0[-j] + rnorm(length(q), 0, .3), q + rnorm(length(q), 0, .15))
      r <- prof_at(j, g[k], st); out[k] <- min(out[k], r$nll); if (r$nll < 1e9) q <- r$q } }
  data.frame(par = n, v = g, nll = out) }
jobs <- do.call(rbind, mclapply(names(GRID), profile_one, mc.cores = ING_CORES))
nll0 <- min(best$nll, jobs$nll); jobs$dll <- jobs$nll - nll0
cat("global nll", best$nll, " min profile nll", min(jobs$nll), "\n")

prof_ci <- function(n) {
  d <- jobs[jobs$par == n, ]; d <- d[is.finite(d$nll), ]; k <- which.min(d$dll)
  cut <- 1.92; lo <- hi <- NA
  if (k > 1) { i <- max(which(d$dll[1:k] > cut), -Inf)
    lo <- if (is.finite(i)) approx(d$dll[i:(i + 1)], d$v[i:(i + 1)], cut)$y else -Inf }
  if (k < nrow(d)) { i <- k - 1 + min(which(d$dll[k:nrow(d)] > cut), Inf)
    hi <- if (is.finite(i)) approx(d$dll[(i - 1):i], d$v[(i - 1):i], cut)$y else Inf }
  c(est = par0[match(n, PN)], lo = lo, hi = hi) }
CI <- t(sapply(names(GRID), prof_ci))

bp <- as.matrix(BT[, paste0(".par", seq_along(PN))]); colnames(bp) <- PN
bp <- bp[stats::complete.cases(bp), , drop = FALSE]
CORR <- cor(bp[, -length(PN)], method = "spearman")
BCI <- t(apply(bp, 2, quantile, c(.025, .975)))
saveRDS(list(jobs = jobs, ci = CI, boot_ci = BCI, corr = CORR, n_boot = nrow(bp)), "out/param_identifiability.rds")

nat <- list(etaP = function(x) exp(2 * x), etaC = function(x) exp(2 * x), a = identity, b = identity,
            `log lambda2` = function(x) exp(x), `log Vmax` = function(x) exp(x))
lab <- c(etaP = "Protein excess:deficit weight $e^{2\\eta_P}$", etaC = "Carbohydrate excess:deficit weight $e^{2\\eta_C}$",
         a = "Metric parameter $a$", b = "Metric parameter $b$",
         `log lambda2` = "Quadratic volume cost $\\lambda_2$", `log Vmax` = "Capacity $V_{\\max}$ ($\\mu$L day$^{-1}$)")
fmt <- function(x) ifelse(is.infinite(x), ifelse(x > 0, "$\\infty$", "0"), formatC(x, digits = 3, format = "g", flag = "#"))
fmtg <- function(x, f) { y <- f(x); if (is.infinite(x)) { if (identical(f, identity)) ifelse(x > 0, "$+\\infty$", "$-\\infty$") else fmt(y) } else formatC(y, digits = 3, format = "g", flag = "#") }
open_end <- function(n, side) { g <- GRID[[n]]; v <- CI[n, side]; is.na(v) || is.infinite(v) }
T18 <- do.call(rbind, lapply(names(GRID), function(n) { f <- nat[[n]]
  lo <- if (open_end(n, "lo")) paste0("$<$", formatC(f(min(GRID[[n]])), digits = 3, format = "g")) else fmtg(CI[n, "lo"], f)
  hi <- if (open_end(n, "hi")) paste0("$>$", formatC(f(max(GRID[[n]])), digits = 3, format = "g")) else fmtg(CI[n, "hi"], f)
  data.frame(par = lab[[n]], est = fmtg(CI[n, "est"], f), prof = paste0("[", lo, ", ", hi, "]"),
             boot = sprintf("[%s, %s]", fmtg(BCI[n, 1], f), fmtg(BCI[n, 2], f)),
             excl = if (n %in% c("etaP", "etaC")) { if (!open_end(n, "lo") && CI[n, "lo"] > 0) "yes ($>1$)" else if (!open_end(n, "hi") && CI[n, "hi"] < 0) "yes ($<1$)" else "no" } else "") }))
print(T18, row.names = FALSE)
write_table(T18, "TableS16_profile_intervals", dir = "tables/supplementary", escape = FALSE, label = "tab:S16",
  header = c("Parameter", "Estimate", "95\\% profile interval", "95\\% bootstrap interval", "Symmetry excluded"),
  caption = paste0("Profile-likelihood and bootstrap intervals for the components of the best model (", ING_BEST,
    ") fitted to the no-choice data of \\citet{Lee2008}. Profiles re-estimate all other parameters at each fixed value; ",
    "the interval is where the log likelihood is within 1.92 of its maximum. An open bound ($<$ or $>$) means that the ",
    "profile had not dropped by 1.92 at the edge of the grid. Symmetry excluded: whether the interval for the ",
    "excess:deficit weight excludes 1 (equal penalties for excess and deficit)."))

dir.create("figures/supplementary", showWarnings = FALSE, recursive = TRUE)
mk <- function() {
  layout(matrix(c(1, 2, 7, 3, 4, 7, 5, 6, 7), 3, byrow = TRUE), widths = c(1, 1, 1.4))
  par(mar = c(4, 4, 1.5, 1), mgp = c(2.3, .7, 0), las = 1, cex = .8)
  xl <- c(etaP = expression(eta[P]), etaC = expression(eta[C]), a = "a", b = "b",
          `log lambda2` = expression(log~lambda[2]), `log Vmax` = expression(log~V[max]))
  for (n in names(GRID)) { d <- jobs[jobs$par == n, ]
    plot(d$v, pmin(d$dll, 10), type = "l", lwd = 2, col = "#1f4e79", xlab = xl[[n]], ylab = expression(Delta~"neg. log L"), ylim = c(0, 10))
    abline(h = 1.92, lty = 2, col = "grey40"); abline(v = par0[match(n, PN)], col = "#b45309", lty = 3)
    if (n %in% c("etaP", "etaC")) abline(v = 0, col = "grey60") }
  C <- CORR; k <- ncol(C); par(mar = c(6.5, 6.5, 1.5, 3))
  image(1:k, 1:k, C[, k:1], zlim = c(-1, 1), col = hcl.colors(41, "Blue-Red 3", rev = TRUE), axes = FALSE, xlab = "", ylab = "")
  axis(1, 1:k, colnames(C), las = 2, cex.axis = .8); axis(2, k:1, colnames(C), las = 1, cex.axis = .8)
  for (i in 1:k) for (j in 1:k) text(i, k + 1 - j, sprintf("%.2f", C[i, j]), cex = .65, col = ifelse(abs(C[i, j]) > .6, "white", "black"))
  box() }
pdf("figures/supplementary/FigS13_profiles_correlation.pdf", 9.5, 7); mk(); dev.off()
png("figures/supplementary/FigS13_profiles_correlation.png", 9.5, 7, units = "in", res = 200); mk(); dev.off()
cc <- which(abs(CORR) > .6 & upper.tri(CORR), arr.ind = TRUE)
if (nrow(cc)) for (r in seq_len(nrow(cc))) cat(sprintf("corr %s ~ %s: %.2f\n", rownames(CORR)[cc[r, 1]], colnames(CORR)[cc[r, 2]], CORR[cc[r, 1], cc[r, 2]]))
cat("S17 done\n")
