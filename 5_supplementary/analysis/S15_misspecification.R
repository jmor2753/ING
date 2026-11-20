# Appendix S5: non-quadratic regulatory losses (Table S17)
source("R/ing_model.R"); source("R/table_utils.R"); library(parallel)
F0 <- readRDS("out/fits_nochoice.rds"); dat <- F0$dat; fits <- F0$fits
args <- commandArgs(TRUE); R_REP <- if (length(args)) as.integer(args[1]) else 100
best <- fits[[ING_BEST]]; p0 <- best$p; tau0 <- best$tau
CP_TRUE <- p0$mu[2] / p0$mu[1]

dev_of <- function(x, mu, eta, scale) {
  r <- sweep(x, 2, mu, "/"); z <- if (scale == "bc05") 2 * (sqrt(r) - 1) else r - 1
  z * ifelse(z > 0, matrix(exp(eta), nrow(z), 2, byrow = TRUE), 1) }
LOSSES <- list(
  quadratic = function(x, p) { d <- dev_of(x, p$mu, p$eta, "prop"); rowSums((d %*% p$Qz) * d) },
  power1    = function(x, p) { d <- dev_of(x, p$mu, p$eta, "prop"); sqrt(rowSums((d %*% p$Qz) * d)) },
  power3    = function(x, p) { d <- dev_of(x, p$mu, p$eta, "prop"); rowSums((d %*% p$Qz) * d)^1.5 },
  huber     = function(x, p) { d <- dev_of(x, p$mu, p$eta, "prop"); r <- sqrt(rowSums((d %*% p$Qz) * d)); ifelse(r <= .5, r^2, .5 * (2 * r - .5)) },
  boxcox05  = function(x, p) { d <- dev_of(x, p$mu, p$eta, "bc05"); rowSums((d %*% p$Qz) * d) },
  linex     = function(x, p) { d <- dev_of(x, p$mu, p$eta, "prop"); w <- diag(p$Qz); rowSums(sweep((exp(2 * d) - 2 * d - 1) / 2, 2, w, "*")) })

rail_gen <- function(loss, p, kappa, comp, Fv, ng = 400) {
  g <- as.vector(comp %*% W_FOOD) / Fv
  sapply(seq_len(nrow(comp)), function(i) {
    tmax <- min(vmax_of(p) / g[i], 6 * max(p$mu / pmax(comp[i, ], 1e-9)))
    J <- function(t) kappa * loss(outer(t, comp[i, ]), p) + p$l1 * g[i] * t + p$l2 * g[i]^2 * t^2
    tt <- seq(1e-3, tmax, length.out = ng); k <- which.min(J(tt))
    o <- optimize(J, c(tt[max(k - 1, 1)], tt[min(k + 1, ng)]), tol = 1e-9)
    c(tt[k], o$minimum)[which.min(c(J(tt[k]), o$objective))] }) }
calib <- function(loss) {
  if (identical(loss, LOSSES$quadratic)) return(c(kappa = 1, s = 1))
  f <- function(q) { p <- p0; p$mu <- exp(q[2]) * p0$mu
    sum((log(rail_gen(loss, p, exp(q[1]), dat$comp, dat$F, 200)) - dat$ybar)^2) }
  o <- optim(c(0, 0), f, control = list(reltol = 1e-8, maxit = 300))
  c(kappa = exp(o$par[1]), s = exp(o$par[2])) }
CAL <- sapply(LOSSES, calib); KAPPA <- CAL["kappa", ]
P_GEN <- lapply(names(LOSSES), function(n) { p <- p0; p$mu <- CAL["s", n] * p0$mu; p }); names(P_GEN) <- names(LOSSES)
TRUE_T <- lapply(names(LOSSES), function(n) rail_gen(LOSSES[[n]], P_GEN[[n]], KAPPA[[n]], dat$comp, dat$F)); names(TRUE_T) <- names(LOSSES)
GEN_RMSE <- sapply(TRUE_T, function(t) sqrt(mean((log(t) - dat$ybar)^2)))
print(rbind(CAL, rmse = GEN_RMSE))

cands <- c(ING_BEST, "Aniso+L2", "Iso+L1")
one_rep <- function(sc, r) {
  set.seed(5000 * match(sc, names(LOSSES)) + r)
  lt <- log(TRUE_T[[sc]])[dat$trt] + rnorm(length(dat$trt), 0, dat$s_i[dat$trt]) + rnorm(length(dat$ybar), 0, tau0)[dat$trt]
  db <- dat; db$lt <- lt; db$ybar <- as.vector(tapply(lt, dat$trt, mean)); db$s2n <- as.vector(tapply(lt, dat$trt, var) / tabulate(dat$trt))
  fs <- lapply(cands, function(n) fit_ing(db, fits[[n]]$spec, n_start = 8, seed = r, start = fits[[n]]$par))
  a <- sapply(fs, `[[`, "AICc")
  c(sel = which.min(a), setNames(sapply(fs, function(f) f$desc[["C_per_P"]]), cands), aicc_best = a[1]) }
out <- list()
for (sc in names(LOSSES)) {
  rows <- mclapply(seq_len(R_REP), function(r) tryCatch(one_rep(sc, r), error = function(e) NULL), mc.cores = ING_CORES)
  out[[sc]] <- do.call(rbind, rows); cat(sc, "done", nrow(out[[sc]]), "\n")
  saveRDS(list(out = out, cal = CAL, gen_rmse = GEN_RMSE, true_t = TRUE_T, cp_true = CP_TRUE), "out/misspecification.rds") }

dist <- sapply(names(LOSSES), function(n) 100 * (exp(sqrt(mean((log(TRUE_T[[n]]) - log(TRUE_T$quadratic))^2))) - 1))
lab <- c(quadratic = "Quadratic (reference)", power1 = "Power 1 (linear in distance)", power3 = "Power 3 (light tails)",
         huber = "Huber (quadratic core, linear tails)", boxcox05 = "Box--Cox deviations, $\\beta=0.5$", linex = "LINEX (exponential tails)")
q <- function(v) sprintf("%.2f [%.2f, %.2f]", median(v), quantile(v, .025), quantile(v, .975))
T17 <- do.call(rbind, lapply(names(out), function(sc) { o <- out[[sc]]
  data.frame(loss = lab[[sc]], dist = round(dist[[sc]], 1), fit = sprintf("%.3f", GEN_RMSE[[sc]]), n = nrow(o),
             best = q(o[, ING_BEST]), aniso = q(o[, "Aniso+L2"]), iso = q(o[, "Iso+L1"]),
             within = round(100 * mean(o[, ING_BEST] >= 1.32 & o[, ING_BEST] <= 1.89)),
             rej4 = round(100 * mean(o[, ING_BEST] < 3)),
             sel_best = round(100 * mean(o[, "sel"] == 1))) }))
print(T17, row.names = FALSE)
write_table(T17, "TableS17_misspecification", dir = "tables/supplementary", escape = FALSE, label = "tab:S17",
  header = c("True regulatory loss", "\\makecell{Array\\\\change (\\%)}", "\\makecell{RMSE vs\\\\data}", "Data sets", "\\makecell{C:P, best model\\\\median [95\\%]}",
             "\\makecell{C:P,\\\\Aniso+L2}", "\\makecell{C:P,\\\\Iso+L1}", "\\makecell{Within 1.32--1.89\\\\(\\%)}",
             "\\makecell{C:P $<3$\\\\(\\%)}", "\\makecell{Best model\\\\selected (\\%)}"),
  caption = paste0("Recovery of the target ratio when the true regulatory loss is not quadratic. Data were simulated on the ",
    "design of \\citet{Lee2008} from the target (C:P ", sprintf("%.2f", CP_TRUE), "), metric, excess--deficit weights, volume cost, ",
    "capacity and noise of the best model, with the loss replaced as indicated, and fitted with three ING models that assume a ",
    "quadratic loss. Each non-quadratic loss was scaled, and its target magnitude (not ratio) adjusted, by least squares on the observed treatment means. Array change: root-mean-square difference between the noise-free array and that of the quadratic loss. RMSE vs data: root-mean-square difference on the log scale between the noise-free array and the observed treatment means. ",
    "Within 1.32--1.89: percentage of estimates inside the 95\\% interval of the empirical estimate. C:P $<3$: percentage of ",
    "estimates below 3, well below the self-selected ratio (3.5). Best model selected: percentage of data sets in ",
    "which the best model had the lowest AICc of the three."))
cat("S15 done\n")
