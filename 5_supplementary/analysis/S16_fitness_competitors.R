# Appendix S5: fitness-based objectives with nutrient costs (Table S18)
source("R/ing_model.R"); source("R/table_utils.R"); suppressMessages(library(nutrigon)); library(parallel)
F0 <- readRDS("out/fits_nochoice.rds"); dat <- F0$dat; fits_ing <- F0$fits
data(lee_drosophila); Lf <- lee_drosophila; Lf$key <- paste(Lf$Ratio, Lf$Food)
stopifnot(all(Lf$key %in% dat$tr$key))
traits <- c("lifespan", "lifetimeegg", "dailyeggs")
T_MAX <- 2 * max(exp(dat$ybar))
set.seed(16)

surf <- function(y, keys = dat$tr$key) {
  d <- Lf[Lf$key %in% keys, ]
  cf <- unname(coef(lm(d[[y]] ~ protein_eaten + carb_eaten + I(protein_eaten^2) + I(carb_eaten^2) + protein_eaten:carb_eaten, d)))
  setNames(cf[2:6] / sd(d[[y]]), c("b1", "b2", "b3", "b4", "b5")) }
ZERO <- c(b1 = 0, b2 = 0, b3 = 0, b4 = 0, b5 = 0)

MODELS <- list(A = c("l1", "ll2"), B = c("l1", "ll2", "aP", "aC"), C = c("l1", "ll2", "aP", "aC", "lbP", "lbC"),
               D = c("l1", "ll2", "aP", "aC", "lbP", "lbC", "lVm"),
               E = c("l1", "aP", "aC", "lbP", "lbC"), F = c("l1", "aP", "aC", "lbP", "lbC", "lVm"))
unpackW <- function(q, m) { p <- setNames(as.list(q), MODELS[[m]])
  list(l1 = p$l1, l2 = if (is.null(p$ll2)) 1 else exp(p$ll2), aP = if (is.null(p$aP)) 0 else p$aP, aC = if (is.null(p$aC)) 0 else p$aC,
       bP = if (is.null(p$lbP)) 0 else exp(p$lbP), bC = if (is.null(p$lbC)) 0 else exp(p$lbC),
       Vm = if (is.null(p$lVm)) Inf else exp(p$lVm)) }
rail <- function(b, p, comp, Fv) sapply(seq_len(nrow(comp)), function(i) {
  c1 <- comp[i, 1]; c2 <- comp[i, 2]; g <- sum(W_FOOD * comp[i, ]) / Fv[i]
  a2 <- -(b[["b3"]] * c1^2 + b[["b4"]] * c2^2 + b[["b5"]] * c1 * c2) + p$l2 * g^2 + p$bP * c1^2 + p$bC * c2^2
  a1 <- -(b[["b1"]] * c1 + b[["b2"]] * c2) + p$l1 * g + p$aP * c1 + p$aC * c2
  tm <- min(T_MAX, p$Vm / g); cand <- c(1e-3, tm)
  if (a2 > 0) cand <- c(cand, min(max(-a1 / (2 * a2), 1e-3), tm))
  cand[which.min(a2 * cand^2 + a1 * cand)] })
nllW <- function(q, m, b, dd) { p <- unpackW(q[-length(q)], m); ts <- rail(b, p, dd$comp, dd$F)
  v <- exp(2 * q[length(q)]) + dd$s2n; 0.5 * sum(log(2 * pi * v) + (dd$ybar - log(ts))^2 / v) }
start_of <- function(m) { s <- c(l1 = 0, ll2 = log(.01), aP = 0, aC = 0, lbP = log(1e-4), lbC = log(1e-4), lVm = log(2))
  s[MODELS[[m]]] }
fitW <- function(m, b, dd, n_start = 24, st = NULL) { best <- NULL
  for (k in seq_len(n_start)) {
    q0 <- if (!is.null(st) && k == 1) st else c(start_of(m) + rnorm(length(MODELS[[m]]), 0, c(l1 = .3, ll2 = 2, aP = .3, aC = .3, lbP = 2, lbC = 2, lVm = .3)[MODELS[[m]]]), log(.15))
    f <- try(optim(q0, nllW, m = m, b = b, dd = dd, control = list(maxit = 6000, reltol = 1e-11)), silent = TRUE)
    if (inherits(f, "try-error")) next
    f2 <- try(optim(f$par, nllW, m = m, b = b, dd = dd, method = "BFGS"), silent = TRUE)
    if (!inherits(f2, "try-error") && f2$value < f$value) f <- f2
    if (is.null(best) || f$value < best$value) best <- f }
  k <- length(best$par); list(par = best$par, nll = best$value, k = k, AICc = 2 * best$value + 2 * k + 2 * k * (k + 1) / (28 - k - 1)) }
implied_target <- function(b, p) {
  H <- 2 * (-matrix(c(b[["b3"]], b[["b5"]] / 2, b[["b5"]] / 2, b[["b4"]]), 2) + diag(c(p$bP, p$bC)))
  if (any(eigen(H, symmetric = TRUE)$values <= 0)) return(c(P = NA, C = NA))
  x <- solve(H, c(b[["b1"]], b[["b2"]]) - c(p$aP, p$aC)); setNames(x, c("P", "C")) }
cv_err <- function(m, y, full_par, fold) {
  groups <- if (fold == "rail") dat$tr$Ratio else dat$tr$F
  e <- c(); for (g in unique(groups)) { tin <- which(groups != g); tout <- which(groups == g)
    b <- if (y == "none") ZERO else surf(y, dat$tr$key[tin]); dd <- subset_trt(dat, tin)
    f <- fitW(m, b, dd, n_start = 8, st = full_par)
    pr <- rail(b, unpackW(f$par[-length(f$par)], m), dat$comp[tout, , drop = FALSE], dat$F[tout])
    e <- c(e, dat$ybar[tout] - log(pr)) }
  100 * (exp(sqrt(mean(e^2))) - 1) }

jobs <- c(lapply(traits, function(y) lapply(c("A", "B", "C", "D"), function(m) list(y = y, m = m))), list(list(list(y = "none", m = "E"), list(y = "none", m = "F"))))
jobs <- unlist(jobs, recursive = FALSE)
res <- mclapply(jobs, function(j) { b <- if (j$y == "none") ZERO else surf(j$y)
  f <- fitW(j$m, b, dat); p <- unpackW(f$par[-length(f$par)], j$m); tg <- implied_target(b, p)
  leak <- if (j$m == "A") { e <- c(); for (g in unique(dat$tr$Ratio)) { tin <- which(dat$tr$Ratio != g); tout <- which(dat$tr$Ratio == g)
      ff <- fitW("A", b, subset_trt(dat, tin), 8, f$par); e <- c(e, dat$ybar[tout] - log(rail(b, unpackW(ff$par[-length(ff$par)], "A"), dat$comp[tout, , drop = FALSE], dat$F[tout]))) }
      100 * (exp(sqrt(mean(e^2))) - 1) } else NA
  list(trait = j$y, model = j$m, fit = f, p = p, target = tg, LORO = cv_err(j$m, j$y, f$par, "rail"),
       LOCO = cv_err(j$m, j$y, f$par, "conc"), LORO_fullsurface = leak) }, mc.cores = ING_CORES)

lab <- c(A = "trait + volume cost", B = "A + linear nutrient costs", C = "B + quadratic nutrient costs", D = "C + capacity",
         E = "nutrient and volume costs only", F = "E + capacity")
tab <- do.call(rbind, lapply(res, function(r) data.frame(trait = r$trait, model = r$model, objective = lab[[r$model]], k = r$fit$k,
  nll = r$fit$nll, AICc = r$fit$AICc, LORO_pct = r$LORO, LOCO_pct = r$LOCO, LORO_fullsurface_pct = r$LORO_fullsurface,
  target_CP = unname(r$target["C"] / r$target["P"]), target_P = unname(r$target["P"]), target_C = unname(r$target["C"]))))
cv <- read.csv("out/table_cv.csv")
ing <- do.call(rbind, lapply(c(ING_BEST, "Aniso+L12", "Aniso+L12+Cap"), function(n) data.frame(trait = "ING", model = n, objective = paste("ING", n),
  k = fits_ing[[n]]$k, nll = fits_ing[[n]]$nll, AICc = fits_ing[[n]]$AICc, LORO_pct = cv$pct_rail[cv$model == n], LOCO_pct = cv$pct_conc[cv$model == n],
  LORO_fullsurface_pct = NA, target_CP = fits_ing[[n]]$desc[["C_per_P"]], target_P = fits_ing[[n]]$p$mu[1], target_C = fits_ing[[n]]$p$mu[2])))
tab <- rbind(tab, ing); tab$dAICc <- tab$AICc - min(tab$AICc)
nest <- do.call(rbind, lapply(traits, function(y) { g <- function(m, t) tab$nll[tab$trait == t & tab$model == m]
  data.frame(trait = y, LR_C_vs_E = 2 * (g("E", "none") - g("C", y)), LR_D_vs_F = 2 * (g("F", "none") - g("D", y))) }))
print(tab[order(tab$AICc), ], row.names = FALSE, digits = 3); print(nest, digits = 3)
saveRDS(list(res = res, tab = tab, nest = nest), "out/fitness_competitors.rds")

f1 <- function(x, d = 1) ifelse(is.na(x), "--", formatC(x, format = "f", digits = d))
tl <- c(lifespan = "lifespan", lifetimeegg = "lifetime eggs", dailyeggs = "egg-production rate", none = "--", ING = "--")
T16 <- tab[order(tab$AICc), ]
T16 <- data.frame(Objective = ifelse(T16$trait == "ING", T16$objective, paste0(T16$objective, ifelse(T16$trait == "none", "", paste0(" (", tl[T16$trait], ")")))),
                  k = T16$k, AICc = f1(T16$AICc), dAICc = f1(T16$dAICc), LORO = f1(T16$LORO_pct), LOCO = f1(T16$LOCO_pct), CP = f1(T16$target_CP, 2))
write_table(T16, "TableS18_fitness_competitors", dir = "tables/supplementary", escape = TRUE, label = "tab:S18",
  header = c("Objective", "$k$", "AICc", "$\\Delta$AICc", "\\makecell{LORO\\\\(\\%)}", "\\makecell{LOCO\\\\(\\%)}", "\\makecell{Implied\\\\target C:P}"),
  caption = paste("Fitness-based objectives with nutrient-specific costs compared with ING on the no-choice array of \\citet{Lee2008}.",
    "Each objective is $-W(\\mathbf{x})$ (a fitted performance surface, in trait standard deviations) plus a volume cost",
    "$\\lambda_1V+\\lambda_2V^2$, optionally with linear and quadratic costs for each nutrient and an intake capacity; the last",
    "two rows without a trait contain the costs only. LORO, LOCO: prediction error of treatment-mean intake for held-out rails",
    "and concentrations, with the performance surface refitted to the training flies in every fold. Implied target: the",
    "minimiser of the nutrient part of the objective (trait surface and nutrient costs), where it is a minimum; for ING, the",
    "latent target. Likelihood-ratio statistics for adding the trait surface to the costs:",
    paste(sprintf("%s %.1f (without capacity) and %.1f (with capacity)", tl[nest$trait], pmax(nest$LR_C_vs_E, 0), pmax(nest$LR_D_vs_F, 0)), collapse = "; "),
    "(values within optimisation tolerance of zero are shown as 0.0).",
    "With a performance surface, nutrient costs and an intake capacity, the objective is a constrained version of symmetric anisotropic",
    "ING with linear and quadratic volume costs and a capacity (Aniso+L12+Cap), and it reached the same likelihood."))
cat("S16 done\n")
