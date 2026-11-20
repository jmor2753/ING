# Non-linear extensions (Box-Cox scale, palatability)
source("R/ing_model.R"); source("R/ing_ext.R"); library(parallel)
F0 <- readRDS("out/fits_nochoice.rds"); dat <- F0$dat; fits <- F0$fits
dat$S <- sucrose_conc(dat$comp, dat$F)
ch <- load_choice(6)
FY <- ch$yeast_gL; FS <- ch$sucrose_gL; lPo <- log(ch$P_ug_day); lCo <- log(ch$C_ug_day); lro <- lCo - lPo

PC <- c("P", "C")
S <- list(
  "Aniso+AsymPC+L2+Cap"        = spec_ext(asym = PC, cap = TRUE),
  "Aniso+AsymPC+L2+Cap+BC"     = spec_ext(asym = PC, cap = TRUE, bc = TRUE),
  "Aniso+AsymPC+L2+Cap+Pal"    = spec_ext(asym = PC, cap = TRUE, pal = TRUE),
  "Aniso+AsymPC+L2+Cap+BC+Pal" = spec_ext(asym = PC, cap = TRUE, bc = TRUE, pal = TRUE),
  "Aniso+AsymPC+L2"         = spec_ext(asym = PC),
  "Aniso+AsymPC+L2+BC"      = spec_ext(asym = PC, bc = TRUE),
  "Aniso+AsymPC+L2+Pal"     = spec_ext(asym = PC, pal = TRUE),
  "Aniso+L2"                = spec_ext(),
  "Aniso+L2+BC"             = spec_ext(bc = TRUE),
  "Aniso+L2+Pal"            = spec_ext(pal = TRUE))
stopifnot(ING_BEST %in% names(S))

s_of <- function(n) { f <- fits[[n]]; s <- f$spec; s$bc <- FALSE; s$pal <- FALSE; list(par = f$par, s = s) }
starts_for <- function(s, from) {
  base <- lapply(from, function(b) embed_par(b$par, b$s, s)); out <- base
  nb <- n_par(s)
  if (s$bc) for (b in base) for (lb in log(c(.5, .2))) { v <- b; v[nb + 1] <- lb; out[[length(out) + 1]] <- v }
  if (s$pal) { j <- nb + s$bc + 1
    for (b in base) for (g in c(.1, -.1)) for (lk in log(c(.5, 2))) { v <- b; v[j] <- g; v[j + 1] <- lk; out[[length(out) + 1]] <- v } }
  out
}
base_from <- function(s) list(s_of(spec(s$aniso, s$asym, s$l1, s$l2, cap = isTRUE(s$cap))$name))

cat("Fitting", length(S), "models\n")
stage1 <- names(S)[!grepl("BC\\+Pal", names(S))]
E <- mclapply(stage1, function(n) { s <- S[[n]]
  fit_ext(dat, s, starts_for(s, base_from(s)), n_random = 8, seed = 7) }, mc.cores = ING_CORES)
names(E) <- stage1
for (n in names(S)[grepl("BC\\+Pal", names(S))]) {
  s <- S[[n]]; a <- sub("\\+Pal$", "", n); b <- sub("\\+BC", "", n)
  par_from <- list(list(par = E[[a]]$par, s = S[[a]]), list(par = E[[b]]$par, s = S[[b]]))
  E[[n]] <- fit_ext(dat, s, starts_for(s, par_from), n_random = 8, seed = 7)
}
E <- E[names(S)]
for (n in names(E)) cat(sprintf("%-26s k=%2d  nll=%8.3f  AICc=%7.2f\n", n, E[[n]]$k, E[[n]]$nll, E[[n]]$AICc))

subset_ext <- function(keep) { d <- subset_trt(dat, keep); d$S <- dat$S[keep]; d }
cv_model <- function(n) {
  f <- E[[n]]; s <- S[[n]]; rail <- as.character(dat$tr$Ratio); conc <- dat$F
  one <- function(test) { tr <- setdiff(seq_along(dat$ybar), test)
    g <- fit_ext(subset_ext(tr), s, list(f$par), n_random = 2, seed = 3, jitter = .2)
    log(pmax(rail_opt_ext(g$p, dat$comp[test, , drop = FALSE], dat$F[test], dat$S[test]), 1e-6)) - dat$ybar[test] }
  eR <- unlist(lapply(unique(rail), function(r) one(which(rail == r))))
  eC <- unlist(lapply(unique(conc), function(cc) one(which(conc == cc))))
  pc <- function(e) 100 * (exp(sqrt(mean(e^2))) - 1)
  nR <- sapply(unique(rail), function(r) sum(rail == r)); no01 <- rep(unique(rail) != "(0:1)", nR)
  c(LORO_pct = pc(eR), LORO_no01_pct = pc(eR[no01]), LOCO_pct = pc(eC))
}
cat("Cross-validation...\n")
CV <- mclapply(names(E), cv_model, mc.cores = ING_CORES); names(CV) <- names(E)

score_choice <- function(p) {
  ss <- function(lr) { pc <- p; pc$mu <- exp(lr) * p$mu; pr <- choice_ext(pc, FY, FS)
    if (any(pr[, 1:2] <= 1e-8)) return(1e10); sum((lPo - log(pr[, 1]))^2 + (lCo - log(pr[, 2]))^2) }
  o <- optimize(ss, c(log(.1), log(3))); pc <- p; pc$mu <- exp(o$minimum) * p$mu
  pr <- choice_ext(pc, FY, FS); lr <- log(pr[, 2] / pr[, 1])
  list(pred = pr, lr = lr, rho = exp(o$minimum), ratio_RMSE = sqrt(mean((lr - lro)^2)), ratio_r = cor(lr, lro),
       slope = unname(coef(lm(lro ~ lr))[2]), CP_S180 = mean(exp(lr[FS == 180])), amount_RMSE = sqrt(o$objective / 18))
}
cat("Choice predictions...\n")
CH <- mclapply(names(E), function(n) score_choice(E[[n]]$p), mc.cores = ING_CORES); names(CH) <- names(E)
chk <- max(abs(choice_ext(E[[ING_BEST]]$p, FY, FS)[, 1:2] / choice_opt_vec(E[[ING_BEST]]$p, FY, FS)[, 1:2] - 1))
cat(sprintf("numerical vs exact choice solver: %.1e\n", chk))

tab <- do.call(rbind, lapply(names(E), function(n) { f <- E[[n]]; d <- f$desc; c_ <- CH[[n]]
  data.frame(model = n, k = f$k, AICc = round(f$AICc, 2), tau = round(f$tau, 3),
    LORO_pct = round(CV[[n]][["LORO_pct"]], 1), LORO_no01_pct = round(CV[[n]][["LORO_no01_pct"]], 1), LOCO_pct = round(CV[[n]][["LOCO_pct"]], 1),
    C_per_P = round(d[["C_per_P"]], 2), omega = round(d[["omega"]], 3),
    exP = round(d[["excess_deficit_P"]], 3), exC = round(d[["excess_deficit_C"]], 3),
    lambda1 = round(d[["lambda1"]], 3), lambda2 = round(d[["lambda2"]], 4),
    beta = round(d[["beta"]], 3), gamma = round(d[["gamma"]], 3), k_pal = round(d[["k_pal"]], 2),
    ch_ratio_RMSE = round(c_$ratio_RMSE, 3), ch_ratio_r = round(c_$ratio_r, 3),
    ch_slope = round(c_$slope, 2), ch_CP_S180 = round(c_$CP_S180, 2)) }))
tab$dAICc <- round(tab$AICc - tab$AICc[tab$model == ING_BEST], 2)
tab <- tab[order(tab$AICc), ]
print(tab, row.names = FALSE)
cat(sprintf("Observed choice C:P at sucrose 180 = %.2f\n", mean(exp(lro[FS == 180]))))
write.csv(tab, "out/table_ext_models.csv", row.names = FALSE)
source("R/table_utils.R")
tS <- tab[, c("model", "k", "AICc", "dAICc", "LORO_pct", "LORO_no01_pct", "LOCO_pct", "C_per_P", "beta", "gamma", "k_pal",
              "ch_ratio_RMSE", "ch_ratio_r", "ch_slope", "ch_CP_S180")]
write_table(tS, "TableS2_extended_models",
  paste("Two non-linear extensions of the best ING model, estimated from the no-choice array: a Box-Cox",
        "regulatory scale (BC; beta = 1 is the proportional scale of the main model) and sucrose palatability",
        "(Pal; reward gamma (S/180)\\textsuperscript{k} per uL of food with sucrose concentration S).",
        "dAICc relative to the best model. LORO\\_no01: leave-one-rail-out error excluding the protein-free",
        "rail, on which BC models extrapolate to zero intake. ch\\_: out-of-sample prediction of the choice",
        "experiment (log C:P; slope of observed on predicted; mean predicted C:P at sucrose 180 g/L, observed 3.60)."),
  "tab:S2", "tables/supplementary", digits = 2)
saveRDS(list(E = E, S = S, CV = CV, CH = CH, tab = tab, dat = dat), "out/ext_fits.rds")
