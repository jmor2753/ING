# Trait-maximisation models compared with ING
source("R/ing_model.R"); suppressMessages(library(nutrigon)); library(parallel)
F0 <- readRDS("out/fits_nochoice.rds"); dat <- F0$dat; fits_ing <- F0$fits
TN <- readRDS("out/traits_nutrigon.rds")
data(lee_drosophila)
traits <- c("lifespan", "lifetimeegg", "dailyeggs")
T_MAX <- 2 * max(exp(dat$ybar))
coefW <- function(y) { cf <- nutrigon_coefs(TN$fits[[paste(y, "lm")]]); s <- sd(lee_drosophila[[y]])
  c(b1 = cf[[".x"]], b2 = cf[[".y"]], b3 = cf[["I(.x^2)"]], b4 = cf[["I(.y^2)"]], b5 = cf[[".x:.y"]]) / s }
rail_W <- function(b, l1, l2, comp, Fv) sapply(seq_len(nrow(comp)), function(i) {
  c1 <- comp[i, 1]; c2 <- comp[i, 2]; g <- sum(W_FOOD * comp[i, ]) / Fv[i]
  a2 <- -(b[["b3"]] * c1^2 + b[["b4"]] * c2^2 + b[["b5"]] * c1 * c2) + l2 * g^2
  a1 <- -(b[["b1"]] * c1 + b[["b2"]] * c2) + l1 * g
  cand <- c(1e-3, T_MAX); if (a2 > 0) cand <- c(cand, min(max(-a1 / (2 * a2), 1e-3), T_MAX))
  cand[which.min(a2 * cand^2 + a1 * cand)] })
nll_W <- function(b, dd, l1, l2, lt) { ts <- rail_W(b, l1, l2, dd$comp, dd$F)
  v <- exp(2 * lt) + dd$s2n; 0.5 * sum(log(2 * pi * v) + (dd$ybar - log(ts))^2 / v) }
fit_given_l2 <- function(b, dd, l2) { best <- NULL
  for (s in c(-3, -1, -.3, 0, .3, 1, 3) * l2 * 1e-2) { f <- optim(c(s, log(.3)), function(q) nll_W(b, dd, q[1], l2, q[2]))
    if (is.null(best) || f$value < best$value) best <- f }
  list(nll = best$value, l1 = best$par[1], lt = best$par[2]) }
grid_l2 <- 10^seq(-2, 6, by = 0.25)
res <- mclapply(traits, function(y) { b <- coefW(y)
  prof <- sapply(grid_l2, function(l2) fit_given_l2(b, dat, l2)$nll)
  k <- which.min(prof); best <- fit_given_l2(b, dat, grid_l2[k])
  err <- c(); for (r in unique(dat$tr$Ratio)) { tin <- which(dat$tr$Ratio != r); tout <- which(dat$tr$Ratio == r)
    dd <- subset_trt(dat, tin); pp <- sapply(grid_l2, function(l2) fit_given_l2(b, dd, l2)$nll)
    l2b <- grid_l2[which.min(pp)]; fb <- fit_given_l2(b, dd, l2b)
    pr <- rail_W(b, fb$l1, l2b, dat$comp[tout, , drop = FALSE], dat$F[tout]); err <- c(err, dat$ybar[tout] - log(pr)) }
  list(trait = y, profile = data.frame(log10_l2 = log10(grid_l2), dnll = prof - min(prof)),
       nll = best$nll, AICc = 2 * best$nll + 2 * 3 + 2 * 3 * 4 / (28 - 4), log10_l2 = log10(grid_l2[k]),
       dnll_trait_dominated = prof[1] - best$nll, LORO = sqrt(mean(err^2))) }, mc.cores = ING_CORES)
names(res) <- traits
cv <- read.csv("out/table_cv.csv")
tab <- do.call(rbind, lapply(res, function(r) data.frame(objective = paste("maximise", r$trait), k = 3,
  AICc = round(r$AICc, 2), log10_lambda2_best = r$log10_l2, dnll_if_trait_dominates = round(r$dnll_trait_dominated, 2),
  LORO_pct = round(100 * (exp(r$LORO) - 1), 1))))
for (n in ING_REP) tab <- rbind(tab, data.frame(objective = paste("ING", n),
  k = fits_ing[[n]]$k, AICc = round(fits_ing[[n]]$AICc, 2), log10_lambda2_best = NA, dnll_if_trait_dominates = NA,
  LORO_pct = cv$pct_rail[cv$model == n]))
print(tab, row.names = FALSE)
write.csv(tab, "out/table_adaptive_test.csv", row.names = FALSE)
saveRDS(res, "out/adaptive_test.rds")
