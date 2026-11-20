# Trait-maximisation test with thin-plate surfaces
source("R/ing_model.R"); suppressMessages(library(mgcv))
F0 <- readRDS("out/fits_nochoice.rds"); dat <- F0$dat; fits <- F0$fits
d <- read_lee("data/lee2008_nochoice.csv"); d <- d[complete.cases(d[, c("protein_eaten","carb_eaten","lifespan","lifetimeegg","dailyeggs")]), ]
d$P <- d$protein_eaten; d$C <- d$carb_eaten
tg <- seq(5, 400, by = 1); g_i <- as.vector(dat$comp %*% W_FOOD) / dat$F
res <- list()
for (y in c("lifespan","lifetimeegg","dailyeggs")) for (s in c("quadratic","thinplate")) {
  m <- if (s == "quadratic") lm(d[[y]] ~ P + C + I(P^2) + I(C^2) + I(P*C), d) else readRDS(paste0("out/gam_", y, ".rds"))
  Wm <- t(sapply(seq_len(nrow(dat$comp)), function(i) predict(m, newdata = data.frame(P = tg * dat$comp[i, 1], C = tg * dat$comp[i, 2]))))
  Wm <- Wm / sd(d[[y]])
  fwd <- function(l1, l2) sapply(seq_len(nrow(Wm)), function(i) { V <- g_i[i] * tg; tg[which.max(Wm[i, ] - l1 * V - l2 * V^2)] })
  nll_at <- function(l1, l2, lt) { ts <- fwd(l1, l2); v <- exp(2 * lt) + dat$s2n; 0.5 * sum(log(2*pi*v) + (dat$ybar - log(ts))^2 / v) }
  prof <- sapply(seq(-2, 5, by = .5), function(lg) { best <- Inf
    for (l1 in c(-3,-1,-.3,0,.3,1,3) * 10^lg) { f <- optimize(function(lt) nll_at(l1, 10^lg, lt), c(-6, 1)); best <- min(best, f$objective) }; best })
  res[[paste(y, s)]] <- data.frame(trait = y, surface = s, best_nll = round(min(prof), 2),
    AICc = round(2 * min(prof) + 6 + 24 / 24, 2), log10_l2_at_best = seq(-2, 5, by = .5)[which.min(prof)],
    nll_W_dominant = round(prof[1], 2))
}
tab <- do.call(rbind, res); rownames(tab) <- NULL; print(tab)
cat("ING best AICc:", round(fits[[ING_BEST]]$AICc, 2), "\n")
write.csv(tab, "out/table_trait_forward_gam.csv", row.names = FALSE)
