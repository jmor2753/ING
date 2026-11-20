# Thin-plate performance surfaces
source("R/ing_model.R"); library(mgcv)
d <- read_lee("data/lee2008_nochoice.csv")
d <- d[complete.cases(d[, c("protein_eaten","carb_eaten","lifespan","lifetimeegg","dailyeggs")]), ]
d$P <- d$protein_eaten; d$C <- d$carb_eaten
traits <- c("lifespan","lifetimeegg","dailyeggs")
gx <- seq(0, 180, length.out = 181); gy <- seq(0, 300, length.out = 181); G <- expand.grid(P = gx, C = gy)
far <- exclude.too.far(G$P, G$C, d$P, d$C, dist = 0.05)
num_hess <- function(f, x, h = c(1, 1)) {
  H <- matrix(0, 2, 2)
  for (i in 1:2) for (j in 1:2) { ei <- diag(2)[i, ] * h[i]; ej <- diag(2)[j, ] * h[j]
    H[i, j] <- (f(x + ei + ej) - f(x + ei - ej) - f(x - ei + ej) + f(x - ei - ej)) / (4 * h[i] * h[j]) }
  H }
res <- list()
for (y in traits) {
  q <- lm(d[[y]] ~ P + C + I(P^2) + I(C^2) + I(P*C), d)
  g <- gam(as.formula(paste(y, "~ s(P, C, k = 40)")), data = d, method = "REML")
  for (nm in c("quadratic", "thinplate")) {
    m <- if (nm == "quadratic") q else g
    pr <- predict(m, newdata = G); pr[far] <- NA
    i <- which.max(pr); x0 <- c(G$P[i], G$C[i])
    f <- function(x) predict(m, newdata = data.frame(P = x[1], C = x[2]))
    o <- optim(x0, function(x) -f(x), method = "L-BFGS-B", lower = pmax(x0 - 15, 0), upper = x0 + 15)
    xs <- o$par; H <- num_hess(f, xs, c(2, 2)); ev <- eigen(-H, symmetric = TRUE)$values
    on_edge <- any(abs(xs - (x0 - 15)) < 1e-3 | abs(xs - (x0 + 15)) < 1e-3)
    nb <- expand.grid(dP = c(-6, 0, 6), dC = c(-10, 0, 10))
    interior <- !any(exclude.too.far(xs[1] + nb$dP, xs[2] + nb$dC, d$P, d$C, dist = 0.05))
    res[[paste(y, nm)]] <- data.frame(trait = y, surface = nm, P = round(xs[1], 1), C = round(xs[2], 1),
      C_per_P = round(xs[2] / xs[1], 2), negH_eig1 = signif(ev[1], 3), negH_eig2 = signif(ev[2], 3),
      local_max = all(ev > 0) && !on_edge, inside_support = interior,
      dev_expl = round(if (nm == "quadratic") summary(q)$r.squared else summary(g)$dev.expl, 3))
  }
  saveRDS(g, paste0("out/gam_", y, ".rds"))
}
tab <- do.call(rbind, res); print(tab, row.names = FALSE)
write.csv(tab, "out/table_trait_audit.csv", row.names = FALSE)
