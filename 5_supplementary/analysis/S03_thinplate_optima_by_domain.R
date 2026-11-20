# Thin-plate optima by domain
source("R/ing_model.R"); suppressMessages(library(mgcv))
d <- read_lee("data/lee2008_nochoice.csv"); d <- d[complete.cases(d[, c("protein_eaten","carb_eaten","lifespan","lifetimeegg","dailyeggs")]), ]
d$P <- d$protein_eaten; d$C <- d$carb_eaten
pip <- function(px, py, h) { n <- nrow(h); ins <- rep(FALSE, length(px)); j <- n
  for (i in 1:n) { cr <- ((h[i,2] > py) != (h[j,2] > py)) & (px < (h[j,1]-h[i,1])*(py-h[i,2])/(h[j,2]-h[i,2]) + h[i,1]); ins <- xor(ins, cr); j <- i }; ins }
hull_of <- function(X) X[chull(X), ]
peel <- function(X, frac = .95) {
  keep <- rep(TRUE, nrow(X)); while (mean(keep) > frac) { idx <- which(keep); h <- idx[chull(X[idx, ])]; keep[h] <- FALSE }; X[keep, ] }
tm <- aggregate(cbind(P, C) ~ Ratio + Food, d, mean)
doms <- list(`treatment means` = hull_of(as.matrix(tm[, c("P","C")])),
             `individuals (95% peeled)` = hull_of(peel(as.matrix(d[, c("P","C")]))),
             `individuals (all)` = hull_of(as.matrix(d[, c("P","C")])))
G <- expand.grid(P = seq(0, 185, by = .5), C = seq(0, 300, by = .5))
out <- list()
for (y in c("lifespan","lifetimeegg","dailyeggs")) {
  q <- lm(d[[y]] ~ P + C + I(P^2) + I(C^2) + I(P*C), d); g <- readRDS(paste0("out/gam_", y, ".rds"))
  for (nm in names(doms)) { ins <- pip(G$P, G$C, doms[[nm]])
    for (s in c("quadratic","thinplate")) { m <- if (s == "quadratic") q else g
      pr <- predict(m, newdata = G[ins, ]); k <- which.max(pr); x <- unlist(G[ins, ][k, ])
      out[[length(out)+1]] <- data.frame(trait = y, domain = nm, surface = s, P = x[1], C = x[2], ratio = sprintf("1:%.1f", x[2]/x[1])) } } }
o <- do.call(rbind, out); rownames(o) <- NULL
print(reshape(o[, c("trait","domain","surface","ratio")], idvar = c("trait","domain"), timevar = "surface", direction = "wide"), row.names = FALSE)
write.csv(o, "out/table_trait_domain.csv", row.names = FALSE)
saveRDS(doms, "out/trait_domains.rds")
