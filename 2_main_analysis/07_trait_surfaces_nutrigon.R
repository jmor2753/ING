# Performance surfaces, optima and angular separations (nutrigon)
source("R/ing_model.R"); suppressMessages(library(nutrigon)); library(parallel)
data(lee_drosophila)
L <- lee_drosophila
traits <- c(lifespan = "Lifespan", lifetimeegg = "Lifetime egg production",
            dailyeggs = "Egg production rate")
XC <- "protein_eaten"; YC <- "carb_eaten"
theta_to_ratio <- function(th) sprintf("1:%.2f", tan(th * pi / 180))
segdist <- function(p, a, b) { ab <- b - a; if (sum(ab^2) == 0) return(sqrt(sum((p - a)^2)))
  t <- max(0, min(1, sum((p - a) * ab) / sum(ab^2))); sqrt(sum((p - a - t * ab)^2)) }
edge_dist <- function(p, hx, hy) { n <- length(hx)
  min(sapply(seq_len(n), function(i) { j <- if (i == n) 1 else i + 1; segdist(p, c(hx[i], hy[i]), c(hx[j], hy[j])) })) }

fit_trait <- function(y, model_type = "lm", ...) {
  set.seed(2025)
  nutrigon_fit(L, XC, YC, y, design_id_col = "treatment", model_type = model_type,
               n_cloud = 4000, n_reps = 10, ...)
}
fits <- list(); rows <- list()
specs <- list(list(y = "lifespan", mt = "lm"), list(y = "lifetimeegg", mt = "lm"), list(y = "dailyeggs", mt = "lm"),
              list(y = "lifetimeegg", mt = "glm", fam = quasipoisson()), list(y = "dailyeggs", mt = "glm", fam = quasipoisson()))
for (sp in specs) {
  f <- if (sp$mt == "glm") fit_trait(sp$y, "glm", family = sp$fam) else fit_trait(sp$y)
  key <- paste(sp$y, sp$mt); fits[[key]] <- f
  cf <- nutrigon_coefs(f)
  H  <- matrix(c(2 * cf[["I(.x^2)"]], cf[[".x:.y"]], cf[[".x:.y"]], 2 * cf[["I(.y^2)"]]), 2)
  ev <- eigen(H, symmetric = TRUE)$values
  shape <- if (all(ev < 0)) "peak" else if (all(ev > 0)) "valley" else "saddle"
  o  <- f$optimum; th <- o$theta
  rg <- f$region; rth <- atan2(rg[[YC]], rg[[XC]]) * 180 / pi
  rows[[key]] <- data.frame(trait = sp$y, model = sp$mt, surface_shape = shape,
    closed_form = !o$used_fallback, P_opt = round(o[[XC]], 1), C_opt = round(o[[YC]], 1),
    theta_deg = round(th, 2), ratio = theta_to_ratio(th),
    region_theta_lo = round(quantile(rth, .025), 1), region_theta_hi = round(quantile(rth, .975), 1),
    region_ratio = paste0(theta_to_ratio(quantile(rth, .025)), " to ", theta_to_ratio(quantile(rth, .975))),
    dist_to_hull_edge = round(edge_dist(c(o[[XC]], o[[YC]]), f$hull$x, f$hull$y), 2),
    elongation = round(nutrigon_elongation(f), 2))
}
tab_opt <- do.call(rbind, rows); rownames(tab_opt) <- NULL
print(tab_opt)
write.csv(tab_opt, "out/table_trait_nutrigon.csv", row.names = FALSE)

pairs <- list(c("lifespan", "lifetimeegg"), c("lifespan", "dailyeggs"), c("lifetimeegg", "dailyeggs"))
cmp <- mclapply(pairs, function(pr) nutrigon_compare(
  L, XC, YC, pr[1], L, XC, YC, pr[2], design_id_col1 = "treatment", design_id_col2 = "treatment",
  paired = TRUE, B = 500, n_cloud = 2000, n_reps_cloud = 5, seed = 2025), mc.cores = ING_CORES)
names(cmp) <- sapply(pairs, paste, collapse = " vs ")
GRID_N <- as.integer(Sys.getenv("ING_GRID_N", "100"))
chk <- lapply(seq_along(pairs), function(i) {
  pr <- pairs[[i]]; f1 <- fits[[paste(pr[1], "lm")]]; f2 <- fits[[paste(pr[2], "lm")]]
  gf <- sprintf("out/validation_grid_%s_%s.rds", pr[1], pr[2])
  if (file.exists(gf) && Sys.getenv("ING_REBUILD_GRID") != "1") g <- readRDS(gf) else {
    g <- nutrigon_build_grid(n_grid_points = GRID_N, n_sim_per_point = 20, B_boot = 50,
           design = L, x_col = XC, y_col = YC, design_id_col = "treatment",
           snr = c(nutrigon_snr(f1), nutrigon_snr(f2)), paired = TRUE, seed = 99, verbose = FALSE)
    saveRDS(g, gf) }
  k <- tryCatch(nutrigon_check(cmp[[i]], elong1 = nutrigon_elongation(f1), elong2 = nutrigon_elongation(f2),
                               validation_grid = g, k = 10, min_angle = 5), error = function(e) NULL)
  if (!is.null(k)) k$grid_support <- attr(g, "angular_support")
  k })
tab_cmp <- do.call(rbind, lapply(seq_along(cmp), function(i) { r <- cmp[[i]]; k <- chk[[i]]
  data.frame(comparison = names(cmp)[i], theta1 = round(r$theta1, 2), theta2 = round(r$theta2, 2),
    separation_deg = round(r$theta_hat, 2), lo95 = round(r$ci[1], 2), hi95 = round(r$ci[2], 2),
    n_boot_failed = sum(is.na(r$boot_diffs)),
    check_flag = if (is.null(k)) NA else k$flag,
    check_coverage = if (is.null(k)) NA else round(k$estimated_coverage, 3),
    coverage_lo = if (is.null(k)) NA else round(k$coverage_ci[1], 3),
    coverage_hi = if (is.null(k)) NA else round(k$coverage_ci[2], 3)) }))
print(tab_cmp, row.names = FALSE)
write.csv(tab_cmp, "out/table_trait_separation.csv", row.names = FALSE)

F0 <- readRDS("out/fits_nochoice.rds"); ing_best <- F0$fits[[ING_BEST]]
draw <- read_lee("data/lee2008_nochoice.csv")
th_ing <- function(f) atan2(f$p$mu[2], f$p$mu[1]) * 180 / pi
opt_theta <- function(dd, y) { set.seed(1)
  f <- nutrigon_fit(dd, XC, YC, y, model_type = "lm", n_cloud = 1000, n_reps = 3); f$optimum$theta }
idx_by <- split(seq_len(nrow(draw)), paste(draw$Ratio, draw$Food))
B <- 200
boot <- do.call(rbind, mclapply(seq_len(B), function(b) { set.seed(5000 + b)
  ii <- unlist(lapply(idx_by, function(ix) ix[sample.int(length(ix), replace = TRUE)]))
  db <- draw[ii, ]
  fi <- fit_ing(prep_intake(db), ing_best$spec, n_start = 8, seed = b, start = ing_best$par)
  dt <- db[complete.cases(db[, c(XC, YC, names(traits))]), ]
  c(ing = th_ing(fi), sapply(names(traits), function(y) opt_theta(dt, y))) }, mc.cores = ING_CORES))
est <- c(ing = th_ing(ing_best), sapply(names(traits), function(y) fits[[paste(y, "lm")]]$optimum$theta))
tab_ing <- do.call(rbind, lapply(names(traits), function(y) { dlt <- boot[, y] - boot[, "ing"]
  data.frame(trait = y, theta_ING = round(est["ing"], 2), theta_trait = round(est[y], 2),
    separation_deg = round(est[y] - est["ing"], 2), lo95 = round(quantile(dlt, .025), 2), hi95 = round(quantile(dlt, .975), 2),
    ING_ratio = theta_to_ratio(est["ing"]), trait_ratio = theta_to_ratio(est[y])) }))
chd <- load_choice(6); CP180 <- mean(chd$C_per_P[chd$sucrose_gL == 180])
tab_ing <- rbind(tab_ing, data.frame(trait = "choice (sucrose 180, observed)", theta_ING = round(est["ing"], 2),
  theta_trait = round(atan(CP180) * 180 / pi, 2), separation_deg = round(atan(CP180) * 180 / pi - est["ing"], 2),
  lo95 = NA, hi95 = NA, ING_ratio = theta_to_ratio(est["ing"]), trait_ratio = sprintf("1:%.2f", CP180)))
cat("\nING theta bootstrap 95% CI:", round(quantile(boot[, "ing"], c(.025, .975)), 2), "\n")
print(tab_ing, row.names = FALSE)
write.csv(tab_ing, "out/table_ing_vs_traits.csv", row.names = FALSE)
saveRDS(list(fits = fits, tab_opt = tab_opt, cmp = cmp, chk = chk, tab_cmp = tab_cmp,
             boot = boot, est = est, tab_ing = tab_ing), "out/traits_nutrigon.rds")
