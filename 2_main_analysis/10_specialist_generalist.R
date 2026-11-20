# ING for the specialist-generalist caterpillar pairs
source("R/ing_model.R"); source("R/table_utils.R")
set.seed(10)
B_BOOT <- as.integer(Sys.getenv("ING_BOOT_LARVAE", "500"))

read_larvae <- function(file) {
  d <- read.csv(file, stringsAsFactors = FALSE)
  d <- d[d$species != "", ]
  ch <- d$experiment == "Choice"
  d$pair <- ifelse(ch, sub("p14:c[0-9]+$", "p14:c28", d$treat), NA)
  d$diet <- ifelse(ch, NA, d$treat)
  d$id <- paste(d$species, d$experiment, ifelse(ch, d$pair, d$diet), d$replicate)
  out <- do.call(rbind, lapply(split(d, d$id), function(z) {
    z <- z[order(z$day2), ]; ok <- which(!is.na(z$protein) & !is.na(z$carb))
    if (!length(ok) || z$day2[max(ok)] < max(z$day2)) return(NULL)
    k <- max(ok)
    data.frame(species = z$species[1], experiment = z$experiment[1], diet = z$diet[1],
               pair = z$pair[1], larva = z$replicate[1], days = z$day2[k],
               P = z$protein[k], C = z$carb[k])
  }))
  out
}
LV <- rbind(read_larvae("data/raw/lee_heliothis_larvae_raw.csv"),
            read_larvae("data/raw/lee_heliothis_hybrid_raw.csv"),
            read_larvae("data/raw/lee_spodoptera_larvae_raw.csv"))
rownames(LV) <- NULL
SPP <- data.frame(
  species = c("Heliothis virescens", "Heliothis subflexa", "Heliothis subflexa x virescens", "Spodoptera littoralis", "Spodoptera exempta"),
  genus = c("Heliothis", "Heliothis", "Heliothis", "Spodoptera", "Spodoptera"),
  diet_breadth = c("generalist", "specialist", "hybrid", "generalist", "specialist"),
  short = c("H. virescens", "H. subflexa", "H. subflexa x virescens", "S. littoralis", "S. exempta"), stringsAsFactors = FALSE)
LV$diet_breadth <- SPP$diet_breadth[match(LV$species, SPP$species)]
write.csv(LV, "data/lee_larvae_endpoints.csv", row.names = FALSE)

pc_of <- function(lbl) { v <- as.numeric(regmatches(lbl, gregexpr("[0-9]+", lbl))[[1]]); v[1:2] }
comp_of <- function(lbl) { v <- pc_of(lbl); v / sum(v) }

lspec <- function(aniso, asym = character(0))
  spec(aniso = aniso, asym = asym, l1 = FALSE, l2 = FALSE)
MODELS <- list(Iso = lspec(FALSE), `Iso+AsymP` = lspec(FALSE, "P"), `Iso+AsymC` = lspec(FALSE, "C"),
               `Iso+AsymPC` = lspec(FALSE, c("P", "C")), Aniso = lspec(TRUE))
prep_sp <- function(dd) {
  diets <- sort(unique(dd$diet)); trt <- match(dd$diet, diets)
  list(lt = log(dd$P + dd$C), trt = trt, diets = diets, n = nrow(dd),
       comp = t(sapply(diets, comp_of)), Fv = rep(1, length(diets)))
}
pred_t <- function(par, s, dat) rail_opt(unpack(par, s), dat$comp, dat$Fv)
nll_ind <- function(par, s, dat) {
  p <- unpack(par, s); if (any(abs(p$eta) > 4)) return(1e10)
  ts <- rail_opt(p, dat$comp, dat$Fv); if (any(!is.finite(ts) | ts <= 1e-8)) return(1e10)
  rss <- sum((dat$lt - log(ts)[dat$trt])^2); dat$n / 2 * (log(2 * pi * rss / dat$n) + 1)
}
fit_l <- function(dat, s, n_start = 25, start = NULL) {
  np <- n_par(s); best <- NULL
  tm <- as.vector(exp(tapply(dat$lt, dat$trt, mean))); mu0 <- log(colSums(dat$comp * tm) / nrow(dat$comp) * 1.2)
  for (k in seq_len(n_start)) {
    p0 <- if (!is.null(start) && k == 1) start else c(mu0 + rnorm(2, 0, .3), rnorm(np - 2, 0, .6))
    f <- try(optim(p0, nll_ind, s = s, dat = dat, control = list(maxit = 5000, reltol = 1e-12)), silent = TRUE)
    if (inherits(f, "try-error") || f$value >= 1e10) next
    f2 <- try(optim(f$par, nll_ind, s = s, dat = dat, method = "BFGS", control = list(maxit = 2000)), silent = TRUE)
    if (!inherits(f2, "try-error") && f2$value <= f$value) f <- f2
    if (is.null(best) || f$value < best$value) best <- f
  }
  k <- np + 1; p <- unpack(best$par, s)
  list(spec = s, par = best$par, nll = best$value, k = k, p = p,
       AICc = 2 * best$value + 2 * k + 2 * k * (k + 1) / (dat$n - k - 1),
       pred = rail_opt(p, dat$comp, dat$Fv))
}
nll_sat <- function(dat) {
  rss <- sum((dat$lt - ave(dat$lt, dat$trt))^2); dat$n / 2 * (log(2 * pi * rss / dat$n) + 1) }

choice_pred <- function(p, pair) {
  f <- strsplit(pair, " vs ")[[1]]; dA <- comp_of(f[1]); dB <- comp_of(f[2])
  ang <- function(x) atan2(x[2], x[1]); a <- sort(c(ang(dA), ang(dB)))
  if (ang(p$mu) >= a[1] && ang(p$mu) <= a[2]) return(p$mu)
  cand <- rbind(rail_opt(p, rbind(dA), 1) * dA, rail_opt(p, rbind(dB), 1) * dB)
  cand[which.min(loss_R(cand, p)), ]
}
summ_p <- function(p) { e <- eigen(p$Qz, symmetric = TRUE); v1 <- e$vectors[, 1]; if (v1[1] < 0) v1 <- -v1
  c(muP = p$mu[[1]], muC = p$mu[[2]], CP = p$mu[[2]] / p$mu[[1]], omega = e$values[1],
    phi = (atan2(v1[2], v1[1]) * 180 / pi) %% 180) }

lin_pred <- function(par, comp) { a <- c(cos(par[1]), sin(par[1])); exp(par[2]) / drop(comp %*% a) }
nll_lin <- function(par, dat) {
  if (par[1] <= 0.01 || par[1] >= pi / 2 - 0.01) return(1e10)
  rss <- sum((dat$lt - log(lin_pred(par, dat$comp))[dat$trt])^2); dat$n / 2 * (log(2 * pi * rss / dat$n) + 1) }
fit_lin <- function(dat, start = NULL) {
  best <- NULL
  for (p0 in c(list(start), lapply(seq(.2, 1.3, length.out = 6), function(ps) c(ps, log(100) )))) {
    if (is.null(p0)) next
    f <- optim(p0, nll_lin, dat = dat, control = list(reltol = 1e-12, maxit = 4000))
    if (is.null(best) || f$value < best$value) best <- f }
  k <- 3; list(par = best$par, nll = best$value, k = k, w = 1 / tan(best$par[1]),
               level = exp(best$par[2]) / sin(best$par[1]),
               AICc = 2 * best$value + 2 * k + 2 * k * (k + 1) / (dat$n - k - 1),
               pred = lin_pred(best$par, dat$comp)) }
lin_dev <- function(par, x) drop(x %*% c(cos(par[1]), sin(par[1]))) / exp(par[2]) - 1
prof_dir <- function(dat, angles = seq(10, 80, 5), n_start = 12) {
  s <- MODELS$Aniso
  sapply(angles, function(a) { th <- a * pi / 180
    f <- function(q) nll_ind(c(q[1] + log(c(cos(th), sin(th))), q[2:3]), s, dat)
    min(sapply(seq_len(n_start), function(k) optim(c(log(runif(1, 50, 500)), rnorm(2, 0, 1)), f)$value)) }) }

loro_err <- function(dat, fitfun) {
  e <- sapply(seq_along(dat$diets), function(i) {
    keep <- dat$trt != i
    d4 <- list(lt = dat$lt[keep], trt = match(dat$trt[keep], setdiff(seq_along(dat$diets), i)),
               n = sum(keep), comp = dat$comp[-i, , drop = FALSE], Fv = dat$Fv[-i])
    log(fitfun(d4, dat$comp[i, , drop = FALSE])) - mean(dat$lt[dat$trt == i]) })
  100 * (exp(sqrt(mean(e^2))) - 1) }
res <- list()
for (sp in SPP$species) {
  cat("\n====", sp, "====\n")
  nc <- LV[LV$species == sp & LV$experiment == "Nochoice", ]; dat <- prep_sp(nc)
  fits <- lapply(MODELS, fit_l, dat = dat); lin <- fit_lin(dat)
  sat <- nll_sat(dat); k_sat <- length(dat$diets) + 1
  aicc <- function(nll, k) 2 * nll + 2 * k + 2 * k * (k + 1) / (dat$n - k - 1)
  tab <- data.frame(model = c(names(fits), "Linear", "Saturated (diet means)"),
                    k = c(sapply(fits, `[[`, "k"), lin$k, k_sat),
                    nll = c(sapply(fits, `[[`, "nll"), lin$nll, sat))
  tab$AICc <- aicc(tab$nll, tab$k); tab$dAICc <- tab$AICc - min(tab$AICc)
  tab$LOF_df <- k_sat - tab$k; tab$LOF_LR <- 2 * (tab$nll - sat)
  tab$LOF_p <- ifelse(tab$LOF_df > 0, pchisq(tab$LOF_LR, tab$LOF_df, lower.tail = FALSE), NA)
  tab$LORO_pct <- c(sapply(names(MODELS), function(m) loro_err(dat, function(d4, ci) {
                      f <- fit_l(d4, MODELS[[m]], n_start = 15, start = fits[[m]]$par); rail_opt(f$p, ci, 1) })),
                    loro_err(dat, function(d4, ci) lin_pred(fit_lin(d4, lin$par)$par, ci)), NA)
  rownames(tab) <- NULL; print(tab, digits = 3)
  prof <- prof_dir(dat); prof <- prof - min(prof)
  ch <- LV[LV$species == sp & LV$experiment == "Choice", ]; pair <- ch$pair[1]
  xo <- cbind(ch$P, ch$C)
  dev_obs <- mean(log1p(lin_dev(lin$par, xo)))
  bt <- t(sapply(seq_len(B_BOOT), function(b) {
    ii <- unlist(lapply(split(seq_len(nrow(nc)), nc$diet), function(v) v[sample.int(length(v), replace = TRUE)]))
    f <- fit_lin(prep_sp(nc[ii, ]), lin$par); jj <- sample.int(nrow(ch), replace = TRUE)
    c(w = f$w, level = f$level, psi = f$par[1], lk = f$par[2],
      dev = mean(log1p(lin_dev(f$par, xo[jj, , drop = FALSE]))),
      obsP = mean(log(ch$P[jj])), obsC = mean(log(ch$C[jj]))) }))
  res[[sp]] <- list(dat = dat, nc = nc, ch = ch, pair = pair, fits = fits, lin = lin, table = tab,
                    profile = data.frame(angle = seq(10, 80, 5), dnll = prof), boot = bt,
                    dev_obs = dev_obs, obs_choice = colMeans(log(xo)))
  cat(sprintf("Linear: w = %.2f [%.2f, %.2f]; choice deviation %.1f%% [%.1f, %.1f]; profile range %.3f\n",
      lin$w, quantile(bt[, "w"], .025), quantile(bt[, "w"], .975), 100 * expm1(dev_obs),
      100 * expm1(quantile(bt[, "dev"], .025)), 100 * expm1(quantile(bt[, "dev"], .975)), max(prof)))
}

tc <- list()
for (f in c("data/raw/lee_heliothis_larvae_raw.csv", "data/raw/lee_heliothis_hybrid_raw.csv", "data/raw/lee_spodoptera_larvae_raw.csv")) {
  d <- read.csv(f, stringsAsFactors = FALSE); d <- d[d$species != "" & d$experiment == "Nochoice", ]
  for (sp in unique(d$species)) for (dy in 1:max(d$day2[d$species == sp])) {
    z <- d[d$species == sp & d$day2 == dy & !is.na(d$protein), ]
    if (length(unique(z$treat)) < 5) next
    dz <- prep_sp(data.frame(diet = z$treat, P = z$protein, C = z$carb)); fz <- fit_lin(dz)
    bz <- replicate(200, { ii <- unlist(lapply(split(seq_len(nrow(z)), z$treat), function(v) v[sample.int(length(v), replace = TRUE)]))
                           fit_lin(prep_sp(data.frame(diet = z$treat[ii], P = z$protein[ii], C = z$carb[ii])), fz$par)$w })
    tc[[length(tc) + 1]] <- data.frame(species = sp, day = dy, w = fz$w, lo = quantile(bz, .025), hi = quantile(bz, .975),
                                       level = fz$level)
  }
}
tc <- do.call(rbind, tc); rownames(tc) <- NULL
saveRDS(list(res = res, SPP = SPP, LV = LV, timecourse = tc), "out/larvae_fits.rds")

contr <- do.call(rbind, lapply(c("Heliothis", "Spodoptera"), function(g) {
  sg <- SPP$species[SPP$genus == g & SPP$diet_breadth == "specialist"]; gg <- SPP$species[SPP$genus == g & SPP$diet_breadth == "generalist"]
  r <- res[[sg]]$boot[, "w"] / res[[gg]]$boot[, "w"]
  data.frame(genus = g, contrast = "specialist / generalist", ratio = res[[sg]]$lin$w / res[[gg]]$lin$w,
             lo = quantile(r, .025), hi = quantile(r, .975)) }))
hy <- "Heliothis subflexa x virescens"; pv <- "Heliothis virescens"; ps <- "Heliothis subflexa"
hc <- function(num, den, lab) { r <- res[[num]]$boot[, "w"] / res[[den]]$boot[, "w"]
  data.frame(genus = "Heliothis", contrast = lab, ratio = res[[num]]$lin$w / res[[den]]$lin$w,
             lo = quantile(r, .025), hi = quantile(r, .975)) }
pos <- with(list(h = log(res[[hy]]$boot[, "w"]), v = log(res[[pv]]$boot[, "w"]), s = log(res[[ps]]$boot[, "w"])), (h - v) / (s - v))
contr <- rbind(contr, hc(hy, pv, "hybrid / H. virescens"), hc(hy, ps, "hybrid / H. subflexa"),
  data.frame(genus = "Heliothis", contrast = "hybrid position (0 = H. virescens, 1 = H. subflexa; log w)",
             ratio = (log(res[[hy]]$lin$w) - log(res[[pv]]$lin$w)) / (log(res[[ps]]$lin$w) - log(res[[pv]]$lin$w)),
             lo = quantile(pos, .025), hi = quantile(pos, .975)))
rownames(contr) <- NULL
print(contr); saveRDS(contr, "out/larvae_contrasts.rds")

cat("\nDone. Tables 6 and S14: 4_main_tables/table6_larvae.R\n")
