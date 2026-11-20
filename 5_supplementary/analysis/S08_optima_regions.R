# Peak regions of the performance surfaces (OptimaRegion)
source("R/ing_model.R"); source("R/fig_theme.R"); source("R/table_utils.R")
if (!requireNamespace("OptimaRegion", quietly = TRUE)) stop("Install the CRAN package OptimaRegion")
dir.create("out/optima", recursive = TRUE, showWarnings = FALSE)
d <- read_lee("data/lee2008_nochoice.csv")
d <- d[complete.cases(d[, c("protein_eaten", "carb_eaten", "lifespan", "lifetimeegg", "dailyeggs")]), ]
d$P <- d$protein_eaten; d$C <- d$carb_eaten
X <- cbind(d$P, d$C); H <- X[rev(chull(X)), ]
F0 <- readRDS("out/fits_nochoice.rds"); best <- F0$fits[[ING_BEST]]
CB <- readRDS("out/cv_boot.rds"); TB <- CB$boots[[ING_BEST]][, c("muP", "muC")]
target <- best$p$mu
ALPHA <- .05; NOSIM <- 2000

r_max <- max(X[, 1] / X[, 2])
tri <- NULL
for (a in seq(max(X[, 2]), 3 * max(X[, 2]), length.out = 400)) {
  ok <- function(s) all(X[, 2] <= a + (s - a) * X[, 1] / (s * r_max) + 1e-9)
  ss <- seq(max(X[, 1]) / r_max, 5 * max(X[, 1]) / r_max, length.out = 400); s <- ss[which(sapply(ss, ok))[1]]
  if (!is.na(s) && (is.null(tri) || a * s < tri$a * tri$s)) tri <- list(a = a, s = s) }
V1 <- c(0.01, tri$a); V2 <- c(tri$s * r_max, tri$s)
inhull <- function(p) { n <- nrow(H); all(sapply(seq_len(n), function(i) { a <- H[i, ]; b <- H[if (i == n) 1 else i + 1, ]
  (b[1] - a[1]) * (p[2] - a[2]) - (b[2] - a[2]) * (p[1] - a[1]) >= -1e-9 })) }
traits <- c(lifespan = "Lifespan", lifetimeegg = "Lifetime eggs", dailyeggs = "Egg-production rate")
OR <- lapply(names(traits), function(y) { set.seed(1)
  invisible(capture.output(o <- OptimaRegion::OptRegionQuad(X, d[[y]], nosim = 500, alpha = ALPHA, LB = c(0, 0), UB = apply(X, 2, max),
    triangularRegion = TRUE, vertex1 = V1, vertex2 = V2, xlab = "Protein eaten (ug/day)", ylab = "Carbohydrate eaten (ug/day)",
    outputPDFFile = file.path("out/optima", paste0("OptimaRegion_", y, ".pdf")))))
  list(mean = o$meanPoint, xin = o$xin, in_hull = mean(apply(o$xin, 1, inhull))) }); names(OR) <- names(traits)

qf <- function(b, P, C) b[1] + b[2] * P + b[3] * C + b[4] * P^2 + b[5] * C^2 + b[6] * P * C
max_on_hull <- function(b) {
  cand <- H; A <- matrix(c(2 * b[4], b[6], b[6], 2 * b[5]), 2)
  st <- tryCatch(solve(A, -b[2:3]), error = function(e) NULL)
  if (!is.null(st) && inhull(st)) cand <- rbind(cand, st)
  n <- nrow(H)
  for (i in seq_len(n)) { a <- H[i, ]; e <- H[if (i == n) 1 else i + 1, ] - a
    c2 <- b[4] * e[1]^2 + b[5] * e[2]^2 + b[6] * e[1] * e[2]
    c1 <- b[2] * e[1] + b[3] * e[2] + 2 * b[4] * a[1] * e[1] + 2 * b[5] * a[2] * e[2] + b[6] * (a[1] * e[2] + a[2] * e[1])
    if (c2 < 0) { s <- -c1 / (2 * c2); if (s > 0 && s < 1) cand <- rbind(cand, a + s * e) } }
  cand[which.max(qf(b, cand[, 1], cand[, 2])), ] }
fit_surface <- function(y, family) {
  f <- as.formula(paste(y, "~ P + C + I(P^2) + I(C^2) + I(P * C)"))
  m <- if (family == "lm") lm(f, d) else glm(f, quasipoisson, d)
  list(b = coef(m), V = vcov(m)) }
cr_hull <- function(y, family) {
  s <- fit_surface(y, family); set.seed(2)
  Bs <- MASS::mvrnorm(NOSIM, s$b, s$V)
  md <- mahalanobis(Bs, s$b, s$V); keep <- md <= qchisq(1 - ALPHA, length(s$b))
  O <- t(apply(Bs[keep, , drop = FALSE], 1, max_on_hull))
  list(opt = max_on_hull(s$b), draws = O, region = O[chull(O), , drop = FALSE], n = sum(keep)) }
surf <- list(lifespan = "lm", lifetimeegg = "lm", dailyeggs = "lm", lifetimeegg = "quasi-Poisson", dailyeggs = "quasi-Poisson")
CRH <- Map(function(y, fam) c(list(trait = y, family = fam), cr_hull(y, fam)), names(surf), surf)
names(CRH) <- paste(names(surf), unlist(surf))

in_poly <- function(p, Pg) { if (nrow(Pg) < 3) return(FALSE); n <- nrow(Pg); s <- sapply(seq_len(n), function(i) { a <- Pg[i, ]; b <- Pg[if (i == n) 1 else i + 1, ]
  (b[1] - a[1]) * (p[2] - a[2]) - (b[2] - a[2]) * (p[1] - a[1]) }); all(s >= -1e-9) || all(s <= 1e-9) }
q3 <- function(v) quantile(v, c(.5, .025, .975))
on_edge <- function(p) { n <- nrow(H); min(sapply(seq_len(n), function(i) { a <- H[i, ]; e <- H[if (i == n) 1 else i + 1, ] - a
  s <- max(0, min(1, sum((p - a) * e) / sum(e^2))); sqrt(sum((p - a - s * e)^2)) })) < 1e-6 }
set.seed(3); nd <- 2000
dist_draws <- function(A, B) { i <- sample(nrow(A), nd, TRUE); j <- sample(nrow(B), nd, TRUE); sqrt(rowSums((A[i, ] - B[j, ])^2)) }
theta <- function(M) atan2(M[, 2], M[, 1]) * 180 / pi
tab12 <- do.call(rbind, lapply(names(CRH), function(k) { r <- CRH[[k]]; th <- theta(r$draws); dt <- dist_draws(as.matrix(TB), r$draws)
  orr <- if (r$family == "lm") OR[[r$trait]] else NULL
  data.frame(trait = traits[[r$trait]], surface = r$family,
    P_opt = r$opt[1], C_opt = r$opt[2], ratio = sprintf("1:%.2f", r$opt[2] / r$opt[1]),
    CR_ratio_range = sprintf("1:%.2f to 1:%.2f", tan(min(th) * pi / 180), tan(max(th) * pi / 180)),
    CR_on_hull_edge_pct = 100 * mean(apply(r$draws, 1, on_edge)),
    target_in_CR = if (in_poly(target, r$region)) "yes" else "no",
    dist_target = sprintf("%.0f [%.0f, %.0f]", q3(dt)[1], q3(dt)[2], q3(dt)[3]),
    OptimaRegion_opt = if (is.null(orr)) "--" else sprintf("(%.0f, %.0f)", orr$mean[1], orr$mean[2]),
    OptimaRegion_in_hull_pct = if (is.null(orr)) NA else 100 * orr$in_hull) }))
rownames(tab12) <- NULL; print(tab12)
lmk <- paste(names(traits), "lm"); pairs <- combn(lmk, 2, simplify = FALSE)
tab13 <- do.call(rbind, lapply(pairs, function(pr) { dd <- dist_draws(CRH[[pr[1]]]$draws, CRH[[pr[2]]]$draws)
  data.frame(comparison = paste(traits[[CRH[[pr[1]]]$trait]], "vs", traits[[CRH[[pr[2]]]$trait]]),
             distance_median = q3(dd)[1], lo95 = q3(dd)[2], hi95 = q3(dd)[3], ratio_ranges_overlap = { t1 <- range(theta(CRH[[pr[1]]]$draws)); t2 <- range(theta(CRH[[pr[2]]]$draws))
                             if (t1[1] <= t2[2] && t2[1] <= t1[2]) "yes" else "no" }) }))
rownames(tab13) <- NULL; print(tab13)
write_table(tab12, "TableS9_optima_regions",
  paste("Confidence regions (95\\%) for the location of each trait optimum. Optimum and confidence region (CR) from the",
        "procedure of del Castillo et al. (OptimaRegion) constrained to the observed intake domain (convex hull of",
        "individual intakes): 2000 coefficient vectors from the sampling distribution of each surface, trimmed to the 95\\%",
        "confidence region of the coefficients, each maximised exactly over the hull. CR\\_on\\_hull\\_edge: percentage of",
        "CR optima on the hull boundary. target\\_in\\_CR: whether the ING target lies inside the CR. dist\\_target: distance",
        "(ug/day) between the ING target and the trait optimum, median [95\\% range] over draws. OptimaRegion\\_opt and",
        "OptimaRegion\\_in\\_hull: mean optimum from OptRegionQuad with its origin-anchored triangular domain, and the",
        "percentage of its CR inside the observed intakes."), "tab:S9", "tables/supplementary", digits = 1)
write_table(tab13, "TableS10_optima_distances",
  "Distances (ug/day) between the optima of pairs of traits (full quadratic surfaces, hull-constrained 95\\% confidence regions): median and 95\\% range over pairs of draws, and whether the ranges of P:C ratio spanned by the two confidence regions overlap (regions on the hull edge are segments, so overlap is assessed on ratio).",
  "tab:S10", "tables/supplementary", digits = 1)
saveRDS(list(OR = OR, CRH = CRH, hull = H, triangle = list(V1 = V1, V2 = V2), tab12 = tab12, tab13 = tab13), "out/optima_regions.rds")

fig("FigS8_optima_regions", 7.2, 7.4, dir = "figures/supplementary", function() {
  par(mfrow = c(2, 2), mar = c(3.8, 3.8, 2.2, 1), oma = c(3.6, 0, 0, 0)); layout(matrix(1:4, 2, byrow = TRUE))
  panels <- list(c("lifespan lm"), c("lifetimeegg lm", "lifetimeegg quasi-Poisson"), c("dailyeggs lm", "dailyeggs quasi-Poisson"))
  for (k in 1:3) { tr <- c("lifespan", "lifetimeegg", "dailyeggs")[k]
    blank(c(0, 200), c(0, 420)); panel_lab(letters[k]); rails_draw(LEE_RATIOS, 600)
    polygon(c(0, V1[1], V2[1]), c(0, V1[2], V2[2]), border = MUTED, lty = 3)
    rect(0, 0, max(X[, 1]), max(X[, 2]), border = MUTED, lty = 3)
    polygon(H, col = adjustcolor(GRID, .35), border = MUTED)
    orr <- OR[[tr]]$xin; hh <- orr[chull(orr), , drop = FALSE]
    if (nrow(unique(round(orr, 2))) >= 3) polygon(hh, col = adjustcolor(MUTED, .55), border = MUTED, lwd = 3) else points(orr, pch = 15, col = MUTED, cex = 1.2)
    for (j in seq_along(panels[[k]])) { r <- CRH[[panels[[k]][j]]]; colj <- c(ACC, ACC2)[j]
      if (nrow(r$region) >= 3) polygon(r$region, col = adjustcolor(colj, .3), border = colj) else points(r$region, col = colj, pch = 19)
      points(r$opt[1], r$opt[2], pch = 24, bg = colj, col = "white", cex = 1.3) }
    points(TB, pch = 16, cex = .3, col = adjustcolor(INK, .25)); points(target[1], target[2], pch = 4, lwd = 2, col = INK)
    panel_title(traits[[tr]])
  }
  LEG <- list(legend = c("observed intakes (hull)", "OptimaRegion domain", "OptimaRegion 95% CR", "hull-constrained 95% CR (lm)",
                 "hull-constrained 95% CR (quasi-Poisson)", "optimum", "ING target (bootstrap draws)"),
      pch = c(15, NA, 15, 15, 15, 24, 4), lty = c(NA, 3, NA, NA, NA, NA, NA), pt.cex = c(1.6, 1, 1.6, 1.6, 1.6, 1, 1),
      col = c(adjustcolor(GRID, .6), MUTED, adjustcolor(MUTED, .6), adjustcolor(ACC, .4), adjustcolor(ACC2, .4), INK, INK), pt.bg = c(NA, NA, NA, NA, NA, MUTED, NA))
  par(mar = c(3.8, 9.5, 2.2, 1)); plot(NA, xlim = c(40, 95), ylim = c(.5, 5.5), yaxt = "n", ylab = "", xlab = expression("Angle from protein axis, " * theta * " (degrees)"), bty = "l")
  nm <- names(CRH); lab <- sapply(CRH, function(r) paste0(traits[[r$trait]], if (r$family == "lm") "" else " (log scale)"))
  axis(2, at = seq_along(nm), labels = lab, las = 1, cex.axis = .7, tick = FALSE)
  for (i in seq_along(nm)) { th <- theta(CRH[[nm[i]]]$draws); segments(min(th), i, max(th), i, col = if (CRH[[nm[i]]]$family == "lm") ACC else ACC2, lwd = 4)
    points(atan2(CRH[[nm[i]]]$opt[2], CRH[[nm[i]]]$opt[1]) * 180 / pi, i, pch = 24, bg = INK, col = "white", cex = 1.2) }
  tq <- quantile(atan2(TB[, 2], TB[, 1]) * 180 / pi, c(.025, .975)); rect(tq[1], 0, tq[2], 6, col = adjustcolor(INK, .12), border = NA)
  abline(v = atan2(target[2], target[1]) * 180 / pi, col = INK, lty = 2)
  mtext("d", side = 3, adj = -0.62, line = 0.6, font = 2, cex = .95); panel_title("Ratio span of the 95% regions")
  do.call(legend_below, c(LEG, list(ncol = 3, cex = .66)))
})
cat("S08 done: Tables S9-S13, Fig. S8\n")
