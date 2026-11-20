# Appendix S2: verification of the slope jump at breakpoints
source("R/ing_model.R"); source("R/table_utils.R")
p <- readRDS("out/fits_nochoice.rds")$fits[[ING_BEST]]$p
kap <- exp(p$eta); q <- p$Qz[1, 2]
Jt <- function(t, c, F) { x <- matrix(t * c, 1); V <- sum(W_FOOD * c) * t / F; loss_R(x, p) + p$l1 * V + p$l2 * V^2 }

al <- seq(0.05, atan(1.9), length.out = 9); C9 <- cbind(sin(al), cos(al)); C9 <- C9 / rowSums(C9)
off <- 1e-2; rows <- list()
for (F in c(360)) for (i in 1:9) for (j in 1:2) {
  c1 <- C9[i, ]; b <- p$mu[j] / c1[j]; k <- 3 - j
  sR <- (-3 * Jt(b, c1, F) + 4 * Jt(b + off, c1, F) - Jt(b + 2 * off, c1, F)) / (2 * off)
  sL <- ( 3 * Jt(b, c1, F) - 4 * Jt(b - off, c1, F) + Jt(b - 2 * off, c1, F)) / (2 * off)
  zk <- (b * c1[k] - p$mu[k]) / p$mu[k]; dk <- zk * (if (zk > 0) kap[k] else 1)
  ana <- 2 * (kap[j] - 1) * c1[j] / p$mu[j] * q * dk
  rows[[length(rows) + 1]] <- data.frame(F_gL = F, rail_PC = round(c1[1] / c1[2], 3), crossing = c("P", "C")[j],
    jump_numeric = signif(sR - sL, 6), jump_formula = signif(ana, 6), relative_error = signif(abs((sR - sL) - ana) / abs(ana), 3),
    kink = ifelse(ana > 0, "convex (corner possible)", "concave (switch possible)"))
}
t1 <- do.call(rbind, rows)
cat("max relative error of Eq. S3:", signif(max(t1$relative_error), 3), "\n")

al <- seq(0.01, atan(1.95), length.out = 2000); D <- cbind(sin(al), cos(al)); D <- D / rowSums(D)
runs <- function(flag, PC, X, F, nutrient) { r <- rle(flag); e <- cumsum(r$lengths); b <- e - r$lengths + 1
  do.call(rbind, lapply(which(r$values), function(k) data.frame(F_gL = F, feature = paste("corner:", nutrient, "held at target"),
    n_rails = r$lengths[k], PC_from = sprintf("1:%.2f", 1 / PC[e[k]]), PC_to = sprintf("1:%.2f", 1 / PC[b[k]])))) }
t2 <- do.call(rbind, lapply(c(45, 90, 180, 360), function(F) {
  X <- rail_opt(p, D, rep(F, nrow(D))) * D; PC <- D[, 1] / D[, 2]
  cC <- abs(X[, 2] - p$mu[2]) < 1e-6 * p$mu[2]; cP <- abs(X[, 1] - p$mu[1]) < 1e-6 * p$mu[1]
  st <- sqrt(rowSums(diff(X)^2)); cand <- which(st > 5 * median(st))
  is_jump <- function(k) { lo <- al[k]; hi <- al[k + 1]; s0 <- st[k]
    xo <- function(a) { d <- c(sin(a), cos(a)); d <- d / sum(d); rail_opt(p, matrix(d, 1), F) * d }
    for (it in 1:40) { m <- (lo + hi) / 2; xl <- xo(lo); xm <- xo(m); xh <- xo(hi)
      if (sum((xm - xl)^2) > sum((xh - xm)^2)) hi <- m else lo <- m }
    sqrt(sum((xo(hi) - xo(lo))^2)) > 0.5 * s0 }
  jmp <- Filter(is_jump, cand)
  Vd <- as.vector(D %*% W_FOOD) / F * rowSums(X); cV <- is.finite(vmax_of(p)) & abs(Vd - vmax_of(p)) < 1e-6 * vmax_of(p)
  rc <- runs(cV, PC, X, F, "x"); if (!is.null(rc)) rc$feature <- sprintf("at intake capacity (V = %.2f uL/day)", vmax_of(p))
  rbind(runs(cC, PC, X, F, "carbohydrate"), runs(cP, PC, X, F, "protein"), rc,
        if (length(jmp)) data.frame(F_gL = F, feature = sprintf("switch: jump of %.2f ug (median step %.3f)", st[jmp], median(st)),
                                    n_rails = NA, PC_from = sprintf("1:%.2f", 1 / PC[jmp + 1]), PC_to = sprintf("1:%.2f", 1 / PC[jmp]))) }))
print(t2, row.names = FALSE)
dir.create("tables/appendix", recursive = TRUE, showWarnings = FALSE)
t1o <- t1; t1o$jump_numeric <- formatC(t1$jump_numeric, format = "e", digits = 3)
t1o$jump_formula <- formatC(t1$jump_formula, format = "e", digits = 3)
t1o$relative_error <- formatC(t1$relative_error, format = "e", digits = 1)
write_table(t1o, "AppendixS2_Table1_slope_jump_check",
  paste("Slope jump of the objective where a rail crosses a target coordinate: numerical value against",
        "Eq. S3, for the best model at each rail and concentration. Upward jumps allow corner solutions;",
        "downward jumps allow abrupt switches."), "tab:S2.1", "tables/appendix", digits = 2)
write_table(t2, "AppendixS2_Table2_corners_and_switches",
  paste("Corner solutions (one nutrient held exactly at its target over a run of rails) and abrupt switches",
        "predicted by the best model, from a dense set of rails at each food concentration."),
  "tab:S2.2", "tables/appendix", digits = 2)
stopifnot(max(t1$relative_error) < 1e-5)
