# Appendix S3: Heliothis hybrid
source("R/ing_model.R"); source("R/table_utils.R"); source("R/fig_theme.R")
L <- readRDS("out/larvae_fits.rds"); res <- L$res; SPP <- L$SPP; CT <- readRDS("out/larvae_contrasts.rds")
HEL <- SPP[SPP$genus == "Heliothis", ]
f2 <- function(x, d = 2) formatC(x, format = "f", digits = d)
neg <- function(x) gsub("(^|\\[|, )-([0-9])", "\\1$-$\\2", x)
ci <- function(est, x, d = 2) { q <- quantile(x, c(.025, .975), na.rm = TRUE); neg(sprintf("%s [%s, %s]", f2(est, d), f2(q[1], d), f2(q[2], d))) }
lab <- function(sp) c("Heliothis virescens" = "\\textit{H.~virescens}", "Heliothis subflexa" = "\\textit{H.~subflexa}",
                      "Heliothis subflexa x virescens" = "hybrid")[[sp]]

T3 <- do.call(rbind, lapply(HEL$species, function(sp) { r <- res[[sp]]; tb <- r$table
  tb <- tb[tb$model != "Saturated (diet means)", ]; tb <- tb[order(tb$AICc), ]
  data.frame(g = lab(sp), n = nrow(r$nc), best = tb$model[1],
             dlin = f2(tb$AICc[tb$model == "Linear"] - tb$AICc[1], 1),
             lof = f2(tb$LOF_p[tb$model == "Linear"], 3), loro = f2(tb$LORO_pct[tb$model == "Linear"], 1),
             w = ci(r$lin$w, r$boot[, "w"]), dev = ci(100 * expm1(r$dev_obs), 100 * expm1(r$boot[, "dev"]), 0),
             pc = f2(exp(r$obs_choice[1] - r$obs_choice[2]), 2)) }))
hc <- CT[CT$genus == "Heliothis" & grepl("hybrid", CT$contrast), ]
write_table(T3, "AppendixS3_Table3_heliothis_hybrid", dir = "tables/appendix", escape = FALSE, label = "tab:A3.3",
  header = c("Genotype", "$n$", "Lowest AICc", "\\makecell{$\\Delta$AICc\\\\linear rule}", "\\makecell{LOF $p$\\\\(linear)}",
             "\\makecell{LORO (\\%)\\\\(linear)}", "$w=a_P/a_C$", "\\makecell{Choice\\\\deviation (\\%)}", "\\makecell{Choice\\\\P:C}"),
  caption = paste0("The \\textit{Heliothis} hybrid (\\textit{H.~subflexa} female $\\times$ \\textit{H.~virescens} male) and its parents, ",
    "analysed with the models of the main text (Table~6, Table~S14). $\\Delta$AICc linear rule: AICc of the linear rule ",
    "$a_PP+a_CC=k$ relative to the lowest-AICc model. LOF: lack of fit of the linear rule against one mean per food; LORO: ",
    "error of the predicted mean intake on a held-out food. $w$: relative weight of protein in the regulated sum. Choice ",
    "deviation: excess of the regulated sum in the choice treatment over the no-choice line. Choice P:C: geometric mean ",
    "ratio selected in the choice treatment. Hybrid relative to its parents (ratio of $w$): ",
    sprintf("%s [%s, %s] times \\textit{H.~virescens} and %s [%s, %s] times \\textit{H.~subflexa}; position on the log scale ",
            f2(hc$ratio[1]), f2(hc$lo[1]), f2(hc$hi[1]), f2(hc$ratio[2]), f2(hc$lo[2]), f2(hc$hi[2])),
    sprintf("between \\textit{H.~virescens} (0) and \\textit{H.~subflexa} (1): %s [%s, %s]. Brackets: 95\\%% bootstrap intervals ",
            f2(hc$ratio[3]), f2(hc$lo[3]), f2(hc$hi[3])),
    "(larvae resampled within treatments)."))

COL <- c(generalist = ACC, specialist = ACC2, hybrid = "#7c3aed")
comp_of <- function(lbl) { v <- as.numeric(regmatches(lbl, gregexpr("[0-9]+", lbl))[[1]])[1:2]; v / sum(v) }
lin_t <- function(psi, lk, d) exp(lk) / drop(d %*% c(cos(psi), sin(psi)))
fig("AppendixS3_Fig2_heliothis_hybrid", 7.2, 3.6, dir = "figures/appendix", function() {
  par(mfrow = c(1, 2), mar = c(3.6, 4.2, 2, 1), oma = c(2.6, 0, 0, 0))
  hy <- "Heliothis subflexa x virescens"; r <- res[[hy]]; nc <- r$nc; ch <- r$ch
  lim <- max(c(nc$P, nc$C, ch$P, ch$C)) * 1.05
  al <- seq(atan2(1, 5), atan2(5, 1), length.out = 300); dense <- cbind(sin(al), cos(al)); dense <- dense / rowSums(dense)
  blank(c(0, lim), c(0, lim), "Protein eaten (mg)", "Carbohydrate eaten (mg)"); panel_lab("a"); panel_title("Hybrid: no-choice array and choice")
  for (dl in unique(nc$diet)) { d <- comp_of(dl); segments(0, 0, 3 * lim * d[1], 3 * lim * d[2], col = GRID, lwd = .8) }
  for (f in strsplit(r$pair, " vs ")[[1]]) { d <- comp_of(f); segments(0, 0, 3 * lim * d[1], 3 * lim * d[2], col = MUTED, lty = 2) }
  points(nc$P, nc$C, pch = 16, cex = .45, col = adjustcolor(INK, .25))
  for (sp in HEL$species) { q <- res[[sp]]; X <- lin_t(q$lin$par[1], q$lin$par[2], dense) * dense
    lines(X, col = COL[HEL$diet_breadth[HEL$species == sp]], lwd = if (sp == hy) 2.2 else 1.2, lty = if (sp == hy) 1 else 3) }
  for (i in seq_along(r$dat$diets)) { y <- r$dat$lt[r$dat$trt == i]; d <- r$dat$comp[i, ]; m <- mean(y); se <- sd(y) / sqrt(length(y))
    segments(exp(m - 1.96 * se) * d[1], exp(m - 1.96 * se) * d[2], exp(m + 1.96 * se) * d[1], exp(m + 1.96 * se) * d[2], col = INK, lwd = 1.2)
    points(exp(m) * d[1], exp(m) * d[2], pch = 21, bg = "white", col = INK, cex = 1.1, lwd = 1.3) }
  lP <- log(ch$P); lC <- log(ch$C); mP <- mean(lP); mC <- mean(lC); sP <- sd(lP) / sqrt(nrow(ch)); sC <- sd(lC) / sqrt(nrow(ch))
  segments(exp(mP - 1.96 * sP), exp(mC), exp(mP + 1.96 * sP), exp(mC), col = INK); segments(exp(mP), exp(mC - 1.96 * sC), exp(mP), exp(mC + 1.96 * sC), col = INK)
  points(exp(mP), exp(mC), pch = 23, bg = COL[["hybrid"]], col = INK, cex = 1.5)
  w <- sapply(HEL$species, function(sp) res[[sp]]$lin$w); q <- sapply(HEL$species, function(sp) quantile(res[[sp]]$boot[, "w"], c(.025, .975)))
  par(mar = c(3.6, 7.4, 2, 1.5)); yy <- 3:1
  plot(NA, xlim = c(.4, 3), ylim = c(.5, 3.5), log = "x", yaxt = "n", ylab = "", bty = "l", xlab = expression("Relative weight of protein, " * italic(w)))
  abline(v = 1, col = GRID, lty = 2); axis(2, yy, c("H. virescens", "H. subflexa", "hybrid"), las = 1, font = 3, tick = FALSE, cex.axis = .8)
  segments(q[1, ], yy, q[2, ], yy, col = COL[HEL$diet_breadth], lwd = 2); points(w, yy, pch = 21, bg = COL[HEL$diet_breadth], col = "white", cex = 1.6)
  panel_lab("b"); panel_title("What is regulated")
  legend_below(ncol = 3, cex = .7, legend = c("diet mean (95% CI)", "linear rule: hybrid", "linear rule: H. virescens",
               "linear rule: H. subflexa", "foods of the choice test", "observed choice (95% CI)"),
               pch = c(21, NA, NA, NA, NA, 23), lty = c(NA, 1, 3, 3, 2, NA), lwd = c(1.3, 2.2, 1.2, 1.2, 1, NA),
               col = c(INK, COL[["hybrid"]], COL[["generalist"]], COL[["specialist"]], MUTED, INK), pt.bg = c("white", NA, NA, NA, NA, COL[["hybrid"]]))
})
cat("Appendix S3 hybrid table and figure written\n")
