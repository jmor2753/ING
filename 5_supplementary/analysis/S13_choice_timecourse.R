# Choice experiment over 24 days (Fig. S12, Table S15)
source("R/ing_model.R"); source("R/table_utils.R"); source("R/fig_theme.R")
F0 <- readRDS("out/fits_nochoice.rds"); p <- F0$fits[[ING_BEST]]$p
set.seed(13)

oos_day <- function(day, sep) {
  ch <- load_choice(day); CH <- choice_CH(ch)
  rr <- optimize(function(lr) { pc <- p; pc$mu <- exp(lr) * p$mu; pr <- choice_opt_vec(pc, CH$FY, CH$FS, sep)
    sum((CH$lP - log(pr[, 1]))^2 + (CH$lC - log(pr[, 2]))^2) }, c(log(.1), log(3)))
  pc <- p; pc$mu <- exp(rr$minimum) * p$mu; pr <- choice_opt_vec(pc, CH$FY, CH$FS, sep)
  lo <- log(ch$C_ug_day / ch$P_ug_day); lp <- log(pr[, 2] / pr[, 1]); s180 <- ch$sucrose_gL == 180
  data.frame(days = sprintf("0--%d", day), cost = if (sep) "food-specific" else "shared", n_flies = sum(ch$n),
             obs_CP_S180 = mean(ch$C_per_P[s180]), pred_CP_S180 = mean(exp(lp[s180])),
             rho = exp(rr$minimum), ratio_cor = cor(lp, lo), ratio_RMSE = sqrt(mean((lp - lo)^2)))
}
A <- do.call(rbind, lapply(c(3, 6, 9, 12), function(d) rbind(oos_day(d, FALSE), oos_day(d, TRUE))))
print(A, row.names = FALSE)

r <- read.csv("data/lee2008_choice_raw.csv")
r$fly <- paste(r$block, r$treatment, r$replicate)
r <- r[order(r$fly, r$day), ]
win <- do.call(rbind, lapply(split(r, r$fly), function(g) {
  P <- c(0, g$P_ug_cum); C <- c(0, g$C_ug_cum); bad <- c(0, g$do_not_use)
  data.frame(fly = g$fly[1], treatment = g$treatment[1], yeast_gL = g$yeast_gL[1], sucrose_gL = g$sucrose_gL[1],
             window = seq_along(g$day), day_end = g$day, dP = diff(P) / 3, dC = diff(C) / 3,
             ok = !is.na(diff(P)) & !is.na(diff(C)) & bad[-1] == 0 & bad[-length(bad)] == 0)
}))
win <- win[win$ok, ]
trt <- unique(win[, c("treatment", "yeast_gL", "sucrose_gL")]); trt <- trt[order(trt$treatment), ]
summ <- function(w) {
  a <- aggregate(cbind(dP, dC) ~ treatment + window, w, mean)
  a$n <- aggregate(dP ~ treatment + window, w, length)$dP
  a$lCP <- log(a$dC / a$dP); a }
stats_of <- function(a) {
  a <- a[is.finite(a$lCP), ]
  slope <- sapply(split(a, a$treatment), function(g) if (nrow(g) > 2) coef(lm(lCP ~ window, g))[2] else NA)
  sdw <- tapply(a$lCP, a$window, sd)
  c(mean_slope = mean(slope, na.rm = TRUE), sd_w1 = sdw[["1"]], sd_w2 = sdw[["2"]], sd_last = sdw[[length(sdw)]],
    slope = slope)
}
S0 <- summ(win); st0 <- stats_of(S0)
B <- t(replicate(500, { wb <- do.call(rbind, lapply(split(win, win$treatment), function(g) {
    f <- unique(g$fly); fb <- sample(f, length(f), TRUE)
    do.call(rbind, lapply(seq_along(fb), function(i) { x <- g[g$fly == fb[i], ]; x$fly <- paste0(x$fly, "#", i); x })) }))
  stats_of(summ(wb)) }))
ci <- function(j) quantile(B[, j], c(.025, .975), na.rm = TRUE)
zero <- aggregate(cbind(zeroP = dP <= 0, zeroC = dC <= 0) ~ window, win, mean)
nw <- aggregate(dP ~ window, win, length); names(nw)[2] <- "n_flies"
cat(sprintf("mean within-treatment drift of log C:P per 3-day window: %.3f [%.3f, %.3f] (x%.2f per window)\n",
            st0[["mean_slope"]], ci("mean_slope")[1], ci("mean_slope")[2], exp(st0[["mean_slope"]])))
cat(sprintf("among-treatment SD of log C:P: window 1 %.2f [%.2f, %.2f], window 2 %.2f, last window %.2f [%.2f, %.2f]\n",
            st0[["sd_w1"]], ci("sd_w1")[1], ci("sd_w1")[2], st0[["sd_w2"]], st0[["sd_last"]], ci("sd_last")[1], ci("sd_last")[2]))
print(merge(nw, zero))
saveRDS(list(A = A, win = win, S = S0, stats = st0, boot = B, n_window = nw, zero = zero, trt = trt), "out/choice_timecourse.rds")

TA <- data.frame(days = A$days, cost = A$cost, n = A$n_flies, obs = A$obs_CP_S180, pred = A$pred_CP_S180,
                 rho = A$rho, r = A$ratio_cor, rmse = A$ratio_RMSE)
write_table(TA, "TableS15_choice_intervals", dir = "tables/supplementary", escape = FALSE, label = "tab:S15",
  header = c("Days", "Volume cost and capacity", "Flies", "\\makecell{Observed C:P\\\\(sucrose 180)}",
             "\\makecell{Predicted C:P\\\\(sucrose 180)}", "$\\rho$", "$r$ (log C:P)", "RMSE (log C:P)"),
  caption = paste("Out-of-sample prediction of the choice experiment of \\citet{Lee2008} from the best model fitted to",
    "the no-choice array (only the appetite scale $\\rho$ estimated), using individual intakes cumulated over",
    "different intervals. The main text uses days 0--6, the interval of the no-choice intakes and of Fig.~3 of",
    "\\citet{Lee2008}. Flies: number with a record at the end of the interval (224 flies were set up)."))

fig("FigS12_choice_timecourse", 7.2, 3.4, function() {
  par(mfrow = c(1, 2), mar = c(4, 3.6, 1.6, .6), oma = c(2.6, 0, 0, 0))
  sc <- c("180" = SEQ[4], "90" = SEQ[3], "45" = SEQ[2]); lt <- c("180" = 1, "90" = 2, "45" = 3)
  a <- S0[is.finite(S0$lCP), ]
  plot(NA, xlim = c(1, 8), ylim = range(a$lCP), xaxt = "n", yaxt = "n", bty = "l",
       xlab = "", ylab = "C:P eaten in the interval")
  mtext("Interval (days)", 1, 2.6, cex = .9)
  axis(1, 1:8, paste0(seq(0, 21, 3), "-", seq(3, 24, 3)), cex.axis = .62, las = 2)
  yt <- c(.5, 1, 2, 4, 8); axis(2, log(yt), yt)
  for (k in trt$treatment) { g <- a[a$treatment == k, ]; tr <- trt[trt$treatment == k, ]
    lines(g$window, g$lCP, col = sc[as.character(tr$sucrose_gL)], lty = lt[as.character(tr$yeast_gL)], lwd = 1.4)
    points(g$window, g$lCP, pch = 19, cex = .5, col = sc[as.character(tr$sucrose_gL)]) }
  abline(h = log(1.57), col = ACC2, lty = 3); panel_lab("a")
  barplot(nw$n_flies, names.arg = paste0(seq(0, 21, 3), "-", seq(3, 24, 3)), col = GRID, border = NA, las = 2,
          cex.names = .62, ylab = "Flies with intake recorded", xlab = "", ylim = c(0, 230))
  mtext("Interval (days)", 1, 2.6, cex = .9)
  panel_lab("b")
  legend_below(legend = c("sucrose 180", "sucrose 90", "sucrose 45", "yeast 180", "yeast 90", "yeast 45", "ING target (C:P 1.57)"),
               col = c(sc, rep(INK, 3), ACC2), lty = c(1, 1, 1, 1, 2, 3, 3), lwd = 1.4, ncol = 4)
}, dir = "figures/supplementary")
