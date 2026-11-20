# Fig. 4: prediction of the choice experiment
source("R/ing_model.R"); source("R/fig_theme.R")
F0 <- readRDS("out/fits_nochoice.rds"); dat <- F0$dat; fits <- F0$fits
best <- fits[[ING_BEST]]; symm <- fits[[ING_SYMM]]
ch <- load_choice(6); CP180 <- mean(ch$C_per_P[ch$sucrose_gL == 180])

pred_choice <- function(p, sep = FALSE) {
  rr <- optimize(function(lr) { pc <- p; pc$mu <- exp(lr) * p$mu; pr <- choice_opt_vec(pc, ch$yeast_gL, ch$sucrose_gL, sep)
    sum((log(ch$P_ug_day) - log(pr[, 1]))^2 + (log(ch$C_ug_day) - log(pr[, 2]))^2) }, c(log(.1), log(3)))
  pc <- p; pc$mu <- exp(rr$minimum) * p$mu; list(pr = choice_opt_vec(pc, ch$yeast_gL, ch$sucrose_gL, sep), rho = exp(rr$minimum)) }
PB <- pred_choice(best$p); PS <- pred_choice(best$p, sep = TRUE); PI <- pred_choice(fits[["Iso+L1"]]$p)
CB <- readRDS("out/cv_boot.rds"); BP <- boot_params(CB, best)
bq <- function(sep) { A <- lapply(BP, function(pb) pred_choice(pb, sep)$pr)
  P <- sapply(A, function(a) a[, 1]); C <- sapply(A, function(a) a[, 2]); R <- log(C / P)
  q <- function(M) apply(M, 1, quantile, c(.025, .975), na.rm = TRUE)
  list(P = q(P), C = q(C), R = q(R)) }
QB <- bq(FALSE); QS <- bq(TRUE)
lo_obs <- log(ch$C_ug_day / ch$P_ug_day)
write.csv(data.frame(sucrose = ch$sucrose_gL, yeast = ch$yeast_gL, obs_logCP = lo_obs,
  pred_shared = log(PB$pr[, 2] / PB$pr[, 1]), lo_shared = QB$R[1, ], hi_shared = QB$R[2, ],
  pred_food = log(PS$pr[, 2] / PS$pr[, 1]), lo_food = QS$R[1, ], hi_food = QS$R[2, ]),
  "out/table_choice_boot_intervals.csv", row.names = FALSE)
fig("Fig6_choice_prediction", 7.2, 7.4, function() {
  par(mfrow = c(2, 2), mar = c(3.6, 3.8, 2, 1), oma = c(4.2, 0, 0, 0)); ycol <- c(`180` = SEQ[4], `90` = SEQ[3], `45` = SEQ[2])
  for (k in 1:3) { S <- c(180, 90, 45)[k]; ii <- which(ch$sucrose_gL == S)
    blank(c(0, 80), c(0, 100), expression("Protein eaten ("*mu*"g day"^-1*")"), expression("Carbohydrate eaten ("*mu*"g day"^-1*")")); panel_lab(letters[k])
    segments(0, 0, 22.5, 90, col = "#9ca3af", lty = 2)
    segments(ch$P_ug_day[ii], ch$C_ug_day[ii], PB$pr[ii, 1], PB$pr[ii, 2], col = GRID)
    segments(QB$P[1, ii], PB$pr[ii, 2], QB$P[2, ii], PB$pr[ii, 2], col = ACC, lwd = 1)
    segments(PB$pr[ii, 1], QB$C[1, ii], PB$pr[ii, 1], QB$C[2, ii], col = ACC, lwd = 1)
    points(PB$pr[ii, 1], PB$pr[ii, 2], pch = 21, bg = "white", col = ACC, lwd = 1.5, cex = 1.3)
    oc <- ycol[as.character(ch$yeast_gL[ii])]
    segments(ch$P_ug_day[ii] * exp(-1.96 * ch$se_lP[ii]), ch$C_ug_day[ii], ch$P_ug_day[ii] * exp(1.96 * ch$se_lP[ii]), ch$C_ug_day[ii], col = oc, lwd = 1.2)
    segments(ch$P_ug_day[ii], ch$C_ug_day[ii] * exp(-1.96 * ch$se_lC[ii]), ch$P_ug_day[ii], ch$C_ug_day[ii] * exp(1.96 * ch$se_lC[ii]), col = oc, lwd = 1.2)
    points(ch$P_ug_day[ii], ch$C_ug_day[ii], pch = 19, col = oc, cex = 1.3)
    text(ch$P_ug_day[ii], ch$C_ug_day[ii], paste0("Y", ch$yeast_gL[ii]), pos = 4, cex = .62, col = INK)
    mtext(paste0("Sucrose ", S, " g/L"), 3, .3, cex = .78)
  }
  lo <- log(ch$C_ug_day / ch$P_ug_day)
  plot(NA, xlim = c(-0.4, 1.6), ylim = c(-0.4, 1.6), xlab = "Predicted log(C:P)", ylab = "Observed log(C:P)", bty = "l"); panel_lab("d")
  abline(0, 1, col = GRID, lty = 2)
  points(log(PI$pr[, 2] / PI$pr[, 1]), lo, pch = 24, col = MUTED, bg = "white", cex = 1.1)
  segments(log(PB$pr[, 2] / PB$pr[, 1]), lo + .012 - 1.96 * ch$se_lCP, log(PB$pr[, 2] / PB$pr[, 1]), lo + .012 + 1.96 * ch$se_lCP, col = GRID, lwd = 1)
  segments(log(PS$pr[, 2] / PS$pr[, 1]), lo - .012 - 1.96 * ch$se_lCP, log(PS$pr[, 2] / PS$pr[, 1]), lo - .012 + 1.96 * ch$se_lCP, col = GRID, lwd = 1)
  segments(QB$R[1, ], lo + .012, QB$R[2, ], lo + .012, col = ACC, lwd = 1)
  segments(QS$R[1, ], lo - .012, QS$R[2, ], lo - .012, col = ACC2, lwd = 1)
  points(log(PB$pr[, 2] / PB$pr[, 1]), lo + .012, pch = 21, col = ACC, bg = ACC, cex = 1.1)
  points(log(PS$pr[, 2] / PS$pr[, 1]), lo - .012, pch = 22, col = ACC2, bg = ACC2, cex = 1.05)
  mtext("Ratio prediction across nine choice treatments", 3, .3, cex = .78)
  legend_below(ncol = 3, pch = c(19, 19, 19, 21, NA, NA, 21, 22, 24), lty = c(NA, NA, NA, 1, 2, NA, NA, NA, NA),
         col = c(ycol, ACC, "#9ca3af", NA, ACC, ACC2, MUTED), pt.bg = c(NA, NA, NA, "white", NA, NA, ACC, ACC2, "white"),
         legend = c("(a-c) observed, yeast 180 g/L", "(a-c) observed, yeast 90 g/L", "(a-c) observed, yeast 45 g/L (95% interval)",
                    "(a-c) ING prediction, 95% interval", "(a-c) P:C 1:4", "",sprintf("(d) best model, shared volume (r = %.2f)", cor(log(PB$pr[, 2] / PB$pr[, 1]), lo)),
                    sprintf("(d) best model, per food (r = %.2f)", cor(log(PS$pr[, 2] / PS$pr[, 1]), lo)),
                    sprintf("(d) isotropic, linear cost (r = %.2f)", cor(log(PI$pr[, 2] / PI$pr[, 1]), lo))))
})
cat("Figure 6 written\n")
