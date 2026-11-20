# Appendix S3: figures and tables
source("R/ing_model.R"); source("R/table_utils.R"); source("R/fig_theme.R")
X <- readRDS("out/external_choice.rds"); JN <- readRDS("out/jensen2015_fits.rds"); PR <- readRDS("out/profile_ratio.rds")
f2 <- function(x, d = 2) formatC(x, format = "f", digits = d)
lab_set <- function(s) gsub("&", "\\\\&", sub(" F$", ", females", sub(" M$", ", males", s)))

H <- X$hyp; P <- X$pred
T1 <- do.call(rbind, lapply(seq_len(nrow(H)), function(i) { h <- H[i, ]; p <- P[P$set == h$set, ]
  data.frame(set = lab_set(h$set), n = h$n_trt, obs = sprintf("%s--%s", f2(min(p$obs_CP)), f2(max(p$obs_CP))),
             pred = sprintf("%s--%s", f2(min(p$pred_CP)), f2(max(p$pred_CP))), rho = f2(p$rho[1]),
             a0 = f2(h$AICc_base, 1), a1 = f2(h$AICc_suc, 1), a2 = f2(h$AICc_cap, 1),
             g = f2(h$gamma), v = f2(h$VAmax)) }))
write_table(T1, "AppendixS3_Table1_external_choice", dir = "tables/appendix", escape = FALSE, sideways = TRUE, label = "tab:A3.1",
  header = c("Data set", "Treatments", "\\makecell{Observed\\\\C:P}", "\\makecell{Predicted\\\\C:P}", "$\\rho$",
             "\\makecell{AICc\\\\(prediction)}", "\\makecell{AICc\\\\(+ sucrose reward)}", "\\makecell{AICc\\\\(+ protein-food limit)}",
             "$\\gamma$", "\\makecell{Limit\\\\(\\textmu L\\,day$^{-1}$)}"),
  caption = paste("Out-of-sample prediction of three choice experiments from the best ING model fitted to the no-choice array",
    "of \\citet{Lee2008}, with only the appetite scale $\\rho$ estimated, and two one-parameter explanations of the",
    "carbohydrate bias fitted to each choice data set: a reward $\\gamma$ per \\textmu L of a pure-sucrose food, and a",
    "volume limit on the protein-richer food. AICc from the squared errors of log protein and log carbohydrate intake.",
    "Ranges of C:P are over treatments. In \\citet{MorimotoWigby2016} no pure-sucrose food was offered, and the flies",
    "drank about 2.5 times the volume of Lee's flies, beyond the capacity estimated from the Lee data, so the full",
    "prediction is not transferable there (see text)."))

T2 <- do.call(rbind, lapply(c("F", "M"), function(sx) { r <- JN[[sx]]; t <- head(r$table, 5); pf <- r$profile
  at <- function(v) f2(approx(pf$r, pf$best, v)$y, 1)
  data.frame(sex = c(if (sx == "F") "Females" else "Males", rep("", nrow(t) - 1)), model = t$model, k = t$k,
             d = f2(t$dAICc), cp = f2(t$CP), om = f2(t$omega), vm = ifelse(is.na(t$Vmax), "--", f2(t$Vmax)),
             p1 = c(at(1.57), rep("", nrow(t) - 1)), p4 = c(at(4), rep("", nrow(t) - 1))) }))
write_table(T2, "AppendixS3_Table2_jensen_nochoice", dir = "tables/appendix", escape = FALSE, label = "tab:A3.2",
  header = c("Sex", "Model", "$k$", "$\\Delta$AICc", "\\makecell{Target\\\\C:P}", "$\\omega$", "\\makecell{$V_{\\max}$\\\\(\\textmu L\\,day$^{-1}$)}",
             "\\makecell{$\\Delta$nLL at\\\\C:P 1.57}", "\\makecell{$\\Delta$nLL at\\\\C:P 4}"),
  caption = paste("ING fitted to the no-choice array of \\citet{Jensen2015} (28 holidic diets, treatment means and standard",
    "errors of each sex): the five lowest-AICc models of the 64 for each sex. Last two columns: profile negative",
    "log-likelihood of the best model at a target C:P of 1.57 (the latent target of the Lee flies) and 4 (the",
    "self-selected ratio), relative to its minimum. Model codes as in Table~1 of the main text."))

SETC <- c("Lee et al. 2008 F" = MUTED, "Jensen et al. 2015 F" = ACC2, "Jensen et al. 2015 M" = ACC, "Morimoto & Wigby 2016 M" = INK)
fig("AppendixS3_Fig1_external_tests", 7.2, 7.8, dir = "figures/appendix", function() {
  par(mfrow = c(2, 2), mar = c(3.8, 4, 2.2, 1), oma = c(4.4, 0, 0, 0))
  plot(NA, xlim = c(.3, 7), ylim = c(.3, 7), log = "xy", xlab = "Predicted C:P (Lee no-choice fit)", ylab = "Observed C:P", bty = "l")
  abline(0, 1, col = GRID, lty = 2); abline(h = 1.57, col = ACC2, lty = 3)
  for (s in names(SETC)) { p <- P[P$set == s, ]
    segments(p$lo_CP, p$obs_CP, p$hi_CP, p$obs_CP, col = adjustcolor(SETC[s], .5))
    points(p$pred_CP, p$obs_CP, pch = 21, bg = SETC[s], col = "white", cex = 1.3) }
  panel_lab("a"); panel_title("Prediction of three choice experiments")
  plot(NA, xlim = c(1, 12), ylim = c(1, 7), log = "xy", xlab = "Predicted C:P", ylab = "Observed C:P", bty = "l")
  abline(0, 1, col = GRID, lty = 2)
  for (s in c("Jensen et al. 2015 F", "Jensen et al. 2015 M")) { p <- P[P$set == s, ]; h <- H[H$set == s, ]
    pc <- as.numeric(strsplit(h$CP_cap, " ")[[1]])
    points(p$pred_CP, p$obs_CP, pch = 1, col = SETC[s], cex = 1.2)
    points(pc, p$obs_CP, pch = 19, col = SETC[s], cex = 1.2)
    text(pc, p$obs_CP, sub(".*@(\\d+) vs 0:1@(\\d+)", "\\1/\\2", p$group), pos = 4, cex = .55, col = SETC[s]) }
  panel_lab("b"); panel_title("Jensen et al. 2015: limit on the protein food")
  plot(NA, xlim = c(40, 400), ylim = c(0, 1.3), log = "x", xaxt = "n", xlab = expression("Food concentration (g L"^-1*")"),
       ylab = expression("Volume eaten ("*mu*"L day"^-1*")"), bty = "l"); axis(1, c(45, 90, 180, 360))
  for (sx in c("F", "M")) { d <- JN[[sx]]$dat; v <- exp(d$ybar) / d$F
    for (r in unique(d$tr$ratio)) { k <- d$tr$ratio == r; o <- order(d$F[k])
      lines(d$F[k][o], v[k][o], col = adjustcolor(c(F = ACC2, M = ACC)[sx], .6)) }
    points(d$F, v, pch = 16, cex = .6, col = c(F = ACC2, M = ACC)[sx])
    abline(h = vmax_of(JN[[sx]]$best$p), lty = 2, col = c(F = ACC2, M = ACC)[sx]) }
  panel_lab("c"); panel_title("Jensen et al. 2015: no-choice volumes (one line per rail)")
  plot(NA, xlim = c(1, 6), ylim = c(0, 40), xlab = "Target C:P ratio", ylab = expression(Delta * " negative log-likelihood"), bty = "l")
  abline(h = qchisq(.95, 1) / 2, col = GRID, lty = 2); abline(v = 4, col = MUTED, lty = 3)
  s <- PR[PR$model == ING_BEST & !PR$joint, ]; s <- s[order(s$r), ]; lines(s$r, s$dnll, col = MUTED, lwd = 2)
  for (sx in c("F", "M")) { pf <- JN[[sx]]$profile; lines(pf$r, pf$best, col = c(F = ACC2, M = ACC)[sx], lwd = 2) }
  panel_lab("d"); panel_title("Profile of the target ratio (no-choice arrays)")
  legend_below(ncol = 2, cex = .68, pt.cex = 1.3, pch = c(21, 21, 21, 21, 1, 19, NA, NA), lty = c(NA, NA, NA, NA, NA, NA, 2, 3),
    pt.bg = c(SETC, NA, NA, NA, NA), col = c(rep("white", 4), INK, INK, INK, MUTED),
    legend = c("Lee et al. 2008 (females)", "Jensen et al. 2015, females", "Jensen et al. 2015, males",
               "Morimoto & Wigby 2016 (males)", "(b) prediction from the Lee fit", "(b) with a protein-food limit",
               "(c) estimated capacity", "(d) C:P = 4"))
})
cat("Appendix S3 outputs written\n")
