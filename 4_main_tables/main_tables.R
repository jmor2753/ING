# Tables 1-2 and Tables S19-S21
source("R/ing_model.R"); source("R/table_utils.R")
ch <- load_choice(6); CP180 <- mean(ch$C_per_P[ch$sucrose_gL == 180])
F0 <- readRDS("out/fits_nochoice.rds"); fits <- F0$fits
cv <- read.csv("out/table_cv.csv")
t1 <- do.call(rbind, lapply(names(fits), function(n) { f <- fits[[n]]; r <- cv[cv$model == n, ]
  data.frame(Model = n, k = f$k, AICc = f$AICc, tau = f$tau, CP_target = f$desc[["C_per_P"]],
             omega = f$desc[["omega"]], LORO_pct = r$pct_rail, LOCO_pct = r$pct_conc) }))
t1$dAICc <- t1$AICc - min(t1$AICc); t1 <- t1[order(t1$AICc), c("Model", "k", "AICc", "dAICc", "tau", "CP_target", "omega", "LORO_pct", "LOCO_pct")]
t1_full <- t1
f_num <- function(x, d = 2) ifelse(is.na(x), "--", formatC(x, format = "f", digits = d))
f_ci  <- function(e, lo, hi, d = 2) ifelse(is.na(e), "--", ifelse(abs(hi - lo) < 1e-12, paste0(f_sig(e), " (fixed)"),
                                    sprintf("%s [%s, %s]", f_sig(e), f_sig(lo), f_sig(hi))))
f_sig <- function(x) ifelse(abs(x) < 1e-3 & x != 0, "$<$0.001", formatC(signif(x, 3), format = "fg", digits = 3, flag = "#"))
f_sig <- function(x) sapply(x, function(v) if (is.na(v)) "--" else if (v != 0 && abs(v) < 1e-3) "$<$0.001" else
  sub("\\.$", "", formatC(signif(v, 3), format = "fg", digits = 3, flag = "#")))
MODEL_KEY <- paste("Model codes: Iso/Aniso, isotropic or anisotropic regulation; AsymP, AsymC, AsymPC, excess--deficit",
  "asymmetry in protein, carbohydrate or both; L1, L2, L12, linear, quadratic or both volume costs; NoAcq, no",
  "acquisition cost; Cap, intake capacity.")
t1 <- t1_full[t1_full$dAICc <= 10 | t1_full$Model %in% c("Iso+NoAcq", "Aniso+NoAcq"), ]
T1 <- data.frame(Model = t1$Model, k = t1$k, AICc = f_num(t1$AICc), dAICc = f_num(t1$dAICc), tau = f_num(t1$tau, 3),
                 CP = f_num(t1$CP_target), omega = f_num(t1$omega), LORO = f_num(t1$LORO_pct, 1), LOCO = f_num(t1$LOCO_pct, 1))
HDR1 <- c("Model", "$k$", "AICc", "$\\Delta$AICc", "$\\tau$", "Target C:P", "$\\omega$", "LORO (\\%)", "LOCO (\\%)")
write_table(T1, "Table1_model_comparison", paste("ING models fitted to the no-choice array of \\citet{Lee2008}: models within 10 AICc units",
  "of the best, and the two models without acquisition costs or capacity (all 64 models in Table~S1). $k$: number",
  "of parameters, including the among-treatment SD $\\tau$. Target C:P: ratio of the latent intake target;",
  "$\\omega$: anisotropy (0.5 = isotropic). LORO and LOCO: leave-one-rail-out and leave-one-concentration-out",
  "prediction error of treatment-mean intake.", MODEL_KEY), "tab:1", dir = "tables/main", header = HDR1)
T10 <- data.frame(Model = t1_full$Model, k = t1_full$k, AICc = f_num(t1_full$AICc), dAICc = f_num(t1_full$dAICc), tau = f_num(t1_full$tau, 3),
                  CP = f_num(t1_full$CP_target), omega = f_num(t1_full$omega), LORO = f_num(t1_full$LORO_pct, 1), LOCO = f_num(t1_full$LOCO_pct, 1))
write_table(T10, "TableS1_full_model_set", paste("All 64 ING models fitted to the no-choice array of \\citet{Lee2008}, ordered by AICc.",
  "Columns as in Table~1.", MODEL_KEY), "tab:S1", dir = "tables/supplementary", header = HDR1)
CB <- readRDS("out/cv_boot.rds")
wrap <- wrap_axis
ci <- do.call(rbind, lapply(names(CB$boots), function(n) { b <- CB$boots[[n]]; b <- b[, !grepl("^\\.par", colnames(b))]; est <- c(fits[[n]]$desc, tau = fits[[n]]$tau)
  b[, "phi_deg"] <- wrap(b[, "phi_deg"]); est["phi_deg"] <- wrap(est["phi_deg"])
  data.frame(model = n, parameter = colnames(b), estimate = est[colnames(b)],
             lo = apply(b, 2, quantile, .025, na.rm = TRUE), hi = apply(b, 2, quantile, .975, na.rm = TRUE)) }))
ci$estimate[ci$model == "Iso+L1" & ci$parameter == "phi_deg"] <- NA
write.csv(ci, "out/table_param_ci.csv", row.names = FALSE)
keep <- c("muP", "muC", "C_per_P", "omega", "phi_deg", "excess_deficit_P", "excess_deficit_C", "lambda1", "lambda2", "Vmax", "tau")
lab2 <- c(muP = "$\\mu_P$ (\\textmu g\\,day$^{-1}$)", muC = "$\\mu_C$ (\\textmu g\\,day$^{-1}$)", C_per_P = "Target C:P",
          omega = "Anisotropy $\\omega$", phi_deg = "Principal axis $\\phi$ ($^\\circ$)",
          excess_deficit_P = "Protein excess:deficit weight", excess_deficit_C = "Carbohydrate excess:deficit weight",
          lambda1 = "Linear volume cost $\\lambda_1$", lambda2 = "Quadratic volume cost $\\lambda_2$",
          Vmax = "Intake capacity $V_{\\max}$ (\\textmu L\\,day$^{-1}$)", tau = "Among-treatment SD $\\tau$")
t2 <- ci[ci$parameter %in% keep, ]; t2$CI <- f_ci(t2$estimate, t2$lo, t2$hi)
t2 <- reshape(t2[, c("model", "parameter", "CI")], idvar = "parameter", timevar = "model", direction = "wide")
names(t2) <- sub("CI.", "", names(t2)); t2 <- t2[match(keep, t2$parameter), ]
mods2 <- setdiff(names(t2), "parameter"); t2[is.na(t2)] <- "--"
T2 <- data.frame(Parameter = lab2[t2$parameter], t2[, mods2], check.names = FALSE)
T2[, mods2] <- lapply(T2[, mods2], function(v) gsub("(^|\\[|, )-([0-9])", "\\1$-$\\2", v))
write_table(T2, "Table2_parameters", paste("Parameter estimates of four ING architectures, with 95\\% parametric-bootstrap",
  "intervals (200 replicates); the best model is Aniso+AsymPC+L2+Cap. Excess:deficit weights are $e^{2\\eta_j}$, the",
  "weight of an excess of nutrient $j$ relative to an equal proportional deficit (1 = symmetric). Fixed: parameter",
  "absent from the model.", MODEL_KEY), "tab:2", dir = "tables/main", header = c("Parameter", mods2), escape = FALSE)
oos <- read.csv("out/table_choice_oos.csv")
sh <- oos[oos$cost == "shared", c("model", "rho", "ratio_cor", "ratio_RMSE", "CP_S180", "p1_ratio_cor", "p1_ratio_RMSE")]
fs <- oos[oos$cost == "separable", c("model", "ratio_cor", "ratio_RMSE", "CP_S180", "p1_ratio_cor", "p1_ratio_RMSE", "p1_CP_S180")]
names(sh)[-1] <- c("rho", "r_shared", "RMSE_shared", "CP_S180_shared", "p1_r_shared", "p1_RMSE_shared")
names(fs)[-1] <- c("r_food", "RMSE_food", "CP_S180_food", "p1_r_food", "p1_RMSE_food", "p1_CP_food")
t3 <- merge(merge(sh, fs, by = "model"), t1_full[, c("Model", "dAICc")], by.x = "model", by.y = "Model")
t3 <- t3[t3$dAICc <= 10 | t3$model %in% ING_REP, ]
t3 <- t3[order(t3$dAICc), ]
T3 <- data.frame(Model = t3$model, dAICc = f_num(t3$dAICc),
                 pr1 = f_num(t3$p1_r_shared), pe1 = f_num(t3$p1_RMSE_shared), pr2 = f_num(t3$p1_r_food), pe2 = f_num(t3$p1_RMSE_food), pc2 = f_num(t3$p1_CP_food),
                 rho = f_num(t3$rho), r1 = f_num(t3$r_shared), e1 = f_num(t3$RMSE_shared),
                 r2 = f_num(t3$r_food), e2 = f_num(t3$RMSE_food), c2 = f_num(t3$CP_S180_food))
write_table(T3, "Table3_choice_prediction", paste0("Prediction of the nine choice treatments from models fitted to ",
  "the no-choice array only, with volume cost and capacity shared between the two foods (shared) or applied to each food ",
  "separately (per food). Prospective: no parameter is informed by the choice data (appetite scale $\\rho = 1$). ",
  "Calibrated: the target is rescaled by a single appetite scale $\\rho$ fitted to the choice amounts, which also changes the predicted ratios. ",
  "$\\Delta$AICc: no-choice fit (Table~1); $r$ and RMSE: correlation and root-mean-square error between observed and predicted log(C:P); ",
  "C:P$_{180}$: mean predicted C:P at sucrose 180~g\\,L$^{-1}$ (observed ", sprintf("%.2f", CP180), "). Models within 10 ",
  "AICc units of the best, plus representative architectures. ", MODEL_KEY), "tab:3", dir = "tables/main",
  header = c("Model", "$\\Delta$AICc",
             "\\makecell{$r$\\\\shared\\\\$\\rho=1$}", "\\makecell{RMSE\\\\shared\\\\$\\rho=1$}", "\\makecell{$r$\\\\per food\\\\$\\rho=1$}",
             "\\makecell{RMSE\\\\per food\\\\$\\rho=1$}", "\\makecell{C:P$_{180}$\\\\per food\\\\$\\rho=1$}",
             "\\makecell{fitted\\\\$\\rho$}", "\\makecell{$r$\\\\shared\\\\fitted $\\rho$}", "\\makecell{RMSE\\\\shared\\\\fitted $\\rho$}",
             "\\makecell{$r$\\\\per food\\\\fitted $\\rho$}", "\\makecell{RMSE\\\\per food\\\\fitted $\\rho$}", "\\makecell{C:P$_{180}$\\\\per food\\\\fitted $\\rho$}"))
tr <- read.csv("out/table_trait_nutrigon.csv")
TRL <- c(lifespan = "Lifespan", lifetimeegg = "Lifetime eggs", dailyeggs = "Egg-production rate")
T4 <- data.frame(Trait = TRL[tr$trait], Surface = ifelse(tr$model == "lm", "quadratic", "quasi-Poisson"), Shape = tr$surface_shape,
                 P = f_num(tr$P_opt, 1), C = f_num(tr$C_opt, 1), theta = f_num(tr$theta_deg, 1), ratio = tr$ratio,
                 region = sprintf("%s--%s", f_num(tr$region_theta_lo, 1), f_num(tr$region_theta_hi, 1)),
                 rr = sub(" to ", "--", tr$region_ratio), edge = f_num(tr$dist_to_hull_edge, 2))
write_table(T4, "Table4_trait_optima", paste("Performance surfaces fitted with \\texttt{nutrigon} to individual flies (full quadratic,",
  "as in \\citealp{Lee2008}; quasi-Poisson with log link as a sensitivity analysis for the egg traits). Optima are",
  "closed-form stationary points when these are maxima inside the intake domain, and otherwise the highest point of",
  "the top-5\\% peak region. $\\theta$: angle from the protein axis (P:C $=1:\\tan\\theta$). Edge distance: distance of",
  "the optimum to the boundary of the intake domain (\\textmu g\\,day$^{-1}$)."), "tab:4", dir = "tables/main",
  header = c("Trait", "Surface", "Shape", "$P_{\\mathrm{opt}}$", "$C_{\\mathrm{opt}}$", "$\\theta$ ($^\\circ$)", "P:C",
             "\\makecell{Peak region\\\\$\\theta$ ($^\\circ$)}", "\\makecell{Peak region\\\\P:C}", "\\makecell{Edge\\\\distance}"), escape = FALSE)
ivt <- read.csv("out/table_ing_vs_traits.csv")
lab5 <- c(TRL, `choice (sucrose 180, observed)` = "Self-selected (choice, sucrose 180~g\\,L$^{-1}$)")
T5 <- data.frame(Optimum = ifelse(ivt$trait %in% names(lab5), lab5[ivt$trait], ivt$trait), theta = f_num(ivt$theta_trait, 1),
                 ratio = ivt$trait_ratio,
                 sep = ifelse(is.na(ivt$lo95), f_num(ivt$separation_deg, 1), sprintf("%s [%s, %s]", f_num(ivt$separation_deg, 1), f_num(ivt$lo95, 1), f_num(ivt$hi95, 1))))
write_table(T5, "Table5_ING_vs_traits", paste0("Angular separation between the latent intake target of the best ING model ",
  "($\\theta=", f_num(ivt$theta_ING[1], 1), "^\\circ$, P:C ", ivt$ING_ratio[1], ") and each trait optimum (full quadratic ",
  "surfaces), with 95\\% intervals from a joint stratified bootstrap (flies resampled within treatments; ING and ",
  "\\texttt{nutrigon} refitted on each resample; $B=200$)."), "tab:5", dir = "tables/main",
  header = c("Optimum", "$\\theta$ ($^\\circ$)", "P:C", "Separation from target ($^\\circ$)"), escape = FALSE)
cat("main tables written to tables/main\n")
