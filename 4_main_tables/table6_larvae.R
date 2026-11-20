# Table 3: caterpillar models
source("R/ing_model.R"); source("R/table_utils.R")
L <- readRDS("out/larvae_fits.rds"); res <- L$res; SPP <- L$SPP
SPP_all <- SPP; SPP <- SPP[SPP$diet_breadth != "hybrid", ]
B <- nrow(res[[1]]$boot)
f2 <- function(x, d = 2) formatC(x, format = "f", digits = d)
neg <- function(x) gsub("(^|\\[|, )-([0-9])", "\\1$-$\\2", x)
row_ci <- function(est, x, d = 2) { q <- quantile(x, c(.025, .975), na.rm = TRUE); neg(sprintf("%s [%s, %s]", f2(est, d), f2(q[1], d), f2(q[2], d))) }
sp_lab <- function(sp) sub("~x~", "} $\\\\times$ \\\\textit{",
  paste0("\\textit{", gsub(" ", "~", SPP_all$short[SPP_all$species == sp]), "}"))
T6 <- do.call(rbind, lapply(SPP$species, function(sp) {
  r <- res[[sp]]; b <- r$boot; tb <- r$table; tb <- tb[tb$model != "Saturated (diet means)", ]; o <- order(tb$AICc)
  bestm <- tb$model[o[1]]
  best_lab <- if (bestm == "Linear") sprintf("Linear (%s over %s)", f2(tb$dAICc[o[2]] - tb$dAICc[o[1]], 1), tb$model[o[2]]) else
    sprintf("%s (%s over Linear)", bestm, f2(tb$AICc[tb$model == "Linear"] - tb$AICc[o[1]], 1))
  data.frame(Species = sp_lab(sp), Breadth = SPP$diet_breadth[SPP$species == sp], n = nrow(r$nc),
    best = best_lab,
    lof = f2(tb$LOF_p[tb$model == "Linear"], 2),
    prof = f2(max(r$profile$dnll), 2),
    w = row_ci(r$lin$w, b[, "w"]),
    level = row_ci(r$lin$level, b[, "level"], 0),
    choice = row_ci(100 * expm1(r$dev_obs), 100 * expm1(b[, "dev"]), 0)) }))
write_table(T6, "Table6_larvae", dir = "tables/main", escape = FALSE,
  header = c("Species", "\\makecell{Diet\\\\breadth}", "$n$", "\\makecell{Lowest AICc\\\\(margin)}", "LOF $p$", "\\makecell{Profile range\\\\($\\Delta$nLL)}",
             "$w=a_P/a_C$", "\\makecell{Level $k/a_C$\\\\(mg)}", "\\makecell{Choice\\\\deviation (\\%)}"),
  caption = paste0("ING applied to two specialist--generalist pairs of caterpillars studied at a single food concentration ",
    "(the \\textit{Heliothis} hybrid: Appendix~S3, Table~S3.3). ",
    "Linear: the rule $a_PP+a_CC=k$, in which a single weighted sum of protein and carbohydrate is regulated (the limit ",
    "$\\omega\\to1$ of the anisotropic metric). Lowest AICc: the best model and its margin, in AICc units, over the next model ",
    "(Table~S14). LOF $p$: lack-of-fit test of the linear rule against one mean per food. Profile range: ",
    "range of the profile negative log-likelihood of the anisotropic model over target directions from P:C 5.7:1 to 1:5.7 ",
    "(a range near zero means that the direction of the target is not identified). $w$: weight of protein relative to ",
    "carbohydrate in the regulated sum; level: the regulated sum in carbohydrate equivalents (cumulative intake to the end of ",
    "feeding). Choice deviation: excess of the regulated sum in the choice treatment over the no-choice level (0: on the ",
    "no-choice line). $n$: larvae in the no-choice treatments (12 per species in the choice treatment). Brackets: 95\\% ",
    "bootstrap intervals (larvae resampled within treatments, $B=", B, "$)."), label = "tab:6")
TS14 <- do.call(rbind, lapply(SPP_all$species, function(sp) {
  t <- res[[sp]]$table; t <- t[order(t$AICc), ]
  data.frame(Species = sp_lab(sp), Model = t$model, k = t$k, dAICc = f2(t$dAICc), LR = f2(t$LOF_LR), df = t$LOF_df,
             p = ifelse(is.na(t$LOF_p), "--", ifelse(t$LOF_p < 1e-4, "$<10^{-4}$", f2(t$LOF_p, 3))),
             LORO = ifelse(is.na(t$LORO_pct), "--", f2(t$LORO_pct, 1))) }))
TS14$Species[duplicated(TS14$Species)] <- ""
write_table(TS14, "TableS14_larvae_models", dir = "tables/supplementary", escape = FALSE,
  header = c("Species", "Model", "$k$", "$\\Delta$AICc", "LOF LR", "LOF df", "LOF $p$", "LORO (\\%)"),
  caption = paste("Model comparison for the caterpillar no-choice arrays, ordered by AICc within species. $k$: parameters,",
    "including the residual SD. Linear: the regulated weighted sum $a_PP+a_CC=k$ (anisotropic limit $\\omega\\to1$).",
    "LOF: likelihood-ratio lack-of-fit test against the saturated model (one mean per food). LORO: root-mean-square",
    "error of the predicted mean intake of a held-out food after refitting to the other four.",
    "Model codes as in Table~1."), label = "tab:S14")
cat("Tables 6 and S14 written\n")
