# Supplementary tables
source("R/ing_model.R"); source("R/table_utils.R"); wrap <- wrap_axis
D <- "tables/supplementary"
fn <- function(x, d = 2) ifelse(is.na(x), "--", formatC(as.numeric(x), format = "f", digits = d))
neg <- function(x) gsub("(^|\\[|, |; )-([0-9])", "\\1$-$\\2", x)
TRL <- c(lifespan = "Lifespan", lifetimeegg = "Lifetime eggs", dailyeggs = "Egg-production rate")
trl <- function(x) { for (k in names(TRL)) x <- gsub(k, TRL[[k]], x, fixed = TRUE); x }
KEY <- "Model codes as in Table~1."
raw_csv <- function(f, key) { r <- file.path("out", paste0("raw_", basename(f)))
  if (file.exists(f) && key %in% names(read.csv(f, nrows = 1, check.names = FALSE))) file.copy(f, r, overwrite = TRUE)
  if (file.exists(r)) read.csv(r) else NULL }

j <- read.csv("out/table_joint_shared_vs_separable.csv"); j <- j[order(j$dAICc), ]
T <- data.frame(Model = j$model, a1 = fn(j$AICc_shared), a2 = fn(j$AICc_separable), d = fn(j$dAICc),
                s1 = fn(j$sig_c_shared), s2 = fn(j$sig_c_sep), cp = fn(j$C_per_P), om = fn(j$omega),
                eP = fn(j$exP), eC = fn(j$exC), rho = fn(j$rho))
write_table(T, "TableS6_joint_models", escape = FALSE, label = "tab:S6", dir = D,
  sideways = TRUE, header = c("Model", "\\makecell{AICc\\\\(shared)}", "\\makecell{AICc\\\\(per food)}", "$\\Delta$AICc", "\\makecell{$\\sigma_c$\\\\(shared)}", "\\makecell{$\\sigma_c$\\\\(per food)}",
             "Target C:P", "$\\omega$", "\\makecell{Protein\\\\ex:def}", "\\makecell{Carbohydrate\\\\ex:def}", "$\\rho$"),
  caption = paste("Joint fits to the no-choice and choice data with common regulation and an appetite scale $\\rho$,",
    "with the volume cost and capacity shared between the two foods or applied to each food separately (per food).",
    "$\\Delta$AICc: relative to the best joint model (per-food form). $\\sigma_c$: among-treatment SD of the choice log",
    "treatment-mean intakes, in addition to their sampling error from the individual flies. The last four models (linear volume",
    "cost) are those that predicted the choice ratios best out of sample (Table~3); models without a quadratic cost have the",
    "same fit under both forms.",
    "Parameter estimates are those of the per-food fit; ex:def, excess:deficit weight.", KEY))

sp <- read.csv("out/table_trait_separation.csv")
T <- data.frame(Comparison = trl(sub(" vs ", " vs.\\ ", sp$comparison)), t1 = fn(sp$theta1, 1), t2 = fn(sp$theta2, 1),
                sep = sprintf("%s [%s, %s]", fn(sp$separation_deg, 1), fn(sp$lo95, 1), fn(sp$hi95, 1)),
                flag = { f <- tolower(sub(":.*", "", sp$check_flag)); paste0(toupper(substr(f, 1, 1)), substring(f, 2)) },
                cov = sprintf("%s [%s, %s]", fn(sp$check_coverage, 3), fn(sp$coverage_lo, 3), fn(sp$coverage_hi, 3)))
write_table(T, "TableS4_trait_separation", escape = FALSE, label = "tab:S4", dir = D,
  header = c("Comparison", "$\\theta_1$ ($^\\circ$)", "$\\theta_2$ ($^\\circ$)", "Separation ($^\\circ$) [95\\% CI]",
             "\\texttt{nutrigon\\_check}", "Coverage [95\\% CI]"),
  caption = paste("Angular separation between trait optima (full quadratic surfaces), estimated with",
    "\\texttt{nutrigon\\_compare} (paired bootstrap of flies, $B=500$). \\texttt{nutrigon\\_check}: reliability of the",
    "interval judged against a validation grid simulated on the Lee design (28 design points, their replication and",
    "intake scatter) at the signal-to-noise ratio of the two surfaces (\\texttt{nutrigon\\_build\\_grid}; 100",
    "configurations, 20 simulated data sets each, $B=50$); coverage: coverage of nominal 95\\% intervals pooled over the",
    "10 nearest configurations."))

ad <- read.csv("out/table_adaptive_test.csv")
if (file.exists("out/table_eggrate_adaptive.csv")) { e <- read.csv("out/table_eggrate_adaptive.csv")
  ad <- rbind(ad, data.frame(objective = e$objective, k = e$k, AICc = e$AICc, log10_lambda2_best = NA,
                             dnll_if_trait_dominates = e$dnll_trait_dominated, LORO_pct = NA)) }
obj <- trl(sub("^maximise ", "Maximise ", ad$objective)); obj <- sub("dailyeggs", "egg-production rate", obj)
obj <- sub("Maximise log egg rate", "Maximise log egg-production rate", obj); obj <- sub("^ING ", "ING, ", obj)
T <- data.frame(Objective = obj, k = ad$k, AICc = fn(ad$AICc), l2 = fn(ad$log10_lambda2_best, 1),
                d = fn(ad$dnll_if_trait_dominates), loro = fn(ad$LORO_pct, 1))
write_table(T, "TableS11_adaptive_test", escape = FALSE, label = "tab:S11", dir = D,
  header = c("Objective", "$k$", "AICc", "$\\log_{10}\\lambda_2$ at best", "$\\Delta$nLL, trait dominant", "LORO (\\%)"),
  caption = paste("Is feeding organised to maximise a fitness component? Models in which intake maximises a",
    "\\texttt{nutrigon} trait surface net of the ING volume costs (Methods, \\textit{Is feeding organised to maximise a fitness component?}), compared with ING.",
    "$\\log_{10}\\lambda_2$ at best: weight of the volume cost at the best fit (6, the upper end of the profile, means",
    "constant-volume feeding, in which the trait has no influence). $\\Delta$nLL, trait dominant: loss of fit when the",
    "trait dominates the objective ($\\lambda_2=10^{-2}$). LORO: leave-one-rail-out prediction error."))

if (file.exists("out/table_eggrate_hessian.csv")) { h <- read.csv("out/table_eggrate_hessian.csv")
  sp2 <- function(x) { a <- strsplit(x, "; boot ")[[1]]; if (length(a) < 2) return(c(a, "--")); a }
  T <- data.frame(Model = h$ING_model, aI = neg(h$axis_ING_deg), aE = neg(h$axis_eggrate_deg),
    d = sapply(h$axis_difference_deg, function(x) sp2(x)[1]), db = sapply(h$axis_difference_deg, function(x) sp2(x)[2]),
    oI = sapply(h$omega_ING, function(x) sp2(x)[1]), oE = sapply(h$omega_eggrate, function(x) sp2(x)[1]),
    oEb = sapply(h$omega_eggrate, function(x) sp2(x)[2]))
  write_table(T, "TableS13_curvature_matching", escape = FALSE, label = "tab:S13", dir = D,
    header = c("ING model", "\\makecell{Axis, ING\\\\($^\\circ$)}", "\\makecell{Axis, egg\\\\rate ($^\\circ$)}", "\\makecell{Axis\\\\difference ($^\\circ$)}",
               "\\makecell{Bootstrap\\\\median [95\\%]}", "$\\omega$, ING", "\\makecell{$\\omega$, egg\\\\rate}", "\\makecell{Bootstrap\\\\median [95\\%]}"),
    caption = paste("Curvature matching between regulation and performance: principal axis and anisotropy of the ING",
      "regulatory metric and of the negative Hessian of log egg-production rate (the only surface with an interior",
      "maximum), both in proportional coordinates around the ING target and normalised to unit trace. Bootstrap: joint",
      "bootstrap in which the ING model and the surface were refitted to the same resamples (replicates with an",
      "interior peak).")) }

if (file.exists("out/table_simulations.csv")) { s <- read.csv("out/table_simulations.csv", check.names = FALSE)
  SM <- readRDS("out/simulations.rds")
  selc <- grep("^pct_sel_", names(s), value = TRUE)
  sel <- sapply(seq_len(nrow(s)), function(i) { g <- s$generating[i]; col <- paste0("pct_sel_", g)
    if (col %in% names(s)) s[i, col] else NA })
  best_alt <- sapply(seq_len(nrow(s)), function(i) { v <- unlist(s[i, selc]); names(v) <- sub("pct_sel_", "", selc)
    v <- v[names(v) != s$generating[i]]; sprintf("%s (%d)", names(which.max(v)), max(v)) })
  T <- data.frame(Scenario = s$scenario, Gen = s$generating, sel = sel, alt = best_alt,
                  cp = sprintf("%s; %s", fn(s$true_CP), s$est_CP), om = sprintf("%s; %s", fn(s$true_omega), s$est_omega),
                  ex = sprintf("%s; %s", fn(s$true_exP), s$est_exP), l2 = sprintf("%s; %s", fn(s$true_l2, 3), s$est_l2))
  write_table(T, "TableS5_simulations", escape = FALSE, label = "tab:S5", dir = D,
    sideways = TRUE, header = c("Scenario", "\\makecell{Generating\\\\model}", "\\makecell{Selected\\\\(\\%)}", "\\makecell{Most frequent\\\\alternative (\\%)}",
               "\\makecell{Target C:P, true;\\\\estimated}", "\\makecell{$\\omega$, true;\\\\estimated}", "\\makecell{Protein ex:def,\\\\true; estimated}",
               "\\makecell{$\\lambda_2$, true;\\\\estimated}"),
    caption = paste("Simulation study on the design of \\citet{Lee2008} (100 data sets per scenario; S1--S7 as in",
      "the Methods). Selected: percentage of data sets in which AICc selected the generating model. Estimates of the",
      "generating model: median [95\\% range]. ex:def, excess:deficit weight.", KEY)) }

if (file.exists("out/table_trait_level.csv")) { o <- read.csv("out/table_trait_level.csv")
  T <- data.frame(Level = sub("individual flies", "Individual flies", sub("treatment means", "Treatment means", o$level)),
                  Trait = TRL[o$trait], P = fn(o$P_opt, 1), C = fn(o$C_opt, 1), th = fn(o$theta_deg, 1), r = o$ratio,
                  cf = ifelse(o$closed_form, "yes", "no"))
  write_table(T, "TableS7_optima_by_level", escape = FALSE, label = "tab:S7", dir = D,
    header = c("Fitting level", "Trait", "$P_{\\mathrm{opt}}$", "$C_{\\mathrm{opt}}$", "$\\theta$ ($^\\circ$)", "P:C", "Interior maximum"),
    caption = paste("Optima of full quadratic surfaces fitted to individual flies or to the 28 treatment means, the",
      "level at which \\citet{Lee2008} reported optima of about 1:16 (lifespan), 1:4 (lifetime eggs) and 1:2 (egg rate).")) }

if (file.exists("out/table_trait_domain.csv")) { o <- read.csv("out/table_trait_domain.csv")
  T <- data.frame(Trait = TRL[o$trait], Domain = gsub("%", "\\\\%", o$domain), Surface = sub("thinplate", "thin-plate", o$surface),
                  P = fn(o$P, 1), C = fn(o$C, 1), r = sub("Inf", "$\\\\infty$", o$ratio), e = fn(o$dist_to_edge, 2))
  T$Trait[duplicated(T$Trait)] <- ""
  write_table(T, "TableS8_optima_by_domain", escape = FALSE, label = "tab:S8", dir = D,
    header = c("Trait", "Domain", "Surface", "$P_{\\mathrm{opt}}$", "$C_{\\mathrm{opt}}$", "P:C", "Distance to edge"),
    caption = paste("Sensitivity of optimum location to surface type (full quadratic or thin-plate regression spline)",
      "and to the definition of the intake domain (convex hull of treatment means, of individuals after removing the",
      "outermost 5\\% by hull peeling, or of all individuals). Distance to edge in \\textmu g\\,day$^{-1}$; every",
      "optimum lies on or next to the edge of its domain.")) }

if (file.exists("out/table_trait_forward_gam.csv")) { o <- read.csv("out/table_trait_forward_gam.csv")
  T <- data.frame(Trait = TRL[o$trait], Surface = sub("thinplate", "thin-plate", o$surface), nll = fn(o$best_nll),
                  AICc = fn(o$AICc), l2 = fn(o$log10_l2_at_best, 1), d = fn(o$nll_W_dominant))
  write_table(T, "TableS12_adaptive_thinplate", escape = FALSE, label = "tab:S12", dir = D,
    header = c("Trait", "Surface", "nLL at best", "AICc", "$\\log_{10}\\lambda_2$ at best", "nLL, trait dominant"),
    caption = paste("Trait-maximisation test (Table~S11) repeated with thin-plate surfaces, evaluated on a fine grid",
      "along each rail. The best ING model has AICc $=-44.24$.")) }

o <- raw_csv(file.path(D, "TableS3_choice_mismatch.csv"), "Extra_parameters")
if (!is.null(o)) {
  sc <- o$Scenario; sc <- gsub("%", "\\\\%", sc); sc <- sub("^Baseline: ", "Baseline, ", sc); sc <- sub("gamma", "$\\\\gamma$", sc)
  T <- data.frame(Scenario = sc, ex = gsub("gamma", "$\\\\gamma$", o$Extra_parameters), rho = fn(o$rho), e = fn(o$Ratio_RMSE),
                  r = fn(o$Ratio_r), p = fn(o$CP_S180_predicted), ob = fn(o$CP_S180_observed), a = fn(o$AICc_choice))
  write_table(T, "TableS3_choice_mismatch", escape = FALSE, label = "tab:S3", dir = D,
    sideways = TRUE, header = c("Scenario", "\\makecell{Extra\\\\parameters}", "$\\rho$", "RMSE", "$r$", "\\makecell{C:P$_{180}$\\\\predicted}", "\\makecell{C:P$_{180}$\\\\observed}", "\\makecell{AICc\\\\(choice)}"),
    caption = paste("Sensitivity of the choice prediction of the best model (shared form) to the assumed yeast",
      "composition (A) and to a sucrose palatability reward estimated from the choice data (C). The no-choice fit",
      "(AICc $=-44.24$) is unchanged in all scenarios. RMSE and $r$: observed against predicted log(C:P) across the",
      "nine treatments; C:P$_{180}$: mean over the three yeast concentrations at sucrose 180~g\\,L$^{-1}$. AICc",
      "(choice): Gaussian AICc of the 18 choice log intakes given the no-choice parameters, counting $\\rho$, the",
      "residual SD and $\\gamma$; not comparable across yeast compositions, which change the observed intakes.")) }

o <- raw_csv(file.path(D, "TableS2_extended_models.csv"), "LORO_no01_pct")
if (!is.null(o)) {
  big <- function(x) ifelse(as.numeric(x) > 1000, "$>10^3$", fn(x, 1))
  T <- data.frame(Model = o$model, k = o$k, a = fn(o$AICc), d = fn(o$dAICc), l1 = big(o$LORO_pct), l2 = fn(o$LORO_no01_pct, 1),
                  lc = big(o$LOCO_pct), cp = fn(o$C_per_P), b = fn(o$beta), g = fn(o$gamma), kp = fn(o$k_pal),
                  r = fn(o$ch_ratio_r), c = fn(o$ch_CP_S180))
  write_table(T, "TableS2_extended_models", escape = FALSE, label = "tab:S2", dir = D,
    sideways = TRUE, header = c("Model", "$k$", "AICc", "$\\Delta$AICc", "\\makecell{LORO\\\\(\\%)}", "\\makecell{LORO$_{-0:1}$\\\\(\\%)}", "\\makecell{LOCO\\\\(\\%)}", "\\makecell{Target\\\\C:P}",
               "$\\beta$", "$\\gamma$", "$k_{\\mathrm{pal}}$", "\\makecell{Choice\\\\$r$}", "C:P$_{180}$"),
    caption = paste("Two non-linear extensions, estimated from the no-choice array: a Box--Cox regulatory scale (BC;",
      "$\\beta=1$ is the proportional scale of the main model) and a sucrose palatability reward (Pal;",
      "$\\gamma(S/180)^{k_{\\mathrm{pal}}}$ per \\textmu L of food with sucrose concentration $S$). LORO$_{-0:1}$:",
      "leave-one-rail-out error excluding the protein-free rail, on which Box--Cox models without capacity extrapolate",
      "to zero intake. Choice $r$ and C:P$_{180}$: out-of-sample prediction of the choice experiment (observed C:P$_{180}$",
      "3.60).", KEY)) }

o <- raw_csv(file.path(D, "TableS9_optima_regions.csv"), "CR_on_hull_edge_pct")
if (!is.null(o)) {
  T <- data.frame(Trait = o$trait, Surface = ifelse(o$surface == "lm", "quadratic", "quasi-Poisson"),
                  opt = sprintf("(%s, %s)", fn(o$P_opt, 0), fn(o$C_opt, 0)), r = o$ratio,
                  cr = sub(" to ", "--", o$CR_ratio_range), edge = fn(o$CR_on_hull_edge_pct, 0), tin = o$target_in_CR,
                  dist = o$dist_target, OR = o$OptimaRegion_opt, ORin = fn(o$OptimaRegion_in_hull_pct, 0))
  write_table(T, "TableS9_optima_regions", escape = FALSE, label = "tab:S9", dir = D,
    sideways = TRUE, header = c("Trait", "Surface", "\\makecell{Optimum\\\\$(P,C)$}", "P:C", "\\makecell{95\\% CR,\\\\P:C range}", "\\makecell{CR on\\\\edge (\\%)}", "\\makecell{Target\\\\in CR}",
               "\\makecell{Distance\\\\to target}", "\\makecell{\\texttt{OptimaRegion}\\\\optimum}", "\\makecell{Inside\\\\hull (\\%)}"),
    caption = paste("95\\% confidence regions (CR) for the location of each trait optimum, from the procedure of",
      "\\citet{delCastillo2022} constrained to the observed intakes (convex hull of individuals): 2,000 coefficient",
      "vectors drawn from the sampling distribution of each surface, trimmed to the 95\\% confidence region of the",
      "coefficients, and each maximised exactly over the hull. CR on edge: percentage of CR optima on the hull boundary.",
      "Distance to target: distance (\\textmu g\\,day$^{-1}$) between the ING target and the optimum, median [95\\% range]",
      "over bootstrap draws of both. Last two columns: optimum from \\texttt{OptRegionQuad}, whose domain is a box",
      "intersected with an origin-anchored triangle, and the percentage of its CR inside the observed intakes.")) }
o <- raw_csv(file.path(D, "TableS10_optima_distances.csv"), "distance_median")
if (!is.null(o)) {
  T <- data.frame(Comparison = sub(" vs ", " vs.\\ ", o$comparison),
                  d = sprintf("%s [%s, %s]", fn(o$distance_median, 1), fn(o$lo95, 1), fn(o$hi95, 1)), ov = o$ratio_ranges_overlap)
  write_table(T, "TableS10_optima_distances", escape = FALSE, label = "tab:S10", dir = D,
    header = c("Comparison", "Distance between optima (\\textmu g\\,day$^{-1}$)", "P:C ranges overlap"),
    caption = paste("Distances between the optima of pairs of traits (full quadratic surfaces, hull-constrained",
      "confidence regions): median [95\\% range] over pairs of draws, and whether the P:C ranges spanned by the two",
      "confidence regions overlap (regions on the hull edge are segments, so overlap is assessed on the ratio).")) }
cat("supplementary tables written to tables/supplementary\n")
