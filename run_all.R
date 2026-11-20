# Reproduces every analysis, figure and table
args <- commandArgs(TRUE); mode <- if (length(args)) args[1] else "all"
dir.create("out/logs", recursive = TRUE, showWarnings = FALSE)

so <- file.path("src", paste0("railopt", .Platform$dynlib.ext))
if (!file.exists(so)) {
  ok <- system2("R", c("CMD", "SHLIB", "src/railopt.c"), stdout = FALSE, stderr = FALSE) == 0
  if (!ok) message("C compiler not available: using the (identical, slower) R solvers.")
}
if (!requireNamespace("nutrigon", quietly = TRUE) || !"design" %in% names(formals(nutrigon::nutrigon_build_grid)))
  install.packages("pkg/nutrigon_0.1.0.tar.gz", repos = NULL, type = "source")
for (p in c("withr", "mgcv")) if (!requireNamespace(p, quietly = TRUE)) stop("Please install R package: ", p)
if (!requireNamespace("OptimaRegion", quietly = TRUE))
  try(install.packages("OptimaRegion", repos = "https://cloud.r-project.org"), silent = TRUE)

run <- function(script, ...) {
  t0 <- Sys.time(); cat(sprintf("  %-58s", script))
  log <- file.path("out/logs", sub("\\.R$", ".log", basename(script)))
  st <- system2("Rscript", c(script, ...), stdout = log, stderr = log)
  cat(sprintf("%s (%.1f min)\n", if (st == 0) "ok" else "FAILED", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  if (st != 0) stop("Stopped at ", script, " -- see ", log, call. = FALSE)
}
stage <- function(title, scripts) { cat(title, "\n"); for (s in scripts) do.call(run, as.list(s)) }

tests      <- c("tests/test_rail_solver.R", "tests/test_choice_solver.R")
schematic  <- c("1_schematic_figures/fig1_geometry_schematic.R", "1_schematic_figures/fig2_evolute_construction.R")
analysis   <- list("2_main_analysis/01_fit_nochoice_models.R",
                   "2_main_analysis/02_crossvalidation_bootstrap.R",
                   c("2_main_analysis/03_simulation_study.R", "100"),
                   "2_main_analysis/04_choice_joint_fits.R",
                   "2_main_analysis/05_target_profile_likelihood.R",
                   "2_main_analysis/06_choice_prediction_and_costs.R",
                   "2_main_analysis/07_trait_surfaces_nutrigon.R",
                   "2_main_analysis/08_curvature_matching.R",
                   "2_main_analysis/09_trait_maximisation_test.R",
                   "2_main_analysis/10_specialist_generalist.R")
supp_an    <- sprintf("5_supplementary/analysis/%s", c("S01_trait_optima_treatment_means.R", "S02_thinplate_surfaces.R",
                "S03_thinplate_optima_by_domain.R", "S04_thinplate_distance_to_edge.R", "S05_trait_maximisation_thinplate.R",
                "S06_choice_mismatch_sensitivity.R", "S07_extended_models.R", "S08_optima_regions.R",
                "S09_larvae_timecourse.R", "S10_jensen2015_nochoice.R", "S11_external_choice_tests.R",
                "S12_appendixS3_outputs.R", "S13_choice_timecourse.R", "S14_heliothis_hybrid.R",
                "S15_misspecification.R", "S16_fitness_competitors.R", "S17_parameter_identifiability.R"))
appendix   <- c("6_appendices/appendixS1_smoothed_evolute.R", "6_appendices/appendixS2_kink_verification.R")
main_figs  <- sprintf("3_main_figures/%s", c("fig3_simulation_study.R", "fig4_nochoice_fit.R",
                "fig5_target_identifiability.R", "fig6_choice_prediction.R", "fig7_trait_surfaces.R",
                "fig8_specialist_generalist.R"))
outputs    <- c(main_figs, "4_main_tables/main_tables.R", "4_main_tables/table6_larvae.R",
                "5_supplementary/supplementary_figures.R", "5_supplementary/supplementary_tables.R",
                "5_supplementary/analysis/S09_larvae_timecourse.R", "5_supplementary/analysis/S12_appendixS3_outputs.R",
                "5_supplementary/analysis/S14_heliothis_hybrid.R")

if (mode == "schematics") { stage("Stage 1: schematic figures", schematic); quit(save = "no") }
if (mode == "all") {
  stage("Stage T: solver tests", tests)
  stage("Stage 1: schematic figures", schematic)
  stage("Stage 2: main analysis", analysis)
  stage("Stage S: supplementary analyses", supp_an)
  stage("Stage A: appendix analyses", appendix)
}
stage("Stages 3-5: figures and tables", outputs)
writeLines(capture.output(sessionInfo()), "out/sessionInfo.txt")
cat("\nDone. figures/{main,supplementary,appendix}  tables/{main,supplementary,appendix}\n")
