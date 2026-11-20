# Inverse Nutritional Geometry (ING): code and data

Code and data for *Inverse Nutritional Geometry: inferring latent regulatory processes from nutritional arrays* (J. Morimoto). A single script reproduces every analysis, figure and table.

## Run

```
Rscript run_all.R              # full run (about 2 h on 2 cores)
Rscript run_all.R outputs      # rebuild figures and tables from the saved results in out/
ING_CORES=8 Rscript run_all.R  # use 8 cores
```

Run from this folder. Requirements: R >= 4.1 (analyses used R 4.3.3) with `withr` and `mgcv`; `OptimaRegion` (CRAN) for script S08 only. `nutrigon` 0.1.0 is included in `pkg/` and installed by `run_all.R` if missing. A C compiler is optional: without one, R versions of the solvers give identical results more slowly. Tests: `Rscript tests/test_rail_solver.R`, `Rscript tests/test_choice_solver.R`.

## Folders

| Folder | Content |
|---|---|
| `data/` | data (below) |
| `R/` | core library: model, exact solvers, likelihood, fitting |
| `src/` | C versions of the exact solvers |
| `1_schematic_figures/` ... `6_appendices/` | analysis, figure and table scripts, in run order |
| `out/` | saved results (`.rds`, `.csv`) and `sessionInfo.txt` |
| `figures/`, `tables/` | outputs (PDF/PNG; CSV/LaTeX) |
| `tools/` | Python scripts that built the CSVs from the source workbooks and figures |
| `tests/` | solver tests |

## Data

| File | Content | Source |
|---|---|---|
| `lee2008_nochoice.csv` | No-choice experiment, 1,008 females, 28 diets (7 P:C rails x 45/90/180/360 g/L). `protein_eaten`, `carb_eaten`: mean daily intake, days 0-6 (ug/day); `lifespan` (days), `lifetimeegg`, `dailyeggs`; `Ratio`, `Food` (g/L) | Lee et al. 2008, *PNAS* 105:2498 |
| `lee2008_choice_raw.csv` | Choice experiment, individual flies: yeast vs sucrose solution (45/90/180 g/L each). `P_ug_cum`, `C_ug_cum`: cumulative intake (ug) from day 0 to `day` (3, 6, ..., 24); `do_not_use = 1`: two cells flagged by the author | K. P. Lee (unpublished records of Lee et al. 2008) |
| `raw/lee_heliothis_larvae_raw.csv`, `raw/lee_spodoptera_larvae_raw.csv` | Caterpillars: cumulative `protein`, `carb` (mg) per larva and day on five foods (P:C 35:7 to 7:35) and one choice pair | K. P. Lee; Lee et al. 2002, 2003, 2006 |
| `raw/lee_heliothis_hybrid_raw.csv` | Same, F1 hybrid *H. subflexa* female x *H. virescens* male | K. P. Lee (unpublished) |
| `lee_larvae_endpoints.csv` | Cumulative intake per larva at the last record (written by script 10) | derived |
| `jensen2015_nochoice_means.csv`, `raw/jensen2015_table*.csv` | Mean (SE) daily P and C intake per diet and sex | Jensen et al. 2015, *Aging Cell* 14:605, Tables S1-S2 |
| `jensen2015_choice_digitised.csv` | Choice experiment, 15-day intake (SE), digitised from their Fig. 4 (`tools/digitise_jensen2015_fig4.py`) | Jensen et al. 2015 |
| `raw/morimoto_wigby2016_choice_raw.csv` | Male choice between P:C 1:16 and 3:1 foods at three concentrations (A, B, C = 45, 90, 180 g/L); cumulative intake (ug) per male and day | Morimoto & Wigby 2016, *Sci. Rep.* 6:27673 |

The K. P. Lee data are shared for review of this manuscript. Please contact the authors before reusing them.

## Outputs and the manuscript

File names follow an earlier numbering. Correspondence with the submitted manuscript:

| Manuscript | File |
|---|---|
| Fig. 1 | `figures/main/Fig1_geometry_schematic` |
| Fig. 2 | `figures/main/Fig4_nochoice_fit` |
| Fig. 3 | `figures/main/Fig5_target_identifiability` |
| Fig. 4 | `figures/main/Fig6_choice_prediction` |
| Fig. 5 | `figures/main/Fig8_specialist_generalist` |
| Fig. S14 | `figures/main/Fig2_evolute_construction` |
| Fig. S15 | `figures/main/Fig3_simulation_study` |
| Fig. S16 | `figures/main/Fig7_trait_surfaces` |
| Figs S1-S13 | `figures/supplementary/FigS*` |
| Tables 1, 2 | `tables/main/Table1_*`, `Table2_*` |
| Table 3 | `tables/main/Table6_larvae` |
| Tables S19, S20, S21 | `tables/main/Table3_*`, `Table4_*`, `Table5_*` |
| Tables S1-S18 | `tables/supplementary/TableS*` |
| Appendix S1-S3 figures and tables | `figures/appendix/`, `tables/appendix/` |
