# Trait optima from treatment means
suppressMessages(library(nutrigon)); data(lee_drosophila)
traits <- c("lifespan", "lifetimeegg", "dailyeggs")
agg <- aggregate(cbind(protein_eaten, carb_eaten, lifespan, lifetimeegg, dailyeggs) ~ Ratio + Food,
                 data = lee_drosophila, FUN = mean)
ind <- readRDS("out/traits_nutrigon.rds")$fits
row <- function(f, level, y) { o <- f$optimum
  data.frame(level = level, trait = y, P_opt = round(o$protein_eaten, 1), C_opt = round(o$carb_eaten, 1),
             theta_deg = round(o$theta, 2), ratio = sprintf("1:%.2f", o$carb_eaten / o$protein_eaten),
             closed_form = !o$used_fallback, elongation = round(f$elongation, 2)) }
tab <- do.call(rbind, c(lapply(traits, function(y) row(ind[[paste(y, "lm")]], "individual flies (n = 973)", y)),
  lapply(traits, function(y) { set.seed(2025)
    row(nutrigon_fit(agg, "protein_eaten", "carb_eaten", y, n_cloud = 4000, n_reps = 10), "treatment means (n = 28)", y) })))
print(tab, row.names = FALSE); write.csv(tab, "out/table_trait_level.csv", row.names = FALSE)
