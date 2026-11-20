# Thin-plate optima: distance to the design edge
source("R/ing_model.R"); suppressMessages(library(mgcv))
d <- read_lee("data/lee2008_nochoice.csv"); d <- d[complete.cases(d[, c("protein_eaten","carb_eaten","lifespan","lifetimeegg","dailyeggs")]), ]
d$P <- d$protein_eaten; d$C <- d$carb_eaten
doms <- readRDS("out/trait_domains.rds"); o <- read.csv("out/table_trait_domain.csv")
segdist <- function(p, a, b) { ab <- b - a; if (sum(ab^2) == 0) return(sqrt(sum((p - a)^2))); t <- max(0, min(1, sum((p - a) * ab) / sum(ab^2))); sqrt(sum((p - a - t * ab)^2)) }
bdist <- function(p, h) { h <- unique(h); n <- nrow(h); min(sapply(1:n, function(i) segdist(p, h[i, ], h[if (i == n) 1 else i + 1, ]))) }
o$dist_to_edge <- round(mapply(function(P, C, dm) bdist(c(P, C), doms[[dm]]), o$P, o$C, o$domain), 2)
print(o, row.names = FALSE)
write.csv(o, "out/table_trait_domain.csv", row.names = FALSE)
