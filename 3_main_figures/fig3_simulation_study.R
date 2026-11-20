# Fig. S15: simulation study
source("R/ing_model.R"); source("R/fig_theme.R")
F0 <- readRDS("out/fits_nochoice.rds"); dat <- F0$dat; fits <- F0$fits
best <- fits[[ING_BEST]]; symm <- fits[[ING_SYMM]]

{
  SM <- readRDS("out/simulations.rds"); out <- SM$out; gen <- SM$gen; cn <- SM$cands
  truth <- lapply(gen, function(n) fits[[n]]$desc)
  sc <- names(out); labs <- c(S1 = "S1 iso", S2 = "S2 aniso", S3 = "S3 iso+asym", S4 = "S4 aniso+asym",
                              S5 = "S5 indiv. targets", S6 = "S6 4 rails", S7 = "S7 2x noise")[sc]
  fig("Fig3_simulation_study", 7.2, 6.6, function() {
    par(mfrow = c(2, 2), mar = c(6.2, 4, 2, 1))
    par(mar = c(6.2, 4, 2, 1))
    rel <- function(v, s, key) v / truth[[s]][[key]]
    bx <- function(key, ylab, lab, ylim, logy = TRUE, which = sc) {
      L <- lapply(which, function(s) rel(out[[s]][, key], s, key)); names(L) <- labs[which]
      boxplot(L, log = if (logy) "y" else "", las = 2, cex.axis = .7, ylab = ylab, ylim = ylim, outline = FALSE,
              col = "#dbeafe", border = ACC, medcol = INK, whisklty = 1, staplewex = 0, boxwex = .55, frame = FALSE)
      abline(h = 1, col = ACC2, lty = 2); panel_lab(lab) }
    bx("C_per_P", "Estimated / true target C:P", "a", c(.6, 1.6)); mtext("Intake-target ratio", 3, .3, cex = .78)
    M <- sapply(sc, function(s) tabulate(out[[s]][, "sel"], length(cn)) / nrow(out[[s]]))
    par(mar = c(6.2, 7.6, 2, 1))
    image(seq_along(sc), seq_along(cn), t(M), col = colorRampPalette(c("#f8fafc", "#93c5fd", "#1e3a8a"))(30), zlim = c(0, 1),
          axes = FALSE, xlab = "", ylab = "")
    axis(1, seq_along(sc), labs, las = 2, cex.axis = .7, tick = FALSE); axis(2, seq_along(cn), cn, las = 1, cex.axis = .7, tick = FALSE)
    for (i in seq_along(sc)) for (j in seq_along(cn)) text(i, j, round(100 * M[j, i]), cex = .65, col = if (M[j, i] > .55) "white" else INK)
    gi <- match(gen, cn); points(seq_along(sc)[!is.na(gi)], gi[!is.na(gi)], pch = 0, cex = 3.2, col = ACC2, lwd = 1.5)
    panel_lab("b"); mtext("Model selected by AICc (%)", 3, .3, cex = .78)
    par(mar = c(6.2, 4, 2, 1))
    bx("lambda2", expression("Estimated / true " * lambda[2]), "c", c(.6, 1.6)); mtext("Quadratic volume cost", 3, .3, cex = .78)
    L <- lapply(sc, function(s) out[[s]][, "omega"]); names(L) <- labs
    boxplot(L, las = 2, cex.axis = .7, ylab = expression("Anisotropy " * omega), ylim = c(.45, 1), outline = FALSE,
            col = "#dbeafe", border = ACC, medcol = INK, whisklty = 1, staplewex = 0, boxwex = .55, frame = FALSE)
    points(seq_along(sc), sapply(sc, function(s) truth[[s]][["omega"]]), pch = 23, bg = ACC2, col = "white", cex = 1.2)
    panel_lab("d"); mtext("Regulatory anisotropy (diamond = truth)", 3, .3, cex = .78)
  })
}
cat("Figure 3 written\n")
