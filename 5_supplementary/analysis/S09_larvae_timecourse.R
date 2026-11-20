# Caterpillar time course
source("R/ing_model.R"); source("R/fig_theme.R")
L <- readRDS("out/larvae_fits.rds"); tc <- L$timecourse; SPP <- L$SPP
COLSP <- c(generalist = ACC, specialist = ACC2, hybrid = "#7c3aed")
fig("FigS11_larvae_timecourse", 7.2, 3.4, dir = "figures/supplementary", function() {
  par(mfrow = c(1, 2), mar = c(3.6, 3.8, 2, 1), oma = c(1.4, 0, 0, 0))
  for (g in c("Heliothis", "Spodoptera")) {
    sp <- SPP[SPP$genus == g, ]
    plot(NA, xlim = c(1, 8), ylim = c(.3, 3), log = "y", xlab = "Day of the experiment", bty = "l",
         ylab = expression("Relative weight of protein, " * w))
    abline(h = 1, col = GRID, lty = 2)
    for (k in seq_len(nrow(sp))) { z <- tc[tc$species == sp$species[k], ]; cl <- COLSP[sp$diet_breadth[k]]
      x <- z$day + (k - (nrow(sp) + 1) / 2) * .12
      segments(x, z$lo, x, z$hi, col = cl); lines(x, z$w, col = cl, lwd = 1.6); points(x, z$w, pch = 19, col = cl, cex = .8) }
    panel_lab(if (g == "Heliothis") "a" else "b"); mtext(bquote(italic(.(g))), 3, .3, cex = .78)
  }
  legend_below(ncol = 3, pch = 19, lwd = 1.6, col = COLSP, legend = c("generalist (H. virescens, S. littoralis)", "specialist (H. subflexa, S. exempta)",
               "hybrid (H. subflexa x virescens)"))
})
cat("Fig. S11 written\n")
