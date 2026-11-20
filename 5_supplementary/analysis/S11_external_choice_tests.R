# Appendix S3: external choice experiments
source("R/ing_model.R"); source("R/table_utils.R"); source("R/fig_theme.R")
fit <- readRDS("out/fits_nochoice.rds")$fits[[ING_BEST]]; p0 <- fit$p
CB <- readRDS("out/cv_boot.rds"); BP <- boot_params(CB, fit)
NB <- as.integer(Sys.getenv("ING_BOOT_EXT", "50")); set.seed(41); BP <- BP[sample(length(BP), min(NB, length(BP)))]

lee <- load_choice(6)
L <- data.frame(study = "Lee et al. 2008", sex = "F", group = paste0("Y", lee$yeast_gL, " S", lee$sucrose_gL),
                FA = lee$yeast_gL, FS = lee$sucrose_gL, cA_P = 0.45, cA_C = 0.24, wP = W_FOOD[1],
                P = lee$P_ug_day, C = lee$C_ug_day, sucrose_food = TRUE, cB_P = 0, cB_C = 1)
je <- read.csv("data/jensen2015_choice_digitised.csv")
rA <- sapply(strsplit(je$foodA_ratio, ":"), function(v) as.numeric(v[1]) / as.numeric(v[2]))
Jd <- data.frame(study = "Jensen et al. 2015", sex = je$sex,
                 group = sprintf("%s@%d vs 0:1@%d", je$foodA_ratio, je$foodA_F_gL, je$sucrose_F_gL),
                 FA = je$foodA_F_gL, FS = je$sucrose_F_gL, cA_P = rA / (1 + rA), cA_C = 1 / (1 + rA), wP = 1,
                 P = je$P_ug_15d / 15, C = je$C_ug_15d / 15, sucrose_food = TRUE, cB_P = 0, cB_C = 1)
mw <- read.csv("data/raw/morimoto_wigby2016_choice_raw.csv"); mw <- mw[mw$Day == 4 & mw$ugP > 0 & mw$ugC > 0, ]
FMW <- c(A = 45, B = 90, C = 180)
mwa <- aggregate(cbind(lP = log(ugP), lC = log(ugC)) ~ Conc, mw, mean)
M <- data.frame(study = "Morimoto & Wigby 2016", sex = "M", group = paste0("1:16 vs 3:1 @", FMW[mwa$Conc]),
                FA = FMW[mwa$Conc], FS = FMW[mwa$Conc], cA_P = 3 / 4, cA_C = 1 / 4, wP = 1,
                P = exp(mwa$lP) / 4, C = exp(mwa$lC) / 4, sucrose_food = FALSE, cB_P = 1 / 17, cB_C = 16 / 17)
D <- rbind(L, Jd, M); D$set <- paste(D$study, D$sex)

solve_pair <- function(p, d, extra = NULL) {
  A <- cbind(d$FA * c(d$cA_P, d$cA_C), d$FS * c(d$cB_P, d$cB_C))
  VM <- vmax_of(p)
  Jm <- function(G) { V <- rowSums(G)
    loss_R(G %*% t(A), p) + p$l1 * V + p$l2 * V^2 + 1e8 * pmax(V - VM, 0)^2 + (if (is.null(extra)) 0 else extra(G, d)) }
  ub <- rep(min(VM, 6), 2); G <- as.matrix(expand.grid(seq(0, ub[1], length.out = 90), seq(0, ub[2], length.out = 90)))
  G <- G[rowSums(G) <= VM, , drop = FALSE]; Jg <- Jm(G); best <- NULL
  for (k in order(Jg)[1:5]) { f <- optim(G[k, ], function(v) Jm(matrix(v, 1)), method = "L-BFGS-B", lower = c(0, 0), upper = ub)
    if (is.null(best) || f$value < best$value) best <- f }
  v <- best$par; if (sum(v) > VM) v <- v * VM / sum(v)
  Jl <- function(s) Jm(cbind(s, VM - s)); ss <- seq(0, VM, length.out = 801); k <- which.min(Jl(ss))
  o <- optimize(Jl, c(ss[max(k - 1, 1)], ss[min(k + 1, 801)])); if (o$objective < Jm(matrix(v, 1))) v <- c(o$minimum, VM - o$minimum)
  x <- A %*% v; c(P = x[[1]], C = x[[2]], VA = v[1], VB = v[2]) }
pred_set <- function(p, d, rho, extra = NULL) { q <- p; q$mu <- rho * p$mu
  t(sapply(seq_len(nrow(d)), function(i) { W_FOOD <<- c(d$wP[i], 1); solve_pair(q, d[i, ], extra) })) }
sse <- function(pr, d) sum((log(pr[, "P"]) - log(d$P))^2 + (log(pr[, "C"]) - log(d$C))^2)
fit_rho <- function(p, d, extra = NULL) optimize(function(lr) sse(pred_set(p, d, exp(lr), extra), d), log(c(.15, 3)))
W0 <- W_FOOD

out <- list()
for (s in unique(D$set)) { d <- D[D$set == s, ]
  o <- fit_rho(p0, d); pr <- pred_set(p0, d, exp(o$minimum))
  bt <- sapply(BP, function(pb) { ob <- fit_rho(pb, d); pp <- pred_set(pb, d, exp(ob$minimum)); log(pp[, "C"] / pp[, "P"]) })
  bt <- matrix(bt, nrow = nrow(d))
  out[[s]] <- data.frame(d, rho = exp(o$minimum), predP = pr[, "P"], predC = pr[, "C"],
                         obs_CP = d$C / d$P, pred_CP = pr[, "C"] / pr[, "P"],
                         lo_CP = exp(apply(bt, 1, quantile, .025)), hi_CP = exp(apply(bt, 1, quantile, .975)))
  cat(sprintf("%-28s rho=%.2f  obs C:P %s | pred %s\n", s, exp(o$minimum), paste(round(d$C / d$P, 2), collapse = " "),
              paste(round(pr[, "C"] / pr[, "P"], 2), collapse = " "))) }
W_FOOD <- W0
OUT <- do.call(rbind, out); rownames(OUT) <- NULL

ext_suc <- function(g) function(G, d) -g * G[, 2] * d$sucrose_food
ext_cap <- function(vc) function(G, d) 1e8 * pmax(G[, 1] - vc, 0)^2 * (d$cA_P > d$cB_P)
hyp <- list()
for (s in unique(D$set)) { d <- D[D$set == s, ]
  n <- 2 * nrow(d)
  base <- fit_rho(p0, d)
  fs <- optim(c(log(.5), log(.2)), function(q) sse(pred_set(p0, d, exp(q[1]), ext_suc(exp(q[2]))), d), control = list(maxit = 300))
  fc <- optim(c(log(.5), log(.4)), function(q) sse(pred_set(p0, d, exp(q[1]), ext_cap(exp(q[2]))), d), control = list(maxit = 300))
  aicc <- function(S, k) n * log(S / n) + 2 * k + 2 * k * (k + 1) / max(n - k - 1, 1)
  prs <- pred_set(p0, d, exp(fs$par[1]), ext_suc(exp(fs$par[2]))); prc <- pred_set(p0, d, exp(fc$par[1]), ext_cap(exp(fc$par[2])))
  hyp[[s]] <- data.frame(set = s, n_trt = nrow(d),
    AICc_base = aicc(base$objective, 1), AICc_suc = aicc(fs$value, 2), AICc_cap = aicc(fc$value, 2),
    gamma = exp(fs$par[2]), VAmax = exp(fc$par[2]),
    rmse_base = sqrt(base$objective / n), rmse_suc = sqrt(fs$value / n), rmse_cap = sqrt(fc$value / n),
    CP_suc = paste(round(prs[, "C"] / prs[, "P"], 2), collapse = " "), CP_cap = paste(round(prc[, "C"] / prc[, "P"], 2), collapse = " "))
  print(hyp[[s]]) }
W_FOOD <- W0
HYP <- do.call(rbind, hyp); rownames(HYP) <- NULL
saveRDS(list(pred = OUT, hyp = HYP, data = D), "out/external_choice.rds")
cat("S11 done\n")
