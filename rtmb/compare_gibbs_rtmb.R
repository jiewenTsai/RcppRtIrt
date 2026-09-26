# compare_gibbs_rtmb.R — Gibbs (smi_gibbs.R) vs RTMB (smi_rtmb.R) on the same simulated datasets:
# targets (gap by G, in the model; gap by Z, not in the model), eta grid, risk-rule choice, time.
if (file.exists("../smi/smi_gibbs.R")) setwd("..")
source("pxda/probit_da_gibbs.R"); source("smi/smi_gibbs.R"); source("smi/sim_smi.R")
source("rtmb/smi_rtmb.R")
etas <- c(0, .1, .25, .5, .75, 1)
risk_pick <- function(m, V) { D <- m - m[1]; which.min(D^2 - V[1] + 2 * V) }
rows <- list()
for (sc in list(c(shift = 0, rep = 1), c(shift = 0, rep = 2), c(shift = .5, rep = 1), c(shift = .5, rep = 2))) {
  dat <- sim_gz(seed = 8000 + sc[["rep"]], shift_Z = sc[["shift"]])
  X <- cbind(G = dat$G - mean(dat$G))
  W <- cbind(G = w_gap(dat$G), Z = w_gap(dat$Z))
  truth <- drop(t(W) %*% dat$theta)
  for (e in etas) {
    t0 <- proc.time()[3]
    fg <- smi_rtirt(dat$Y, dat$logT, eta = e, n_iter = 3000, n_burn = 500, seed = 1,
                    temper = "marginal", Xth = X, Xtau = if (e > 0) X else NULL)
    tg <- proc.time()[3] - t0
    dg <- fg$theta %*% W
    t0 <- proc.time()[3]
    fr <- smi_rtmb(dat$Y, dat$logT, eta = e, X = X, W = W)
    tr <- proc.time()[3] - t0
    for (k in 1:2) rows[[length(rows) + 1]] <- data.frame(shift = sc[["shift"]], rep = sc[["rep"]],
      target = colnames(W)[k], eta = e, truth = truth[k],
      gibbs = mean(dg[, k]), gibbs_sd = sd(dg[, k]), rtmb = fr$est[k], rtmb_se = fr$se[k],
      t_gibbs = tg, t_rtmb = tr, conv = fr$conv)
  }
}
res <- do.call(rbind, rows)
saveRDS(res, "rtmb/compare_gibbs_rtmb_results.rds")
options(width = 200)
print(format(res, digits = 3), row.names = FALSE)
cat("\nrisk-rule eta: Gibbs vs RTMB\n")
pk <- do.call(rbind, lapply(split(res, list(res$shift, res$rep, res$target), drop = TRUE), function(r) {
  r <- r[order(r$eta), ]
  data.frame(shift = r$shift[1], rep = r$rep[1], target = r$target[1],
             eta_gibbs = r$eta[risk_pick(r$gibbs, r$gibbs_sd^2)],
             eta_rtmb = r$eta[risk_pick(r$rtmb, r$rtmb_se^2)])
}))
print(pk, row.names = FALSE)
cat("\nmax |gibbs - rtmb| estimate:", round(max(abs(res$gibbs - res$rtmb)), 4),
    "; SE ratio rtmb/gibbs (median):", round(median(res$rtmb_se / res$gibbs_sd), 3), "\n")
cat("mean time per fit: gibbs", round(mean(res$t_gibbs), 2), "s, rtmb", round(mean(res$t_rtmb), 2), "s\n")
