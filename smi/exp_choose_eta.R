# exp_choose_eta.R — data-driven choice of eta for the group gap Delta = mean(theta_B) - mean(theta_A),
# using the cut posterior (eta = 0) as the unbiased reference.
#   D(eta)  = Delta_hat(eta) - Delta_hat(0)
#   V(eta)  = posterior variance of Delta under eta
#   Hausman: H(eta) = D^2 / (V(0) - V(eta));  rule "hausman" = largest eta with H < 3.84
#   Risk:    R(eta) = D^2 - V(0) + 2 V(eta)  (unbiased for MSE when Var(D) ≈ V(0) - V(eta));
#            rule "risk" = argmin R
suppressMessages(library(parallel))
if (file.exists("../smi/smi_gibbs.R")) setwd("..")
source("pxda/probit_da_gibbs.R"); source("smi/smi_gibbs.R"); source("smi/sim_smi.R")
etas <- c(0, .1, .25, .4, .55, .75, 1)
jobs <- expand.grid(eta = etas, shift = c(0, .25, .5), rep = 1:6)
fit_one <- function(k) {
  jb <- jobs[k, ]
  dat <- sim_smi(500, 10, seed = 300 + jb$rep, gam_A = .5, gam_B = .5, gam_sd = .1, v = .15,
                 tau_shift_B = jb$shift)
  fit <- smi_rtirt(dat$Y, dat$logT, eta = jb$eta, n_iter = 3000, n_burn = 500, seed = jb$rep,
                   temper = "marginal")
  B <- dat$grp == 1
  g <- rowMeans(fit$theta[, B]) - rowMeans(fit$theta[, !B])
  list(eta = jb$eta, shift = jb$shift, rep = jb$rep, gap_mean = mean(g), gap_var = var(g),
       theta_mean = colMeans(fit$theta), true_gap = mean(dat$theta[B]) - mean(dat$theta[!B]),
       theta = dat$theta)
}
fits <- mclapply(seq_len(nrow(jobs)), function(k) tryCatch(fit_one(k), error = function(e) {
  message("job ", k, ": ", conditionMessage(e)); NULL }), mc.cores = 4)
saveRDS(fits, "smi/exp_choose_eta_fits.rds")

rows <- list()
for (sh in unique(jobs$shift)) for (r in unique(jobs$rep)) {
  F <- Filter(function(f) !is.null(f) && f$shift == sh && f$rep == r, fits)
  F <- F[order(sapply(F, `[[`, "eta"))]
  e <- sapply(F, `[[`, "eta"); m <- sapply(F, `[[`, "gap_mean"); V <- sapply(F, `[[`, "gap_var")
  D <- m - m[1]; dV <- pmax(V[1] - V, 1e-12)
  H <- ifelse(e == 0, 0, D^2 / dV)
  R <- D^2 - V[1] + 2 * V
  pick_risk <- which.min(R)
  ok <- which(H < qchisq(.95, 1)); pick_h <- max(ok[ok == seq_along(ok)])   # largest eta before first rejection
  truth <- F[[1]]$true_gap; th <- F[[1]]$theta
  for (rule in c("cut", "full", "risk", "hausman")) {
    idx <- switch(rule, cut = 1, full = length(F), risk = pick_risk, hausman = pick_h)
    rows[[length(rows) + 1]] <- data.frame(shift = sh, rep = r, rule = rule, eta = e[idx],
      gap_err = m[idx] - truth, theta_rmse = sqrt(mean((F[[idx]]$theta_mean - th)^2)))
  }
}
out <- do.call(rbind, rows)
summ <- aggregate(cbind(eta, gap_err, abs_gap_err = abs(gap_err), theta_rmse) ~ shift + rule, data = out, FUN = mean)
summ$rmse_gap <- aggregate(gap_err ~ shift + rule, data = out, FUN = function(x) sqrt(mean(x^2)))$gap_err
summ <- summ[order(summ$shift, match(summ$rule, c("cut", "full", "risk", "hausman"))), ]
options(width = 200); print(format(summ, digits = 3), row.names = FALSE)
cat("\nchosen eta by dataset\n"); print(xtabs(eta ~ shift + rep + rule, data = out)[, , c("risk", "hausman")])
saveRDS(out, "smi/exp_choose_eta_results.rds")
