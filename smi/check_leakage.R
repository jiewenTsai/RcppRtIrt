# check_leakage.R — does lambda(eta) = eta 1'Omega^-1 gamma / (I_irt + eta gamma'Omega^-1 gamma) predict
# the bias of the theta gap along a speed-shifted variable not in the model?
# Omega = diag(sigma^2) + v 11' (marginal RT covariance), I_irt = mean posterior precision of theta
# from the responses and the prior (estimated as 1 / posterior variance of theta at eta = 0).
# Predicted bias for a Z gap = lambda * (tau gap on Z not absorbed by the model) * (-1) (tau lower = slower).
if (file.exists("../smi/smi_gibbs.R")) setwd("..")
source("pxda/probit_da_gibbs.R"); source("smi/smi_gibbs.R"); source("smi/sim_smi.R")
lam <- function(eta, gam, s2, v, I_irt) {
  Oi <- solve(diag(s2) + v)                       # v added to every entry = v 11'
  eta * sum(Oi %*% gam) / (I_irt + eta * drop(t(gam) %*% Oi %*% gam))
}
res <- NULL
for (r in 1:4) {
  dat <- sim_gz(seed = 4000 + r, shift_Z = .5)
  X <- cbind(G = dat$G - mean(dat$G))
  f0 <- smi_rtirt(dat$Y, dat$logT, eta = 0, n_iter = 2000, n_burn = 500, seed = r, Xth = X)
  I_irt <- 1 / mean(apply(f0$theta, 2, var))
  gz0 <- mean(rowMeans(f0$theta[, dat$Z == 1]) - rowMeans(f0$theta[, dat$Z == 0]))
  # tau gap on Z left after regressing the true shift on G (what the conditioning model cannot absorb)
  shift <- -.5 * dat$Z; resid <- resid(lm(shift ~ dat$G))
  dtau <- mean(resid[dat$Z == 1]) - mean(resid[dat$Z == 0])
  for (e in c(.25, .5, 1)) {
    f <- smi_rtirt(dat$Y, dat$logT, eta = e, n_iter = 2000, n_burn = 500, seed = r,
                   temper = "marginal", Xth = X, Xtau = X)
    gz <- mean(rowMeans(f$theta[, dat$Z == 1]) - rowMeans(f$theta[, dat$Z == 0]))
    # plug-in lambda from the fitted RT module (stage-1 posterior means)
    L <- lam(e, colMeans(f$gamma_stage1), f$s2_stage1_mean, mean(f$v_stage1), I_irt)
    res <- rbind(res, data.frame(rep = r, eta = e, lambda = L, dtau = dtau,
                                 predicted = -L * dtau, observed = gz - gz0))
  }
}
print(format(res, digits = 3), row.names = FALSE)
cat("\nmean observed / predicted by eta:\n")
print(aggregate(cbind(predicted, observed) ~ eta, res, mean))
