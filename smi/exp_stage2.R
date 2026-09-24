# exp_stage2.R — SMI stage 2 (RT module given theta, full likelihood)
#  1. rt_sweeps_cpp() matches rt_block() in distribution (fixed theta, long runs)
#  2. cost per sweep, R vs C++
#  3. nested MCMC with K inner sweeps (warm start) vs a reference that, for each thinned stage-1
#     theta draw, runs an independent 300-sweep inner chain from a fixed start
if (file.exists("../smi/smi_gibbs.R")) setwd("..")
source("pxda/probit_da_gibbs.R"); source("smi/smi_gibbs.R"); source("smi/sim_smi.R")
Rcpp::sourceCpp("smi/rt_sweeps.cpp")
dat <- sim_smi(500, 10, seed = 401, gam_A = .5, gam_B = .5, gam_sd = .1, v = .15, tau_shift_B = 0)
n <- 500; p <- 10; logT <- dat$logT

## 1-2: R vs C++ at theta fixed at the truth -----------------------------------------------
env <- environment(smi_rtirt)
rt_block_R <- function(th, st, n_sweeps) {
  # re-create rt_block with the closure variables smi_rtirt uses
  P0xg <- diag(c(1 / 100, 1)); m0xg <- c(4, 0)
  out <- matrix(NA, n_sweeps, 2)
  for (k in 1:n_sweeps) {
    w <- 1 / st$s2; prec_tau <- 1 / st$v + sum(w)
    rhs <- drop((matrix(st$xi, n, p, byrow = TRUE) + outer(th, st$gam) - logT) %*% w)
    st$tau <- rhs / prec_tau + rnorm(n) / sqrt(prec_tau)
    W <- cbind(1, th); WtW <- crossprod(W); Wty <- crossprod(W, logT + st$tau)
    for (j in 1:p) {
      Vj <- solve(P0xg + WtW / st$s2[j]); mj <- Vj %*% (P0xg %*% m0xg + Wty[, j] / st$s2[j])
      xg <- mj + t(chol(Vj)) %*% rnorm(2); st$xi[j] <- xg[1]; st$gam[j] <- xg[2]
    }
    E <- logT - (matrix(st$xi, n, p, byrow = TRUE) - st$tau + outer(th, st$gam))
    st$s2 <- 1 / rgamma(p, 1 + n / 2, 1 + colSums(E^2) / 2)
    st$v <- 1 / rgamma(1, 1 + n / 2, 1 + sum(st$tau^2) / 2)
    out[k, ] <- c(mean(st$gam), st$v)
  }
  out
}
st0 <- list(tau = rep(0, n), xi = colMeans(logT), gam = rep(0, p), s2 = rep(.25, p), v = 1)
set.seed(1); t0 <- proc.time()[3]; oR <- rt_block_R(dat$theta, st0, 5000); tR <- proc.time()[3] - t0
set.seed(2); oC <- matrix(NA, 5000, 2); st <- st0; t0 <- proc.time()[3]
for (k in 1:5000) { st <- rt_sweeps_cpp(dat$theta, logT, st$tau, st$xi, st$gam, st$s2, st$v, 1)
  oC[k, ] <- c(mean(st$gam), st$v) }
tC <- proc.time()[3] - t0
set.seed(3); t0 <- proc.time()[3]; invisible(rt_sweeps_cpp(dat$theta, logT, st0$tau, st0$xi, st0$gam, st0$s2, st0$v, 5000))
tC_batch <- proc.time()[3] - t0
qs <- c(.05, .5, .95)
cat("== R vs C++ inner sweeps, theta fixed (5000 sweeps, first 500 dropped)\n")
print(round(rbind(R_gbar = quantile(oR[-(1:500), 1], qs), Cpp_gbar = quantile(oC[-(1:500), 1], qs),
                  R_v = quantile(oR[-(1:500), 2], qs), Cpp_v = quantile(oC[-(1:500), 2], qs)), 4))
cat(sprintf("ms per sweep: R %.3f | C++ one call per sweep %.3f | C++ batched %.3f\n",
            1000 * tR / 5000, 1000 * tC / 5000, 1000 * tC_batch / 5000))

## 3: nested K vs reference, eta = 0 (cut) --------------------------------------------------
res <- list()
for (K in c(1, 5, 20)) {
  t0 <- proc.time()[3]
  f <- smi_rtirt(dat$Y, logT, eta = 0, n_iter = 3000, n_burn = 500, K_inner = K, use_cpp = TRUE, seed = 7)
  res[[paste0("K", K)]] <- list(gbar = rowMeans(f$gamma_stage2), v = f$v_stage2, g1 = f$gamma_stage2[, 1],
                                secs = proc.time()[3] - t0, theta = f$theta)
}
TH <- res$K20$theta[seq(1, 2500, by = 5), ]            # 500 thinned stage-1 draws
set.seed(11); t0 <- proc.time()[3]
ref <- t(apply(TH, 1, function(th) { s <- rt_sweeps_cpp(th, logT, st0$tau, st0$xi, st0$gam, st0$s2, st0$v, 300)
  c(mean(s$gam), s$v, s$gam[1]) }))
t_ref <- proc.time()[3] - t0
cat("\n== stage 2 under the cut (eta = 0): nested K vs reference (fresh 300-sweep chain per theta draw)\n")
tab <- rbind(reference = c(quantile(ref[, 1], qs), quantile(ref[, 2], qs), quantile(ref[, 3], qs), secs = t_ref),
             t(sapply(res, function(r) c(quantile(r$gbar, qs), quantile(r$v, qs), quantile(r$g1, qs), secs = r$secs))))
colnames(tab) <- c(paste0("gbar_", qs), paste0("v_", qs), paste0("g1_", qs), "secs")
print(round(tab, 4))
# sim_smi() draws gamma_j separately for the two groups (same mean), so the pooled truth is their average
gpool <- (dat$gA + dat$gB) / 2
cat(sprintf("truth (pooled over groups): gbar = %.3f, v = 0.15, g1 = %.3f\n", mean(gpool), gpool[1]))
