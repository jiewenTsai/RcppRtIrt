# exp_stage2_eta.R — stage-2 inference on the RT module (gamma, v) across eta, no speed shift.
# Under the cut, theta ~ p(theta | Y) is not conditioned on RT, so regressing log T on theta draws
# attenuates gamma by the reliability of theta (errors-in-variables) and inflates v.
suppressMessages(library(parallel))
if (file.exists("../smi/smi_gibbs.R")) setwd("..")
source("pxda/probit_da_gibbs.R"); source("smi/smi_gibbs.R"); source("smi/sim_smi.R")
Rcpp::sourceCpp("smi/rt_sweeps.cpp")
jobs <- expand.grid(eta = c(0, .25, .5, .75, 1), rep = 1:3)
one <- function(k) {
  jb <- jobs[k, ]
  dat <- sim_smi(500, 10, seed = 400 + jb$rep, gam_A = .5, gam_B = .5, gam_sd = .1, v = .15)
  f <- smi_rtirt(dat$Y, dat$logT, eta = jb$eta, n_iter = 3000, n_burn = 500, K_inner = 1,
                 use_cpp = TRUE, temper = "marginal", seed = jb$rep)
  gpool <- (dat$gA + dat$gB) / 2
  post_var_theta <- mean(apply(f$theta, 2, var))
  c(eta = jb$eta, rep = jb$rep, gbar_hat = mean(f$gamma_stage2), gbar_true = mean(gpool),
    v_hat = mean(f$v_stage2), reliability = 1 - post_var_theta,
    theta_rmse = sqrt(mean((colMeans(f$theta) - dat$theta)^2)))
}
R <- as.data.frame(do.call(rbind, mclapply(seq_len(nrow(jobs)), one, mc.cores = 4)))
R$ratio <- R$gbar_hat / R$gbar_true
agg <- aggregate(. ~ eta, data = R[, setdiff(names(R), "rep")], FUN = mean)
options(width = 200); print(round(agg, 3), row.names = FALSE)
cat("truth: v = 0.15; attenuation predicted under the cut ≈ reliability\n")
saveRDS(R, "smi/exp_stage2_eta_results.rds")
