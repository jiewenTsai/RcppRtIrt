# exp_temper.R — fine eta grid (likelihood tempering) with and without a speed shift in group B
suppressMessages(library(parallel))
if (file.exists("../smi/smi_gibbs.R")) setwd("..")
source("pxda/probit_da_gibbs.R"); source("smi/smi_gibbs.R"); source("smi/sim_smi.R")
etas <- c(0, .02, .05, .1, .25, .5, 1)
jobs <- expand.grid(eta = etas, shift = c(0, .5), rep = 1:3)
one <- function(k) {
  jb <- jobs[k, ]
  dat <- sim_smi(500, 10, seed = 200 + jb$rep, gam_A = .5, gam_B = .5, gam_sd = .1, v = .15,
                 tau_shift_B = jb$shift)
  fit <- smi_rtirt(dat$Y, dat$logT, eta = jb$eta, n_iter = 3000, n_burn = 500, seed = jb$rep)
  m <- colMeans(fit$theta); err <- m - dat$theta; B <- dat$grp == 1
  c(eta = jb$eta, shift = jb$shift, rep = jb$rep,
    gap_bias = (mean(m[B]) - mean(m[!B])) - (mean(dat$theta[B]) - mean(dat$theta[!B])),
    rmse = sqrt(mean(err^2)), post_sd = mean(apply(fit$theta, 2, sd)),
    max_abs_a = max(abs(fit$a)), max_abs_theta = max(abs(fit$theta)))
}
res <- mclapply(seq_len(nrow(jobs)), one, mc.cores = 4)
R <- as.data.frame(do.call(rbind, res))
agg <- aggregate(cbind(gap_bias, rmse, post_sd, max_abs_a, max_abs_theta) ~ shift + eta,
                 data = R, FUN = mean)
agg <- agg[order(agg$shift, agg$eta), ]
options(width = 200); print(format(agg, digits = 3), row.names = FALSE)
saveRDS(R, "smi/exp_temper_results.rds")
