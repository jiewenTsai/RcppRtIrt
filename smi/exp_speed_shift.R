# exp_speed_shift.R — group B is slower (tau mean shifted) but equally able; model assumes one
# speed distribution. How does eta trade bias in group B against precision?
suppressMessages(library(parallel))
if (file.exists("../smi/smi_gibbs.R")) setwd("..")
source("pxda/probit_da_gibbs.R"); source("smi/smi_gibbs.R"); source("smi/sim_smi.R")
etas <- c(0, .25, .5, .75, 1); shifts <- c(0, .25, .5); reps <- 1:3
jobs <- expand.grid(eta = etas, shift = shifts, rep = reps)
one <- function(k) {
  jb <- jobs[k, ]
  dat <- sim_smi(500, 10, seed = 200 + jb$rep, gam_A = .5, gam_B = .5, gam_sd = .1, v = .15,
                 tau_shift_B = jb$shift)
  fit <- smi_rtirt(dat$Y, dat$logT, eta = jb$eta, n_iter = 3000, n_burn = 500, seed = jb$rep)
  m <- colMeans(fit$theta); sdv <- apply(fit$theta, 2, sd)
  lo <- apply(fit$theta, 2, quantile, .025); hi <- apply(fit$theta, 2, quantile, .975)
  cov <- dat$theta >= lo & dat$theta <= hi; err <- m - dat$theta; B <- dat$grp == 1
  c(eta = jb$eta, shift = jb$shift, rep = jb$rep,
    bias_A = mean(err[!B]), bias_B = mean(err[B]), gap_BA = mean(m[B]) - mean(m[!B]),
    rmse = sqrt(mean(err^2)), rmse_B = sqrt(mean(err[B]^2)),
    cover_A = mean(cov[!B]), cover_B = mean(cov[B]), post_sd = mean(sdv),
    true_gap = mean(dat$theta[B]) - mean(dat$theta[!B]))
}
R <- as.data.frame(do.call(rbind, mclapply(seq_len(nrow(jobs)), one, mc.cores = 4)))
agg <- aggregate(. ~ shift + eta, data = R[, setdiff(names(R), "rep")], FUN = mean)
agg <- agg[order(agg$shift, agg$eta), ]
options(width = 200); print(round(agg, 3), row.names = FALSE)
saveRDS(R, "smi/exp_speed_shift_results.rds")
