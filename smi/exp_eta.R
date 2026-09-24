# exp_eta.R — theta recovery by group across eta, correct vs misspecified RT module
suppressMessages(library(parallel))
if (file.exists("../smi/smi_gibbs.R")) setwd("..")
source("pxda/probit_da_gibbs.R"); source("smi/smi_gibbs.R"); source("smi/sim_smi.R")
etas <- c(0, .25, .5, .75, 1)
scen <- list(correct = c(gam_A = .3, gam_B = .3), misspec = c(gam_A = .4, gam_B = -.2))
reps <- 1:3
jobs <- expand.grid(eta = etas, scen = names(scen), rep = reps, stringsAsFactors = FALSE)
one <- function(k) {
  jb <- jobs[k, ]; sc <- scen[[jb$scen]]
  dat <- sim_smi(500, 15, seed = 100 + jb$rep, gam_A = sc["gam_A"], gam_B = sc["gam_B"])
  fit <- smi_rtirt(dat$Y, dat$logT, eta = jb$eta, n_iter = 3000, n_burn = 500, seed = jb$rep)
  m <- colMeans(fit$theta); sdv <- apply(fit$theta, 2, sd)
  lo <- apply(fit$theta, 2, quantile, .025); hi <- apply(fit$theta, 2, quantile, .975)
  cov <- dat$theta >= lo & dat$theta <= hi; err <- m - dat$theta
  B <- dat$grp == 1
  c(eta = jb$eta, rep = jb$rep,
    bias_A = mean(err[!B]), bias_B = mean(err[B]),
    rmse_A = sqrt(mean(err[!B]^2)), rmse_B = sqrt(mean(err[B]^2)),
    cover_A = mean(cov[!B]), cover_B = mean(cov[B]), post_sd = mean(sdv),
    rmse_a = sqrt(mean((colMeans(fit$a) - dat$a)^2)),
    gam_hat = mean(fit$gamma_stage1))
}
res <- mclapply(seq_len(nrow(jobs)), one, mc.cores = 4)
R <- cbind(scen = jobs$scen, as.data.frame(do.call(rbind, res)))
agg <- aggregate(. ~ scen + eta, data = R[, setdiff(names(R), "rep")], FUN = mean)
agg <- agg[order(agg$scen, agg$eta), ]
print(format(agg, digits = 3), row.names = FALSE)
saveRDS(R, "smi/exp_eta_results.rds")
