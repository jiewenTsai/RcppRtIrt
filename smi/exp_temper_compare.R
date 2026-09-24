# exp_temper_compare.R — likelihood vs marginal tempering on an eta grid, with and without a
# speed shift in group B. info_frac = share of the cut-to-full gain in mean posterior precision of theta.
suppressMessages(library(parallel))
if (file.exists("../smi/smi_gibbs.R")) setwd("..")
source("pxda/probit_da_gibbs.R"); source("smi/smi_gibbs.R"); source("smi/sim_smi.R")
etas <- c(0, .05, .1, .25, .5, .75, 1)
jobs <- expand.grid(eta = etas, temper = c("likelihood", "marginal"), shift = c(0, .5), rep = 1:3,
                    stringsAsFactors = FALSE)
one <- function(k) {
  jb <- jobs[k, ]
  dat <- sim_smi(500, 10, seed = 200 + jb$rep, gam_A = .5, gam_B = .5, gam_sd = .1, v = .15,
                 tau_shift_B = jb$shift)
  fit <- smi_rtirt(dat$Y, dat$logT, eta = jb$eta, n_iter = 3000, n_burn = 500, seed = jb$rep,
                   temper = jb$temper)
  m <- colMeans(fit$theta); err <- m - dat$theta; B <- dat$grp == 1
  c(eta = jb$eta, shift = jb$shift, rep = jb$rep,
    gap_bias = (mean(m[B]) - mean(m[!B])) - (mean(dat$theta[B]) - mean(dat$theta[!B])),
    rmse = sqrt(mean(err^2)), post_prec = mean(1 / apply(fit$theta, 2, var)))
}
res <- mclapply(seq_len(nrow(jobs)), function(k) tryCatch(one(k), error = function(e) {
  message("job ", k, ": ", conditionMessage(e)); rep(NA_real_, 6) }), mc.cores = 4)
R <- cbind(temper = jobs$temper, as.data.frame(do.call(rbind, res)))
agg <- aggregate(cbind(gap_bias, rmse, post_prec) ~ temper + shift + eta, data = R, FUN = mean)
agg <- agg[order(agg$temper, agg$shift, agg$eta), ]
agg$info_frac <- ave(agg$post_prec, agg$temper, agg$shift, FUN = function(x) (x - x[1]) / (x[length(x)] - x[1]))
options(width = 200); print(format(agg, digits = 3), row.names = FALSE)
saveRDS(R, "smi/exp_temper_compare_results.rds")
