# bench_collapse.R — expanded probit DA: standard vs partially collapsed vs z-scale PX-DA
suppressMessages({library(coda); library(parallel)})
if (file.exists("../pxda/probit_da_gibbs.R")) setwd("..")
source("pxda/probit_da_gibbs.R"); source("pxda/sim_data.R")
one <- function(seed, variant) {
  dat <- sim_pxda(500, 15, seed = seed)
  t0 <- proc.time()[3]
  S <- probit_da_gibbs(dat$Y, dat$logT, 6000, 1000, "expanded", collapse_items = variant == 1,
                       zscale_px = variant == 2, seed = seed)
  secs <- unname(proc.time()[3] - t0)
  idd <- identified_draws(S)
  c(seed = seed, variant = variant, secs = secs,
    ESS_kappa = unname(effectiveSize(idd$summary[, "kappa"])),
    ESS_rho = unname(effectiveSize(idd$summary[, "rho_star"])),
    ESS_a_med = median(effectiveSize(mcmc(idd$a_id))),
    ESS_a_min = min(effectiveSize(mcmc(idd$a_id))),
    ESS_lamdev_med = median(effectiveSize(mcmc(idd$lamdev_id))),
    corr_a = cor(colMeans(idd$a_id), dat$a))
}
jobs <- expand.grid(seed = 1:4, variant = 0:2)   # 0 standard, 1 collapsed theta, 2 z-scale PX
res <- do.call(rbind, mclapply(seq_len(nrow(jobs)), function(k) one(jobs$seed[k], jobs$variant[k]), mc.cores = 4))
res <- as.data.frame(res)
print(round(res, 3))
agg <- aggregate(. ~ variant, data = res[, -1], FUN = mean)
agg$ESSps_a_med <- agg$ESS_a_med / agg$secs; agg$ESSps_rho <- agg$ESS_rho / agg$secs
cat("\nmean over datasets\n"); print(round(agg, 2))
