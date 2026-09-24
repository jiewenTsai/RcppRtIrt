# exp_hier_vs_cut.R — what does the cut add over a hierarchical (conditioning) model?
# Data: G = reporting group with a true ability gap (0.3 SD); Z = a second variable (e.g. language
# background, accommodation) correlated with G (P(Z=1|G) = .3/.7) that carries a speed shift
# (Z = 1 slower by 0.5 in tau) and no ability effect of its own. Targets: theta gap by G and by Z.
# Methods (eta = 1 "full", eta = 0 "cut"), covariates entering the theta mean (Xth) and tau mean (Xtau):
#   full        : none                 cut        : none
#   full+G      : G in theta and tau   cut+G      : G in theta
#   full+GZ     : G, Z in both         cut+GZ     : G, Z in theta
suppressMessages(library(parallel))
if (file.exists("../smi/smi_gibbs.R")) setwd("..")
source("pxda/probit_da_gibbs.R"); source("smi/smi_gibbs.R"); source("smi/sim_smi.R")

methods <- data.frame(name = c("full", "full+G", "full+GZ", "cut", "cut+G", "cut+GZ"),
                      eta = c(1, 1, 1, 0, 0, 0), cov = c("", "G", "GZ", "", "G", "GZ"))
jobs <- expand.grid(m = seq_len(nrow(methods)), shift = c(0, .5), rep = 1:20)
fit_one <- function(k) {
  jb <- jobs[k, ]; me <- methods[jb$m, ]
  dat <- sim_gz(seed = 700 + jb$rep, shift_Z = jb$shift)
  X <- switch(me$cov, G = cbind(G = dat$G - mean(dat$G)),
              GZ = cbind(G = dat$G - mean(dat$G), Z = dat$Z - mean(dat$Z)), NULL)
  fit <- smi_rtirt(dat$Y, dat$logT, eta = me$eta, n_iter = 2500, n_burn = 500, seed = jb$rep,
                   temper = "likelihood", Xth = X, Xtau = if (me$eta > 0) X else NULL)
  gap <- function(b) rowMeans(fit$theta[, b == 1]) - rowMeans(fit$theta[, b == 0])
  tg <- function(b) mean(dat$theta[b == 1]) - mean(dat$theta[b == 0])
  gG <- gap(dat$G); gZ <- gap(dat$Z)
  data.frame(method = me$name, shift = jb$shift, rep = jb$rep,
             errG = mean(gG) - tg(dat$G), errZ = mean(gZ) - tg(dat$Z),
             sdG = sd(gG), sdZ = sd(gZ),
             covG = abs(mean(gG) - tg(dat$G)) < 1.96 * sd(gG),
             covZ = abs(mean(gZ) - tg(dat$Z)) < 1.96 * sd(gZ),
             theta_rmse = sqrt(mean((colMeans(fit$theta) - dat$theta)^2)))
}
t0 <- proc.time()[3]
res <- do.call(rbind, mclapply(seq_len(nrow(jobs)), function(k) tryCatch(fit_one(k), error = function(e) {
  message("job ", k, ": ", conditionMessage(e)); NULL }), mc.cores = 4))
cat("elapsed", round(proc.time()[3] - t0), "s\n")
saveRDS(res, "smi/exp_hier_vs_cut_results.rds")
summ <- aggregate(cbind(biasG = errG, biasZ = errZ, covG, covZ, theta_rmse) ~ method + shift, res, mean)
summ$rmseG <- aggregate(errG ~ method + shift, res, function(x) sqrt(mean(x^2)))$errG
summ$rmseZ <- aggregate(errZ ~ method + shift, res, function(x) sqrt(mean(x^2)))$errZ
summ <- summ[order(summ$shift, match(summ$method, methods$name)), ]
options(width = 200); print(format(summ, digits = 3), row.names = FALSE)
