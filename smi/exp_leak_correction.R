# exp_leak_correction.R — post-hoc leakage correction for secondary analyses (leakage.R).
# The producer fits the full joint model once (eta = 1, G in theta and tau means) and releases draws
# of theta and of the speed residual u = tau - m_tau, plus the leakage ratio c. A secondary analyst
# estimates a target on a variable the model did not condition on:
#   binary Z (sim_gz, shift_Z in {0, .25, .5}): theta gap by Z
#   continuous Zc (sim_cz, none / linear / threshold): OLS slope of theta on Zc
# raw = from theta draws; corrected = target(theta) - c * target(u), per draw; cut = cut+G fit.
suppressMessages(library(parallel))
if (file.exists("../smi/smi_gibbs.R")) setwd("..")
source("pxda/probit_da_gibbs.R"); source("smi/smi_gibbs.R"); source("smi/sim_smi.R"); source("smi/leakage.R")
scen <- data.frame(name = c("Z0", "Z25", "Z50", "Cnone", "Clinear", "Cthreshold"),
                   kind = rep(c("gz", "cz"), each = 3),
                   shift = c(0, .25, .5, NA, NA, NA),
                   shape = c(NA, NA, NA, "none", "linear", "threshold"), stringsAsFactors = FALSE)
n_rep <- 50
jobs <- expand.grid(s = seq_len(nrow(scen)), rep = 1:n_rep)
target_fun <- function(z, binary) {
  if (binary) function(M) rowMeans(M[, z == 1, drop = FALSE]) - rowMeans(M[, z == 0, drop = FALSE])
  else { zc <- z - mean(z); function(M) drop(M %*% zc) / sum(zc^2) }
}
fit_one <- function(k) {
  jb <- jobs[k, ]; sc <- scen[jb$s, ]
  if (nzchar(Sys.getenv("SMI_PROGRESS"))) cat(k, "\n", file = Sys.getenv("SMI_PROGRESS"), append = TRUE)
  dat <- if (sc$kind == "gz") sim_gz(seed = 6000 + jb$rep, shift_Z = sc$shift) else
    sim_cz(seed = 6000 + jb$rep, shape = sc$shape)
  z <- if (sc$kind == "gz") dat$Z else dat$Zc
  tf <- target_fun(z, sc$kind == "gz")
  truth <- tf(matrix(dat$theta, 1))
  X <- cbind(G = dat$G - mean(dat$G))
  full <- smi_rtirt(dat$Y, dat$logT, eta = 1, n_iter = 2000, n_burn = 500, seed = jb$rep,
                    temper = "marginal", Xth = X, Xtau = X, keep_tau = TRUE)
  cut <- smi_rtirt(dat$Y, dat$logT, eta = 0, n_iter = 2000, n_burn = 500, seed = jb$rep, Xth = X)
  L <- leakage(full)
  draws <- list(raw = tf(full$theta), corrected = tf(full$theta) - L$c * tf(full$tau_resid),
                cut = tf(cut$theta))
  do.call(rbind, lapply(names(draws), function(m) {
    x <- draws[[m]]; q <- quantile(x, c(.025, .975))
    data.frame(scen = sc$name, rep = jb$rep, method = m, err = mean(x) - truth, sd = sd(x),
               cover = q[1] <= truth & truth <= q[2], c = L$c, lambda = L$lambda)
  }))
}
t0 <- proc.time()[3]
res <- do.call(rbind, mclapply(seq_len(nrow(jobs)), function(k) tryCatch(fit_one(k), error = function(e) {
  message("job ", k, ": ", conditionMessage(e)); NULL }), mc.cores = 4))
cat("elapsed", round(proc.time()[3] - t0), "s\n")
saveRDS(res, "smi/exp_leak_correction_results.rds")
summ <- aggregate(cbind(bias = err, sd, cover) ~ method + scen, res, mean)
summ$rmse <- aggregate(err ~ method + scen, res, function(x) sqrt(mean(x^2)))$err
summ$n <- aggregate(err ~ method + scen, res, length)$err
summ <- summ[order(match(summ$scen, scen$name), match(summ$method, c("raw", "corrected", "cut"))), ]
options(width = 200); print(format(summ, digits = 3), row.names = FALSE)
cat("\nleakage ratio c and lambda (mean by scenario)\n")
print(aggregate(cbind(c, lambda) ~ scen, subset(res, method == "raw"), mean))
