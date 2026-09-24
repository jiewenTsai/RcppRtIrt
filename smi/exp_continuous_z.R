# exp_continuous_z.R — continuous speed-related covariate Zc (not in the model) and a placebo W.
# Conditioning model: G in the theta and tau means at every eta. Two tools:
#   1. screen: tau_screen() on the full fit (eta = 1), ordered by Zc and by W (size check);
#   2. risk rule on the target "slope of theta on Zc" (OLS per posterior draw), reference cut+G.
# Shapes (sim_cz): none, linear (0.25 per SD of Zc), threshold (0.5 when Zc > 0.5). 50 datasets each.
suppressMessages(library(parallel))
if (file.exists("../smi/smi_gibbs.R")) setwd("..")
source("pxda/probit_da_gibbs.R"); source("smi/smi_gibbs.R"); source("smi/sim_smi.R")
source("smi/tau_screen.R")
etas <- c(0, .1, .25, .5, .75, 1)
shapes <- c("none", "linear", "threshold")
n_rep <- 50
jobs <- expand.grid(eta = etas, shape = shapes, rep = 1:n_rep, stringsAsFactors = FALSE)
fit_one <- function(k) {
  jb <- jobs[k, ]
  dat <- sim_cz(seed = 3000 + jb$rep, shape = jb$shape)
  X <- cbind(G = dat$G - mean(dat$G))
  fit <- smi_rtirt(dat$Y, dat$logT, eta = jb$eta, n_iter = 2000, n_burn = 500, seed = jb$rep,
                   temper = "marginal", Xth = X, Xtau = if (jb$eta > 0) X else NULL)
  zc <- dat$Zc - mean(dat$Zc)
  slope <- drop(fit$theta %*% zc) / sum(zc^2)                   # OLS slope of theta on Zc per draw
  out <- list(eta = jb$eta, shape = jb$shape, rep = jb$rep, m = mean(slope), V = var(slope),
              q = quantile(slope, c(.025, .975)), truth = sum(dat$theta * zc) / sum(zc^2),
              theta_rmse = sqrt(mean((colMeans(fit$theta) - dat$theta)^2)))
  if (jb$eta == 1) { out$screen_Z <- tau_screen(fit, dat$Zc, X); out$screen_W <- tau_screen(fit, dat$W, X) }
  out
}
t0 <- proc.time()[3]
fits <- mclapply(seq_len(nrow(jobs)), function(k) tryCatch(fit_one(k), error = function(e) {
  message("job ", k, ": ", conditionMessage(e)); NULL }), mc.cores = 4)
cat("elapsed", round(proc.time()[3] - t0), "s; failed", sum(sapply(fits, is.null)), "\n")
saveRDS(fits, "smi/exp_continuous_z_fits.rds")

rows <- list(); scr <- list()
for (sh in shapes) for (r in 1:n_rep) {
  F <- Filter(function(f) !is.null(f) && f$shape == sh && f$rep == r, fits)
  if (length(F) < length(etas)) next
  F <- F[order(sapply(F, `[[`, "eta"))]
  e <- sapply(F, `[[`, "eta"); m <- sapply(F, `[[`, "m"); V <- sapply(F, `[[`, "V")
  D <- m - m[1]; R <- D^2 - V[1] + 2 * V
  full <- F[[length(F)]]
  scr[[length(scr) + 1]] <- data.frame(shape = sh, rep = r,
    rejZ_DM = full$screen_Z["p_DM"] < .05, rejZ_LM2 = full$screen_Z["p_LM2"] < .05,
    rejW_DM = full$screen_W["p_DM"] < .05, rejW_LM2 = full$screen_W["p_LM2"] < .05)
  for (rule in c("cut", "full", "risk")) {
    i <- switch(rule, cut = 1, full = length(F), risk = which.min(R))
    q <- F[[i]]$q; tr <- F[[1]]$truth
    rows[[length(rows) + 1]] <- data.frame(shape = sh, rep = r, rule = rule, eta = e[i],
      err = m[i] - tr, cover = q[1] <= tr & tr <= q[2], theta_rmse = F[[i]]$theta_rmse)
  }
}
out <- do.call(rbind, rows); scr <- do.call(rbind, scr)
saveRDS(list(rules = out, screen = scr), "smi/exp_continuous_z_results.rds")
summ <- aggregate(cbind(eta, bias = err, cover, theta_rmse) ~ shape + rule, out, mean)
summ$rmse <- aggregate(err ~ shape + rule, out, function(x) sqrt(mean(x^2)))$err
summ <- summ[order(match(summ$shape, shapes), match(summ$rule, c("cut", "full", "risk"))), ]
options(width = 200); print(format(summ, digits = 3), row.names = FALSE)
cat("\nscreen rejection rates at .05 (Zc: power, except under 'none'; W: size)\n")
print(aggregate(cbind(rejZ_DM, rejZ_LM2, rejW_DM, rejW_LM2) ~ shape, scr, mean))
