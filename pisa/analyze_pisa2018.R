# analyze_pisa2018.R — empirical illustration of SMI RT-IRT on a prepared PISA 2018 block.
# Conditioning model: gender (female) in the theta and tau means. Targets:
#   gap_female  theta gap female - male        (in the model)
#   gap_ml      theta gap other-language-at-home - test-language   (not in the model)
#   slope_escs  OLS slope of theta on ESCS      (not in the model)
# For each target: estimates across eta, the risk rule's eta, and the cut / full contrast.
# Diagnostics on the full fit: tau-score screen by ml and ESCS, item-level RT residual gaps by ml
# (differential response time), leakage coefficient lambda.
#
# Usage: Rscript pisa/analyze_pisa2018.R data/pisa2018/prepared_USA_CM.rds [n_iter]
suppressMessages(library(parallel))
source("pxda/probit_da_gibbs.R"); source("smi/smi_gibbs.R"); source("smi/tau_screen.R")
source("smi/leakage.R")
args <- commandArgs(trailingOnly = TRUE)
f <- if (length(args) >= 1) args[1] else "data/pisa2018/prepared_USA_CM.rds"
n_iter <- if (length(args) >= 2) as.integer(args[2]) else 4000
d <- readRDS(f)
n <- nrow(d$Y); p <- ncol(d$Y)
cat(sprintf("%s %s: n = %d, p = %d\n", d$country, d$domain, n, p))
X <- cbind(female = d$female - mean(d$female))
etas <- c(0, .1, .25, .5, .75, 1)

fits <- mclapply(etas, function(e) smi_rtirt(d$Y, d$logT, eta = e, n_iter = n_iter,
  n_burn = n_iter %/% 4, seed = 1, temper = "marginal", Xth = X,
  Xtau = if (e > 0) X else NULL, keep_tau = e == 1), mc.cores = min(4, length(etas)))
names(fits) <- etas

gap <- function(M, b) { ok <- !is.na(b); rowMeans(M[, ok & b == 1, drop = FALSE]) - rowMeans(M[, ok & b == 0, drop = FALSE]) }
slope <- function(M, z) { ok <- !is.na(z); zc <- z[ok] - mean(z[ok]); drop(M[, ok] %*% zc) / sum(zc^2) }
targets <- list(gap_female = function(M) gap(M, d$female), gap_ml = function(M) gap(M, d$ml),
                slope_escs = function(M) slope(M, d$escs))
est <- do.call(rbind, lapply(names(targets), function(tn) do.call(rbind, lapply(seq_along(etas), function(k) {
  x <- targets[[tn]](fits[[k]]$theta)
  data.frame(target = tn, eta = etas[k], mean = mean(x), sd = sd(x),
             lo = quantile(x, .025), hi = quantile(x, .975),
             theta_psd = mean(apply(fits[[k]]$theta, 2, sd)))
}))))
rownames(est) <- NULL
est <- do.call(rbind, lapply(split(est, est$target), function(e) {
  e <- e[order(e$eta), ]; D <- e$mean - e$mean[1]
  e$D <- D; e$risk <- D^2 - e$sd[1]^2 + 2 * e$sd^2; e$chosen <- seq_len(nrow(e)) == which.min(e$risk); e
}))

full <- fits[["1"]]
screen <- rbind(ml = tau_screen(full, d$ml, X), escs = tau_screen(full, d$escs, X))
# item-level RT residuals at the full fit: logT - xi + tau - gamma theta, then gap by ml per item
th_hat <- colMeans(full$theta); g_hat <- colMeans(full$gamma_stage1)
R <- d$logT - matrix(full$xi_mean, n, p, byrow = TRUE) + full$tau_mean - outer(th_hat, g_hat)
drt <- t(sapply(seq_len(p), function(j) {
  cf <- summary(lm(R[, j] ~ d$ml + d$female))$coefficients
  c(gamma = g_hat[j], gap_ml = cf[2, 1], t = cf[2, 3])
}))
rownames(drt) <- d$items
L <- leakage(full)

res <- list(data = d$log, estimates = est, screen = screen, drt = drt, leakage = L[c("lambda", "c")],
            rt = list(gamma = colMeans(full$gamma_stage1), v = mean(full$v_stage1),
                      beta_tau = colMeans(full$beta_tau), beta_theta = colMeans(full$beta_theta)))
out <- sub("prepared_", "results_", f)
saveRDS(res, out)
options(width = 200)
cat("\nestimates by eta (risk-rule choice marked)\n")
print(format(est[, c("target", "eta", "mean", "sd", "lo", "hi", "D", "risk", "chosen", "theta_psd")], digits = 3), row.names = FALSE)
cat("\ntau-score screen (full fit)\n"); print(round(screen, 4))
cat("\nitem RT residual gap by ml (differential response time), gamma\n"); print(round(drt, 3))
cat("\nleakage lambda =", round(L$lambda, 3), ", c =", round(L$c, 3), "\n")
cat("RT module: v =", round(res$rt$v, 3), ", beta_tau(female) =", round(res$rt$beta_tau, 3),
    ", beta_theta(female) =", round(res$rt$beta_theta, 3), "\n")
cat("saved", out, "\n")
