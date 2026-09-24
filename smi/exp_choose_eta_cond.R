# exp_choose_eta_cond.R — eta rules with a conditioning model (theta and tau means include G), so
# the reference is cut+G, which is unbiased for a real group gap (exp_hier_vs_cut.R).
# Data (sim_gz): true ability gap 0.3 SD by G; Z correlated with G, no ability effect.
#   none   : no speed shift
#   shiftG : G = 1 slower by 0.5 (the conditioning model can absorb it)
#   shiftZ : Z = 1 slower by 0.5 (not in the model)
# Rules are applied separately for each target: the gap by G and the gap by Z.
#   D = gap(eta) - gap(0),  V = posterior variance of the gap
#   hausman: largest eta before the first H = D^2 / (V(0) - V(eta)) >= 3.84;  risk: argmin D^2 - V(0) + 2 V(eta)
suppressMessages(library(parallel))
if (file.exists("../smi/smi_gibbs.R")) setwd("..")
source("pxda/probit_da_gibbs.R"); source("smi/smi_gibbs.R"); source("smi/sim_smi.R")
etas <- c(0, .1, .25, .4, .55, .75, 1)
scen <- list(none = c(0, 0), shiftG = c(.5, 0), shiftZ = c(0, .5))     # (shift_G, shift_Z)
jobs <- expand.grid(eta = etas, scen = names(scen), rep = 1:20, stringsAsFactors = FALSE)
fit_one <- function(k) {
  jb <- jobs[k, ]; sh <- scen[[jb$scen]]
  dat <- sim_gz(seed = 900 + jb$rep, shift_G = sh[1], shift_Z = sh[2])
  X <- cbind(G = dat$G - mean(dat$G))
  fit <- smi_rtirt(dat$Y, dat$logT, eta = jb$eta, n_iter = 2500, n_burn = 500, seed = jb$rep,
                   temper = "marginal", Xth = X, Xtau = if (jb$eta > 0) X else NULL)
  gap <- function(b) rowMeans(fit$theta[, b == 1]) - rowMeans(fit$theta[, b == 0])
  tg <- function(b) mean(dat$theta[b == 1]) - mean(dat$theta[b == 0])
  gG <- gap(dat$G); gZ <- gap(dat$Z)
  list(eta = jb$eta, scen = jb$scen, rep = jb$rep,
       m = c(G = mean(gG), Z = mean(gZ)), V = c(G = var(gG), Z = var(gZ)),
       truth = c(G = tg(dat$G), Z = tg(dat$Z)),
       theta_rmse = sqrt(mean((colMeans(fit$theta) - dat$theta)^2)))
}
t0 <- proc.time()[3]
fits <- mclapply(seq_len(nrow(jobs)), function(k) tryCatch(fit_one(k), error = function(e) {
  message("job ", k, ": ", conditionMessage(e)); NULL }), mc.cores = 4)
cat("elapsed", round(proc.time()[3] - t0), "s; failed", sum(sapply(fits, is.null)), "\n")
saveRDS(fits, "smi/exp_choose_eta_cond_fits.rds")

rows <- list()
for (sc in names(scen)) for (r in unique(jobs$rep)) {
  F <- Filter(function(f) !is.null(f) && f$scen == sc && f$rep == r, fits)
  F <- F[order(sapply(F, `[[`, "eta"))]
  e <- sapply(F, `[[`, "eta"); rm <- sapply(F, `[[`, "theta_rmse")
  for (tgt in c("G", "Z")) {
    m <- sapply(F, function(f) f$m[tgt]); V <- sapply(F, function(f) f$V[tgt])
    D <- m - m[1]; H <- ifelse(e == 0, 0, D^2 / pmax(V[1] - V, 1e-12)); R <- D^2 - V[1] + 2 * V
    ok <- which(H < qchisq(.95, 1)); pick_h <- max(ok[ok == seq_along(ok)])
    truth <- F[[1]]$truth[tgt]
    for (rule in c("cut", "full", "risk", "hausman")) {
      i <- switch(rule, cut = 1, full = length(F), risk = which.min(R), hausman = pick_h)
      rows[[length(rows) + 1]] <- data.frame(scen = sc, rep = r, target = tgt, rule = rule, eta = e[i],
        err = m[i] - truth, cover = abs(m[i] - truth) < 1.96 * sqrt(V[i]), theta_rmse = rm[i])
    }
  }
}
out <- do.call(rbind, rows)
summ <- aggregate(cbind(eta, bias = err, cover, theta_rmse) ~ target + scen + rule, out, mean)
summ$rmse <- aggregate(err ~ target + scen + rule, out, function(x) sqrt(mean(x^2)))$err
summ <- summ[order(summ$target, match(summ$scen, names(scen)), match(summ$rule, c("cut", "full", "risk", "hausman"))), ]
options(width = 200); print(format(summ, digits = 3), row.names = FALSE)
saveRDS(out, "smi/exp_choose_eta_cond_results.rds")
