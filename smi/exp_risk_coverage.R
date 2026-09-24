# exp_risk_coverage.R — coverage study of the risk rule (100 datasets per scenario) with the
# conditioning model (G in the theta and tau means), reference cut+G.
# Scenarios (sim_gz, true G gap 0.3 SD): none; Z = 1 slower by 0.25; Z = 1 slower by 0.5.
# Targets: gap by G (in the model) and gap by Z (not in the model).
# For each target: cut, full, risk rule, Hausman rule; bias, RMSE, 95% interval coverage of the
# posterior at the chosen eta (post-selection, as an analyst would report it).
suppressMessages(library(parallel))
if (file.exists("../smi/smi_gibbs.R")) setwd("..")
source("pxda/probit_da_gibbs.R"); source("smi/smi_gibbs.R"); source("smi/sim_smi.R")
etas <- c(0, .1, .25, .5, .75, 1)
scen <- c(none = 0, shiftZ25 = .25, shiftZ50 = .5)
n_rep <- 100
jobs <- expand.grid(eta = etas, scen = names(scen), rep = 1:n_rep, stringsAsFactors = FALSE)
fit_one <- function(k) {
  jb <- jobs[k, ]
  dat <- sim_gz(seed = 2000 + jb$rep, shift_Z = scen[[jb$scen]])
  X <- cbind(G = dat$G - mean(dat$G))
  fit <- smi_rtirt(dat$Y, dat$logT, eta = jb$eta, n_iter = 2000, n_burn = 500, seed = jb$rep,
                   temper = "marginal", Xth = X, Xtau = if (jb$eta > 0) X else NULL)
  gap <- function(b) rowMeans(fit$theta[, b == 1]) - rowMeans(fit$theta[, b == 0])
  tg <- function(b) mean(dat$theta[b == 1]) - mean(dat$theta[b == 0])
  gG <- gap(dat$G); gZ <- gap(dat$Z)
  list(eta = jb$eta, scen = jb$scen, rep = jb$rep,
       m = c(G = mean(gG), Z = mean(gZ)), V = c(G = var(gG), Z = var(gZ)),
       q = rbind(G = quantile(gG, c(.025, .975)), Z = quantile(gZ, c(.025, .975))),
       truth = c(G = tg(dat$G), Z = tg(dat$Z)),
       theta_rmse = sqrt(mean((colMeans(fit$theta) - dat$theta)^2)))
}
t0 <- proc.time()[3]
fits <- mclapply(seq_len(nrow(jobs)), function(k) tryCatch(fit_one(k), error = function(e) {
  message("job ", k, ": ", conditionMessage(e)); NULL }), mc.cores = 4)
cat("elapsed", round(proc.time()[3] - t0), "s; failed", sum(sapply(fits, is.null)), "\n")
saveRDS(fits, "smi/exp_risk_coverage_fits.rds")

rows <- list()
for (sc in names(scen)) for (r in 1:n_rep) {
  F <- Filter(function(f) !is.null(f) && f$scen == sc && f$rep == r, fits)
  if (length(F) < length(etas)) next
  F <- F[order(sapply(F, `[[`, "eta"))]
  e <- sapply(F, `[[`, "eta"); rm <- sapply(F, `[[`, "theta_rmse")
  for (tgt in c("G", "Z")) {
    m <- sapply(F, function(f) f$m[tgt]); V <- sapply(F, function(f) f$V[tgt])
    D <- m - m[1]; H <- ifelse(e == 0, 0, D^2 / pmax(V[1] - V, 1e-12)); R <- D^2 - V[1] + 2 * V
    ok <- which(H < qchisq(.95, 1)); pick_h <- max(ok[ok == seq_along(ok)])
    truth <- F[[1]]$truth[tgt]
    for (rule in c("cut", "full", "risk", "hausman")) {
      i <- switch(rule, cut = 1, full = length(F), risk = which.min(R), hausman = pick_h)
      qi <- F[[i]]$q[tgt, ]
      rows[[length(rows) + 1]] <- data.frame(scen = sc, rep = r, target = tgt, rule = rule, eta = e[i],
        err = m[i] - truth, cover = qi[1] <= truth & truth <= qi[2], width = qi[2] - qi[1],
        theta_rmse = rm[i])
    }
  }
}
out <- do.call(rbind, rows)
saveRDS(out, "smi/exp_risk_coverage_results.rds")
summ <- aggregate(cbind(eta, bias = err, cover, width, theta_rmse) ~ target + scen + rule, out, mean)
summ$rmse <- aggregate(err ~ target + scen + rule, out, function(x) sqrt(mean(x^2)))$err
summ$n <- aggregate(err ~ target + scen + rule, out, length)$err
summ <- summ[order(summ$target, match(summ$scen, names(scen)), match(summ$rule, c("cut", "full", "risk", "hausman"))), ]
options(width = 200); print(format(summ, digits = 3), row.names = FALSE)
cat("\nMC SE of coverage at .95 with 100 datasets: ", round(sqrt(.95 * .05 / n_rep), 3), "\n")
cat("\nrisk-rule eta distribution (Z target)\n")
print(with(subset(out, target == "Z" & rule == "risk"), table(scen, eta)))
