source("proto/sim_drift.R")
library(parallel)
run_fast <- function(dat) {
  fit <- gibbs_rtirt_null(dat$Y, dat$logT, n_iter = 1000, n_burnin = 500, verbose = FALSE)
  gh <- gh_nodes(9)
  ml <- bhhh_mml(psi_from_gibbs(fit), dat$Y, dat$logT, gh, max_iter = 5)
  S <- casewise_scores(ml$psi, dat$Y, dat$logT, gh)
  score_instability(S, 4 * dat$p + 1, dat$z)
}
conds <- list(
  H0         = list(f = function(z) rep(.35, length(z)),          R = 100),
  H1_linear  = list(f = function(z) .55 - .40 * z,                R = 50),
  H1_step70  = list(f = function(z) ifelse(z < .7, .45, .15),     R = 50))
out <- list()
for (nm in names(conds)) {
  cc <- conds[[nm]]
  res <- mclapply(seq_len(cc$R), function(r) {
    d <- sim_drift(1000, 15, cc$f, seed = 1000 * match(nm, names(conds)) + r)
    tryCatch(run_fast(d), error = function(e) rep(NA_real_, 4))
  }, mc.cores = 4)
  m <- do.call(rbind, res)
  out[[nm]] <- c(R = sum(complete.cases(m)),
                 rej_DM = mean(m[, 2] < .05, na.rm = TRUE),
                 rej_LM_median = mean(m[, 4] < .05, na.rm = TRUE))
  cat(nm, "done\n")
  saveRDS(out, "proto/sim_out.rds")
}
print(round(do.call(rbind, out), 3))
