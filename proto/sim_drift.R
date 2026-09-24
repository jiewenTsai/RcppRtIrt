suppressMessages(library(Rcpp))
if (file.exists("../aghq_score.R")) setwd("..")
sourceCpp("RtIrtGibbs.cpp"); source("aghq_score.R")

# DGP: Sigma_11 = 1, Sigma_22 = 0.4; Sigma_12 depends on covariate z ~ U(0,1)
sim_drift <- function(n, p, s12_fun, seed) {
  set.seed(seed)
  a <- runif(p, .7, 1.5); b <- rnorm(p, 0, .6); xi <- rnorm(p, 4, .2); sig <- runif(p, .3, .6)
  z <- runif(n); s12 <- s12_fun(z); S22 <- 0.4
  th <- rnorm(n); ta <- s12 * th + rnorm(n) * sqrt(S22 - s12^2)
  Y <- matrix(rbinom(n * p, 1, plogis(outer(th, a) - rep(1, n) %o% (a * b))), n)
  logT <- matrix(xi, n, p, byrow = TRUE) - ta + matrix(rnorm(n * p), n) %*% diag(sig)
  list(Y = Y, logT = logT, z = z, p = p)
}

run_one <- function(dat, K = 9, verbose = FALSE) {
  t0 <- proc.time()[3]
  fit <- gibbs_rtirt_null(dat$Y, dat$logT, n_iter = 2000, n_burnin = 1000, verbose = FALSE)
  t1 <- proc.time()[3]
  gh <- gh_nodes(K)
  psi0 <- psi_from_gibbs(fit)
  ml <- bhhh_mml(psi0, dat$Y, dat$logT, gh, verbose = verbose)
  S <- casewise_scores(ml$psi, dat$Y, dat$logT, gh)
  t2 <- proc.time()[3]
  k12 <- 4 * dat$p + 1
  c(score_instability(S, k12, dat$z),
    s12_gibbs = psi0[k12], s12_mml = ml$psi[k12], bhhh_iter = ml$iter,
    max_abs_score_mean = max(abs(colMeans(S))),
    sec_gibbs = unname(t1 - t0), sec_aghq = unname(t2 - t1))
}
