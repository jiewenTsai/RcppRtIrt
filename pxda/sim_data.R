# sim_data.R — probit RT-IRT with cross-loadings, generated on the identified scale
sim_pxda <- function(n = 500, p = 15, seed = 1, S12 = 0.3, S22 = 0.4, lam_mean = 0.15, lam_sd = 0.1) {
  set.seed(seed)
  a <- runif(p, .7, 1.5); d <- rnorm(p, 0, .6); xi <- rnorm(p, 4, .2); s <- runif(p, .3, .6)
  lam <- rnorm(p, lam_mean, lam_sd)
  pers <- MASS::mvrnorm(n, c(0, 0), matrix(c(1, S12, S12, S22), 2))
  Y <- matrix(rbinom(n * p, 1, pnorm(outer(pers[, 1], a) - rep(1, n) %o% d)), n)
  logT <- matrix(xi, n, p, byrow = TRUE) - pers[, 2] + outer(pers[, 1], lam) +
    matrix(rnorm(n * p), n) %*% diag(s)
  lb <- mean(lam)
  kappa <- S12 - lb; vstar <- S22 - 2 * lb * S12 + lb^2
  list(Y = Y, logT = logT, a = a, d = d, lam = lam,
       truth = c(kappa = kappa, vstar = vstar, rho_star = kappa / sqrt(vstar)))
}
