# sim_smi.R — two groups with equal ability distributions. gamma may differ by group and group B
# may be slower overall (tau shifted down); the fitted model assumes one gamma and one speed
# distribution for everyone.
sim_smi <- function(n = 500, p = 15, seed = 1, gam_A = 0.3, gam_B = 0.3, gam_sd = 0.1, v = 0.3,
                    tau_shift_B = 0) {
  set.seed(seed)
  a <- runif(p, .7, 1.5); d <- rnorm(p, 0, .6); xi <- rnorm(p, 4, .2); s <- runif(p, .3, .6)
  gA <- rnorm(p, gam_A, gam_sd); gB <- rnorm(p, gam_B, gam_sd)
  grp <- rbinom(n, 1, .5)                                   # 1 = group B
  th <- rnorm(n); tau <- rnorm(n, 0, sqrt(v)) - tau_shift_B * grp   # group B slower on average
  G <- t(sapply(grp, function(g) if (g == 1) gB else gA))   # n x p person-specific gamma
  Y <- matrix(rbinom(n * p, 1, pnorm(outer(th, a) - rep(1, n) %o% d)), n)
  logT <- matrix(xi, n, p, byrow = TRUE) - tau + G * th + matrix(rnorm(n * p), n) %*% diag(s)
  list(Y = Y, logT = logT, theta = th, a = a, grp = grp, gA = gA, gB = gB)
}
