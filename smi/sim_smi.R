# sim_smi.R — two groups with equal ability distributions. gamma may differ by group and group B
# may be slower overall (tau shifted down); the fitted model assumes one gamma and one speed
# distribution for everyone. gamma_j is drawn separately for each group (gam_sd around the group
# mean), so even with gam_A == gam_B the two groups differ item by item (mild heterogeneity).
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

# sim_gz — reporting group G with a true ability gap gap_G; Z correlated with G (P(Z=1|G) = .3/.7)
# with no ability effect. Speed shifts: Z = 1 slower by shift_Z, G = 1 slower by shift_G (in tau).
sim_gz <- function(n = 500, p = 10, seed = 1, gap_G = .3, shift_Z = .5, gam = .5, v = .15,
                   shift_G = 0) {
  set.seed(seed)
  a <- runif(p, .7, 1.5); d <- rnorm(p, 0, .6); xi <- rnorm(p, 4, .2); s <- runif(p, .3, .6)
  g <- rnorm(p, gam, .1)
  G <- rbinom(n, 1, .5); Z <- rbinom(n, 1, ifelse(G == 1, .7, .3))
  th <- rnorm(n) + gap_G * (G - .5)
  tau <- rnorm(n, 0, sqrt(v)) - shift_Z * Z - shift_G * G
  Y <- matrix(rbinom(n * p, 1, pnorm(outer(th, a) - rep(1, n) %o% d)), n)
  logT <- matrix(xi, n, p, byrow = TRUE) - tau + outer(th, g) + matrix(rnorm(n * p), n) %*% diag(s)
  list(Y = Y, logT = logT, theta = th, G = G, Z = Z)
}

# sim_cz — like sim_gz but the speed-related variable Zc is continuous, Zc = 0.8 (G - 1/2) + N(0, 1),
# with no ability effect. Speed shift in tau: "none", "linear" (Zc slower by 0.25 per SD), or
# "threshold" (slower by 0.5 when Zc > 0.5). W ~ N(0, 1) is an unrelated placebo covariate.
sim_cz <- function(n = 500, p = 10, seed = 1, gap_G = .3, shape = c("none", "linear", "threshold"),
                   gam = .5, v = .15) {
  shape <- match.arg(shape)
  set.seed(seed)
  a <- runif(p, .7, 1.5); d <- rnorm(p, 0, .6); xi <- rnorm(p, 4, .2); s <- runif(p, .3, .6)
  g <- rnorm(p, gam, .1)
  G <- rbinom(n, 1, .5); Zc <- .8 * (G - .5) + rnorm(n); W <- rnorm(n)
  th <- rnorm(n) + gap_G * (G - .5)
  shift <- switch(shape, none = 0, linear = .25 * Zc, threshold = .5 * (Zc > .5))
  tau <- rnorm(n, 0, sqrt(v)) - shift
  Y <- matrix(rbinom(n * p, 1, pnorm(outer(th, a) - rep(1, n) %o% d)), n)
  logT <- matrix(xi, n, p, byrow = TRUE) - tau + outer(th, g) + matrix(rnorm(n * p), n) %*% diag(s)
  list(Y = Y, logT = logT, theta = th, G = G, Zc = Zc, W = W)
}
