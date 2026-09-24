# smi_gibbs.R — semi-modular inference (SMI; Carmona & Nicholls 2020) for probit RT-IRT
#
# Identified parameterisation (the shear direction is removed by fixing tau ⟂ theta; the whole
# ability-to-RT relation is carried by gamma_j):
#   module 1 (trusted):  z_ij = a_j theta_i - d_j + e,  Y_ij = 1{z_ij > 0},  theta_i ~ N(0, 1)
#   module 2 (suspect):  log T_ij = xi_j - tau_i + gamma_j theta_i + eps,  eps ~ N(0, sigma_j^2),
#                        tau_i ~ N(0, v)
#   phi = (theta, a, d),  psi = (tau, xi, gamma, sigma, v)
#
# Stage 1 targets p_eta(phi, psi~) ∝ p(Y | phi) p(T | phi, psi~)^eta p(phi) p(psi~):
#   only the RT likelihood is tempered; the prior of tau is not. eta = 0 gives the cut posterior
#   for phi (tau~ then follows its prior and carries no information about theta), eta = 1 the
#   full posterior. Every conditional stays conjugate; the RT precision terms are multiplied by eta.
# Not offered: raising the whole suspect module, [p(T | tau, theta, ...) p(tau | v)]^eta, to the
#   power eta. Integrating tau out of [N(tau; 0, v)]^eta leaves a factor v^((1 - eta)/2) per person,
#   so the auxiliary target is improper in v for eta < 1 (v drifts to infinity). A graded version
#   has to temper the marginal RT likelihood p(T | theta, psi)^eta instead; with tau as a latent
#   variable that needs an extra MH correction |Omega|^(n (1 - eta) / 2) for (sigma^2, v).
# Stage 2 draws psi ~ p(psi | phi, T) (eta = 1) by K_inner Gibbs sweeps per outer iteration,
#   warm-started from the previous psi (nested MCMC; exact as K_inner -> infinity).
#
# Priors: a ~ N(1, 1), d ~ N(0, 2^2), xi ~ N(4, 10^2), gamma ~ N(0, 1), sigma^2 ~ IG(1, 1),
#         v ~ IG(1, 1).

# Requires rtnorm_side() and zscale_move(); run from the repository root.
if (!exists("rtnorm_side")) source("pxda/probit_da_gibbs.R")

smi_rtirt <- function(Y, logT, eta = 1, n_iter = 3000, n_burn = 500, K_inner = 0,
                      zscale_px = TRUE, seed = 1, keep_theta = TRUE) {
  set.seed(seed)
  n <- nrow(Y); p <- ncol(Y)
  a <- rep(1, p); d <- rep(0, p); th <- rep(0, n)
  tau <- rep(0, n); xi <- colMeans(logT); gam <- rep(0, p); s2 <- rep(.25, p); v <- 1
  tau2 <- tau; xi2 <- xi; gam2 <- gam; s22 <- s2; v2 <- v          # stage-2 state
  P0ad <- diag(c(1, 1 / 4)); m0ad <- c(1, 0)
  P0xg <- diag(c(1 / 100, 1)); m0xg <- c(4, 0)
  n_save <- n_iter - n_burn
  TH <- if (keep_theta) matrix(NA_real_, n_save, n) else NULL
  A <- matrix(NA_real_, n_save, p); G1 <- G2 <- matrix(NA_real_, n_save, p)
  V1 <- V2 <- numeric(n_save)

  rt_block <- function(th, tau, xi, gam, s2, v, w_eta, w_prior = 1) {
    # one sweep over (tau, (xi, gamma), sigma^2, v) with the RT likelihood raised to w_eta and the
    # tau prior N(0, v) raised to w_prior
    w <- 1 / s2
    prec_tau <- w_prior / v + w_eta * sum(w)
    rhs <- w_eta * drop((matrix(xi, n, p, byrow = TRUE) + outer(th, gam) - logT) %*% w)
    tau <- rhs / prec_tau + rnorm(n) / sqrt(prec_tau)
    W <- cbind(1, th); WtW <- crossprod(W); Wty <- crossprod(W, logT + tau)
    for (j in 1:p) {
      Vj <- solve(P0xg + w_eta * WtW / s2[j])
      mj <- Vj %*% (P0xg %*% m0xg + w_eta * Wty[, j] / s2[j])
      xg <- mj + t(chol(Vj)) %*% rnorm(2)
      xi[j] <- xg[1]; gam[j] <- xg[2]
    }
    E <- logT - (matrix(xi, n, p, byrow = TRUE) - tau + outer(th, gam))
    s2 <- 1 / rgamma(p, 1 + w_eta * n / 2, 1 + w_eta * colSums(E^2) / 2)
    v <- 1 / rgamma(1, 1 + w_prior * n / 2, 1 + w_prior * sum(tau^2) / 2)
    list(tau = tau, xi = xi, gam = gam, s2 = s2, v = v)
  }

  for (it in 1:n_iter) {
    ## ---- stage 1 --------------------------------------------------------------
    Z <- matrix(rtnorm_side(outer(th, a) - matrix(d, n, p, byrow = TRUE), Y), n, p)
    # theta | z, items, tempered RT
    w <- 1 / s2
    prec_th <- 1 + sum(a^2) + eta * sum(gam^2 * w)
    Rth <- logT - matrix(xi, n, p, byrow = TRUE) + tau                 # = gamma theta + eps
    b_th <- drop((Z + matrix(d, n, p, byrow = TRUE)) %*% a) + eta * drop(Rth %*% (gam * w))
    th <- b_th / prec_th + rnorm(n) / sqrt(prec_th)
    # items (a, d) with optional z-scale PX-DA
    X <- cbind(th, -1); Vad <- solve(P0ad + crossprod(X))
    if (zscale_px) for (j in 1:p) Z[, j] <- Z[, j] * zscale_move(Z[, j], X, P0ad, m0ad, Vad)
    Mad <- Vad %*% (drop(P0ad %*% m0ad) + crossprod(X, Z))
    AD <- Mad + t(chol(Vad)) %*% matrix(rnorm(2 * p), 2)
    a <- AD[1, ]; d <- AD[2, ]
    # auxiliary RT parameters psi~ under the tempered likelihood
    rb <- rt_block(th, tau, xi, gam, s2, v, eta)
    tau <- rb$tau; xi <- rb$xi; gam <- rb$gam; s2 <- rb$s2; v <- rb$v

    ## ---- stage 2: psi | theta, T (full likelihood), nested ---------------------
    if (K_inner > 0) {
      for (k in 1:K_inner) {
        rb2 <- rt_block(th, tau2, xi2, gam2, s22, v2, 1)
        tau2 <- rb2$tau; xi2 <- rb2$xi; gam2 <- rb2$gam; s22 <- rb2$s2; v2 <- rb2$v
      }
    }
    if (it > n_burn) {
      k <- it - n_burn
      if (keep_theta) TH[k, ] <- th
      A[k, ] <- a; G1[k, ] <- gam; V1[k] <- v
      if (K_inner > 0) { G2[k, ] <- gam2; V2[k] <- v2 }
    }
  }
  list(theta = TH, a = A, gamma_stage1 = G1, v_stage1 = V1,
       gamma_stage2 = if (K_inner > 0) G2 else NULL, v_stage2 = if (K_inner > 0) V2 else NULL,
       eta = eta, K_inner = K_inner)
}
