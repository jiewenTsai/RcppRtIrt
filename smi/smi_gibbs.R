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
# temper = "marginal" raises the marginal RT likelihood L_i = N(T_i; xi + gamma theta_i, Omega),
#   Omega = diag(sigma^2) + v 11' (tau integrated out), to the power eta. It is sampled with the
#   augmentation tau~_i ~ N(0, v / eta), T_i | tau~ ~ N(xi - tau~ + gamma theta_i, diag(sigma^2) / eta),
#   whose tau~-integral is N(T_i; mu_i, Omega / eta) = const * L_i^eta * |Omega|^(-(1 - eta)/2)
#   (check_marginal_identity.R), so the augmented target carries the correction |Omega|^(n (1-eta)/2),
#   |Omega| = prod sigma_j^2 (1 + v sum_j 1/sigma_j^2). tau~, theta, (xi, gamma) stay conjugate;
#   sigma_j^2 and v are drawn by univariate slice sampling on the log scale. eta = 0 is the cut.
# Not offered: raising the whole module with the tau prior as a density, [p(T | tau, ...) p(tau | v)]^eta.
#   Integrating tau out of [N(tau; 0, v)]^eta leaves v^((1 - eta)/2) per person, so the auxiliary
#   target is improper in v for eta < 1.
# Stage 2 draws psi ~ p(psi | phi, T) (eta = 1) by K_inner Gibbs sweeps per outer iteration,
#   warm-started from the previous psi (nested MCMC; exact as K_inner -> infinity).
#
# Priors: a ~ N(1, 1), d ~ N(0, 2^2), xi ~ N(4, 10^2), gamma ~ N(0, 1), sigma^2 ~ IG(1, 1),
#         v ~ IG(1, 1).
# Optional person covariates (hierarchical / conditioning model; no intercept column, centre them):
#   Xth:  theta_i ~ N(Xth_i beta_th, 1);   Xtau: tau_i ~ N(Xtau_i beta_tau, v);   beta ~ N(0, 10 I).
#   Xtau needs K_inner = 0. Under temper = "marginal" the augmentation becomes
#   tau~_i ~ N(Xtau_i beta_tau, v / eta), so the marginal is L_i^eta with mean xi - Xtau_i beta_tau
#   + gamma theta_i and the |Omega| correction is unchanged; beta_tau | tau~ uses variance v / eta.

# Requires rtnorm_side() and zscale_move(); run from the repository root.
if (!exists("rtnorm_side")) source("pxda/probit_da_gibbs.R")

smi_rtirt <- function(Y, logT, eta = 1, n_iter = 3000, n_burn = 500, K_inner = 0,
                      zscale_px = TRUE, seed = 1, keep_theta = TRUE,
                      temper = c("likelihood", "marginal"), use_cpp = FALSE,
                      Xth = NULL, Xtau = NULL) {
  temper <- match.arg(temper)
  if (!is.null(Xtau) && K_inner > 0) stop("Xtau needs K_inner = 0")
  if (use_cpp && !exists("rt_sweeps_cpp")) Rcpp::sourceCpp("smi/rt_sweeps.cpp")
  set.seed(seed)
  n <- nrow(Y); p <- ncol(Y)
  a <- rep(1, p); d <- rep(0, p); th <- rep(0, n)
  tau <- rep(0, n); xi <- colMeans(logT); gam <- rep(0, p); s2 <- rep(.25, p); v <- 1
  tau2 <- tau; xi2 <- xi; gam2 <- gam; s22 <- s2; v2 <- v          # stage-2 state
  m_th <- m_tau <- rep(0, n)                                         # covariate means
  if (!is.null(Xth)) { Xth <- as.matrix(Xth); bet_th <- rep(0, ncol(Xth)) }
  if (!is.null(Xtau)) { Xtau <- as.matrix(Xtau); b_tau <- rep(0, ncol(Xtau)) }
  draw_beta <- function(X, y, s2y) {                                 # y ~ N(X b, s2y), b ~ N(0, 10 I)
    Vb <- solve(crossprod(X) / s2y + diag(ncol(X)) / 10)
    drop(Vb %*% crossprod(X, y) / s2y + t(chol(Vb)) %*% rnorm(ncol(X)))
  }
  B_th <- B_tau <- NULL
  S_tau <- rep(0, n)                          # running sum of the casewise tau-location score
  P0ad <- diag(c(1, 1 / 4)); m0ad <- c(1, 0)
  P0xg <- diag(c(1 / 100, 1)); m0xg <- c(4, 0)
  n_save <- n_iter - n_burn
  TH <- if (keep_theta) matrix(NA_real_, n_save, n) else NULL
  A <- matrix(NA_real_, n_save, p); G1 <- G2 <- matrix(NA_real_, n_save, p)
  V1 <- V2 <- numeric(n_save)
  if (!is.null(Xth)) B_th <- matrix(NA_real_, n_save, ncol(Xth))
  if (!is.null(Xtau)) B_tau <- matrix(NA_real_, n_save, ncol(Xtau))

  rt_block <- function(th, tau, xi, gam, s2, v, w_eta, w_prior = 1, m_tau = 0) {
    # one sweep over (tau, (xi, gamma), sigma^2, v) with the RT likelihood raised to w_eta and the
    # tau prior N(m_tau, v) raised to w_prior
    w <- 1 / s2
    prec_tau <- w_prior / v + w_eta * sum(w)
    rhs <- w_eta * drop((matrix(xi, n, p, byrow = TRUE) + outer(th, gam) - logT) %*% w) +
      w_prior * m_tau / v
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
    v <- 1 / rgamma(1, 1 + w_prior * n / 2, 1 + w_prior * sum((tau - m_tau)^2) / 2)
    list(tau = tau, xi = xi, gam = gam, s2 = s2, v = v)
  }

  rt_block_marginal <- function(th, tau, xi, gam, s2, v, eta, m_tau = 0) {
    # one sweep under the marginally tempered RT module (see header)
    w <- 1 / s2
    prec_tau <- eta / v + eta * sum(w)
    rhs <- eta * drop((matrix(xi, n, p, byrow = TRUE) + outer(th, gam) - logT) %*% w) + eta * m_tau / v
    tau <- rhs / prec_tau + rnorm(n) / sqrt(prec_tau)
    W <- cbind(1, th); WtW <- crossprod(W); Wty <- crossprod(W, logT + tau)
    for (j in 1:p) {
      Vj <- solve(P0xg + eta * WtW / s2[j])
      mj <- Vj %*% (P0xg %*% m0xg + eta * Wty[, j] / s2[j])
      xg <- mj + t(chol(Vj)) %*% rnorm(2)
      xi[j] <- xg[1]; gam[j] <- xg[2]
    }
    E <- logT - (matrix(xi, n, p, byrow = TRUE) - tau + outer(th, gam))
    SSR <- colSums(E^2)
    kpow <- n * (1 - eta) / 2
    for (j in 1:p) {        # sigma_j^2: IG(1 + eta n/2, 1 + eta SSR_j/2) x (1 + v sum_k 1/sigma_k^2)^kpow
      S_other <- sum(1 / s2[-j])
      lf <- function(u) -(1 + eta * n / 2) * u - (1 + eta * SSR[j] / 2) * exp(-u) +
        kpow * log1p(v * (S_other + exp(-u)))
      s2[j] <- exp(slice1(log(s2[j]), lf, width = 0.5))
    }
    S <- sum(1 / s2)        # v: IG(1 + n/2, 1 + eta sum tau~^2 / 2) x (1 + v S)^kpow
    ss_tau <- sum((tau - m_tau)^2)
    lv <- function(u) -(1 + n / 2) * u - (1 + eta * ss_tau / 2) * exp(-u) + kpow * log1p(exp(u) * S)
    v <- exp(slice1(log(v), lv, width = 0.5))
    list(tau = tau, xi = xi, gam = gam, s2 = s2, v = v)
  }

  for (it in 1:n_iter) {
    ## ---- stage 1 --------------------------------------------------------------
    Z <- matrix(rtnorm_side(outer(th, a) - matrix(d, n, p, byrow = TRUE), Y), n, p)
    # theta | z, items, tempered RT
    w <- 1 / s2
    prec_th <- 1 + sum(a^2) + eta * sum(gam^2 * w)
    Rth <- logT - matrix(xi, n, p, byrow = TRUE) + tau                 # = gamma theta + eps
    b_th <- drop((Z + matrix(d, n, p, byrow = TRUE)) %*% a) + eta * drop(Rth %*% (gam * w)) + m_th
    th <- b_th / prec_th + rnorm(n) / sqrt(prec_th)
    if (!is.null(Xth)) { bet_th <- draw_beta(Xth, th, 1); m_th <- drop(Xth %*% bet_th) }
    # items (a, d) with optional z-scale PX-DA
    X <- cbind(th, -1); Vad <- solve(P0ad + crossprod(X))
    if (zscale_px) for (j in 1:p) Z[, j] <- Z[, j] * zscale_move(Z[, j], X, P0ad, m0ad, Vad)
    Mad <- Vad %*% (drop(P0ad %*% m0ad) + crossprod(X, Z))
    AD <- Mad + t(chol(Vad)) %*% matrix(rnorm(2 * p), 2)
    a <- AD[1, ]; d <- AD[2, ]
    # auxiliary RT parameters psi~ under the tempered module (skipped at eta = 0: pure cut)
    if (eta > 0) {
      rb <- if (temper == "marginal") rt_block_marginal(th, tau, xi, gam, s2, v, eta, m_tau) else
        rt_block(th, tau, xi, gam, s2, v, eta, m_tau = m_tau)
      tau <- rb$tau; xi <- rb$xi; gam <- rb$gam; s2 <- rb$s2; v <- rb$v
      if (!is.null(Xtau)) { b_tau <- draw_beta(Xtau, tau, if (temper == "marginal") v / eta else v); m_tau <- drop(Xtau %*% b_tau) }
    }

    ## ---- stage 2: psi | theta, T (full likelihood), nested ---------------------
    if (K_inner > 0) {
      if (use_cpp) {
        rb2 <- rt_sweeps_cpp(th, logT, tau2, xi2, gam2, s22, v2, K_inner)
        tau2 <- rb2$tau; xi2 <- rb2$xi; gam2 <- rb2$gam; s22 <- rb2$s2; v2 <- rb2$v
      } else {
        for (k in 1:K_inner) {
          rb2 <- rt_block(th, tau2, xi2, gam2, s22, v2, 1)
          tau2 <- rb2$tau; xi2 <- rb2$xi; gam2 <- rb2$gam; s22 <- rb2$s2; v2 <- rb2$v
        }
      }
    }
    if (it > n_burn) {
      k <- it - n_burn
      if (keep_theta) TH[k, ] <- th
      A[k, ] <- a; G1[k, ] <- gam; V1[k] <- v
      if (!is.null(Xth)) B_th[k, ] <- bet_th
      if (!is.null(Xtau)) B_tau[k, ] <- b_tau
      if (eta > 0) S_tau <- S_tau + (tau - m_tau) / v
      if (K_inner > 0) { G2[k, ] <- gam2; V2[k] <- v2 }
    }
  }
  list(theta = TH, a = A, gamma_stage1 = G1, v_stage1 = V1,
       gamma_stage2 = if (K_inner > 0) G2 else NULL, v_stage2 = if (K_inner > 0) V2 else NULL,
       beta_theta = B_th, beta_tau = B_tau,
       tau_score = if (eta > 0) S_tau / n_save else NULL,   # posterior mean of (tau_i - m_i) / v
       eta = eta, K_inner = K_inner)
}
