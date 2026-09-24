# probit_da_gibbs.R — probit (Albert–Chib) DA Gibbs for the RT-IRT model with cross-loadings,
# with optional PX group moves (scale, shear). Pure R, vectorised over persons.
#
# Model (expanded parameterisation, method = "expanded"):
#   z_ij ~ N(a_j theta_i - d_j, 1),  Y_ij = 1{z_ij > 0}
#   log T_ij ~ N(xi_j - tau_i + lambda_j theta_i, sigma_j^2)
#   (theta_i, tau_i) ~ N(0, Sigma),  Sigma ~ IW(nu0, I)
#   a_j ~ N(1, 1), d_j ~ N(0, 2^2), xi_j ~ N(4, 10^2), lambda_j ~ N(0, 1), sigma_j^2 ~ IG(1, 1)
# The scale of theta and the common part of lambda are not identified; report the identified
# functions in identified_draws().
#
# method = "identified": Sigma_11 = 1 is fixed; tau | theta ~ N(beta theta, v) with
#   beta ~ N(0, 1), v ~ IG(1.5, 0.5), both conjugate. This is a different prior on the
#   identified quantities than the expanded model induces.
#
# Group moves (expanded only), each a symmetric random-walk MH on the group element:
#   scale: theta, Sigma row/col 1 × alpha; a, lambda / alpha; d unchanged
#          Jacobian exponent n - 2p + 3
#   shear: tau += c theta, lambda += c, Sigma -> L Sigma L^T; Jacobian 1
# z and log T are unchanged by both moves, so only prior terms enter the acceptance ratio.

rtnorm_side <- function(mu, y) {           # z ~ N(mu, 1) truncated to z > 0 (y = 1) or z < 0
  # inverse CDF on the log scale, so extreme mu does not underflow to +-Inf
  lu <- log(runif(length(mu)))
  s <- ifelse(y == 1, 1, -1)
  mu - s * qnorm(lu + pnorm(s * mu, log.p = TRUE), log.p = TRUE)
}

slice1 <- function(x0, logf, width = 1, max_steps = 50) {   # univariate slice sampler (Neal 2003)
  ly <- logf(x0) - rexp(1)
  L <- x0 - runif(1) * width; U <- L + width
  k <- max_steps
  while (k > 0 && logf(L) > ly) { L <- L - width; k <- k - 1 }
  k <- max_steps
  while (k > 0 && logf(U) > ly) { U <- U + width; k <- k - 1 }
  repeat {
    x1 <- runif(1, L, U)
    if (logf(x1) > ly) return(x1)
    if (x1 < x0) L <- x1 else U <- x1
  }
}

# z-scale PX-DA move for one probit regression block z ~ N(X beta, I), beta ~ N(m0, P0^{-1}):
# z -> g z along the orbit, with beta integrated out. Target along the orbit (w.r.t. dg):
#   pi(g) ∝ g^(n-1) exp(-(A g^2 - 2 B g) / 2),
#   A = z'z - z'X V X'z,  B = z'c - z'X V X'c,  c = X m0,  V = (P0 + X'X)^{-1}.
# Independence proposal g^2 ~ Gamma(n/2, A/2) is exact when B = 0; MH weight exp(B (g - 1)),
# since the current state is g = 1. sign(g z) = sign(z), so the data are unaffected.
zscale_move <- function(z, X, P0, m0, V) {
  n <- length(z); cc <- drop(X %*% m0)
  Xz <- crossprod(X, z); Xc <- crossprod(X, cc)
  A <- sum(z^2) - drop(t(Xz) %*% V %*% Xz)
  B <- sum(z * cc) - drop(t(Xz) %*% V %*% Xc)
  if (!is.finite(A) || !is.finite(B) || A <= 0) return(1)
  g_new <- sqrt(rgamma(1, n / 2, A / 2))
  if (log(runif(1)) < B * (g_new - 1)) g_new else 1
}

riwish <- function(df, S) solve(rWishart(1, df, solve(S))[, , 1])

log_iw <- function(Sig, df, Psi) {         # IW log density up to a constant
  d <- nrow(Sig)
  -(df + d + 1) / 2 * as.numeric(determinant(Sig)$modulus) - 0.5 * sum(diag(Psi %*% solve(Sig)))
}

log_mvn0 <- function(P, Sig) {             # sum_i log N(P[i, ]; 0, Sig) up to a constant
  n <- nrow(P)
  -n / 2 * as.numeric(determinant(Sig)$modulus) - 0.5 * sum((P %*% solve(Sig)) * P)
}

probit_da_gibbs <- function(Y, logT, n_iter = 6000, n_burn = 1000,
                            method = c("expanded", "identified"),
                            scale_move = FALSE, shear_move = FALSE,
                            collapse_items = FALSE, zscale_px = FALSE,
                            nu0 = 4, seed = 1) {
  method <- match.arg(method)
  if (method == "identified" && (scale_move || shear_move))
    stop("group moves are defined on the expanded model only")
  set.seed(seed)
  n <- nrow(Y); p <- ncol(Y)
  # initial values
  a <- rep(1, p); d <- rep(0, p); xi <- colMeans(logT); lam <- rep(0, p); s2 <- rep(0.25, p)
  pers <- matrix(0, n, 2)
  Sig <- diag(2); beta <- 0; v <- 1
  su <- 0.05; sc <- 0.05                   # proposal sds, adapted during burn-in
  acc_u <- acc_c <- 0
  out <- matrix(NA_real_, n_iter - n_burn, 4 + 2 * p,
                dimnames = list(NULL, c("S11", "S12", "S22", "lbar",
                                        paste0("a", 1:p), paste0("lam", 1:p))))
  P0ad <- diag(c(1, 1 / 4)); m0ad <- c(1, 0)
  P0xl <- diag(c(1 / 100, 1)); m0xl <- c(4, 0)

  for (it in 1:n_iter) {
    # 1. latent z
    mu <- outer(pers[, 1], a) - matrix(d, n, p, byrow = TRUE)
    Z <- matrix(rtnorm_side(mu, Y), n, p)

    # 2-3. persons and items (a_j, d_j)
    R <- logT - matrix(xi, n, p, byrow = TRUE)                 # = lam theta - tau + e
    w <- 1 / s2
    Prt <- matrix(c(sum(lam^2 * w), -sum(lam * w), -sum(lam * w), sum(w)), 2)   # RT part
    Brt <- cbind(R %*% (lam * w), -R %*% w)
    if (collapse_items) {
      # Partially collapsed: draw each (a_j, d_j) from p(a_j, d_j | z, log T, other items)
      # with persons integrated out, then persons from their full conditional.
      Zd <- Z + matrix(d, n, p, byrow = TRUE)
      P <- solve(Sig) + Prt + matrix(c(sum(a^2), 0, 0, 0), 2)
      B <- Brt + cbind(Zd %*% a, 0)
      for (j in 1:p) {
        Pm <- P - matrix(c(a[j]^2, 0, 0, 0), 2)                 # leave item j out
        Bm <- B; Bm[, 1] <- Bm[, 1] - a[j] * Zd[, j]
        Vm <- solve(Pm); m1 <- drop(Bm %*% Vm[, 1]); v11 <- Vm[1, 1]
        zj <- Z[, j]
        logpost_a <- function(aa) {                           # d_j integrated out, prior N(0, 4)
          s2a <- 1 + aa^2 * v11; u <- aa * m1 - zj; ub <- mean(u)
          dnorm(aa, 1, 1, log = TRUE) - n / 2 * log(s2a) - sum((u - ub)^2) / (2 * s2a) +
            0.5 * log(s2a / n) + dnorm(ub, 0, sqrt(4 + s2a / n), log = TRUE)
        }
        a[j] <- slice1(a[j], logpost_a, width = 0.5)
        s2a <- 1 + a[j]^2 * v11; u <- a[j] * m1 - zj
        prec <- n / s2a + 1 / 4
        d[j] <- rnorm(1, sum(u) / s2a / prec, sqrt(1 / prec))
        Zd[, j] <- Z[, j] + d[j]
        P <- Pm + matrix(c(a[j]^2, 0, 0, 0), 2)
        B <- Bm; B[, 1] <- B[, 1] + a[j] * Zd[, j]
      }
      V <- solve(P)
      pers <- B %*% V + matrix(rnorm(2 * n), n) %*% chol(V)
    } else {
      A <- Prt + matrix(c(sum(a^2), 0, 0, 0), 2)
      V <- solve(solve(Sig) + A)
      Bm <- Brt + cbind((Z + matrix(d, n, p, byrow = TRUE)) %*% a, 0)
      pers <- Bm %*% V + matrix(rnorm(2 * n), n) %*% chol(V)
      th <- pers[, 1]
      X <- cbind(th, -1)
      Vad <- solve(P0ad + crossprod(X))
      if (zscale_px) {                 # Liu-Wu PX-DA on each item's latent scale (z_j, a_j, d_j)
        for (j in 1:p) Z[, j] <- Z[, j] * zscale_move(Z[, j], X, P0ad, m0ad, Vad)
      }
      Mad <- Vad %*% (drop(P0ad %*% m0ad) + crossprod(X, Z))
      AD <- Mad + t(chol(Vad)) %*% matrix(rnorm(2 * p), 2)
      a <- AD[1, ]; d <- AD[2, ]
    }
    th <- pers[, 1]; ta <- pers[, 2]

    # 4. (xi_j, lambda_j): log T_j + tau = xi_j + lambda_j theta + e
    W <- cbind(1, th); WtW <- crossprod(W); Wty <- crossprod(W, logT + ta)
    for (j in 1:p) {
      Vj <- solve(P0xl + WtW / s2[j])
      mj <- Vj %*% (P0xl %*% m0xl + Wty[, j] / s2[j])
      xl <- mj + t(chol(Vj)) %*% rnorm(2)
      xi[j] <- xl[1]; lam[j] <- xl[2]
    }

    # 5. sigma_j^2
    E <- logT - (matrix(xi, n, p, byrow = TRUE) - ta + outer(th, lam))
    s2 <- 1 / rgamma(p, 1 + n / 2, 1 + colSums(E^2) / 2)

    # 6. person covariance
    if (method == "expanded") {
      Sig <- riwish(nu0 + n, diag(2) + crossprod(pers))
    } else {
      vb <- 1 / (1 + sum(th^2) / v); beta <- rnorm(1, vb * sum(th * ta) / v, sqrt(vb))
      v <- 1 / rgamma(1, 1.5 + n / 2, 0.5 + sum((ta - beta * th)^2) / 2)
      Sig <- matrix(c(1, beta, beta, beta^2 + v), 2)
    }

    # 7. scale move: theta, Sigma row/col 1 × alpha; a, lambda / alpha
    if (scale_move) {
      u <- rnorm(1, 0, su); al <- exp(u)
      D <- diag(c(al, 1))
      pers1 <- pers; pers1[, 1] <- pers1[, 1] * al
      Sig1 <- D %*% Sig %*% D; a1 <- a / al; lam1 <- lam / al
      lr <- (log_mvn0(pers1, Sig1) + log_iw(Sig1, nu0, diag(2)) +
               sum(dnorm(a1, 1, 1, log = TRUE)) + sum(dnorm(lam1, 0, 1, log = TRUE))) -
            (log_mvn0(pers, Sig) + log_iw(Sig, nu0, diag(2)) +
               sum(dnorm(a, 1, 1, log = TRUE)) + sum(dnorm(lam, 0, 1, log = TRUE))) +
            (n - 2 * p + 3) * u
      if (log(runif(1)) < lr) {
        pers <- pers1; Sig <- Sig1; a <- a1; lam <- lam1; acc_u <- acc_u + 1
      }
    }

    # 8. shear move: tau += c theta, lambda += c, Sigma -> L Sigma L^T
    if (shear_move) {
      cc <- rnorm(1, 0, sc)
      L <- matrix(c(1, cc, 0, 1), 2)
      pers1 <- pers; pers1[, 2] <- pers1[, 2] + cc * pers1[, 1]
      Sig1 <- L %*% Sig %*% t(L); lam1 <- lam + cc
      lr <- (log_mvn0(pers1, Sig1) + log_iw(Sig1, nu0, diag(2)) +
               sum(dnorm(lam1, 0, 1, log = TRUE))) -
            (log_mvn0(pers, Sig) + log_iw(Sig, nu0, diag(2)) +
               sum(dnorm(lam, 0, 1, log = TRUE)))
      if (log(runif(1)) < lr) {
        pers <- pers1; Sig <- Sig1; lam <- lam1; acc_c <- acc_c + 1
      }
    }

    # adapt proposal sds during burn-in only (then fixed, so the saved chain is time-homogeneous)
    if (it <= n_burn && it %% 50 == 0) {
      if (scale_move) { su <- su * exp(acc_u / 50 - 0.44); acc_u <- 0 }
      if (shear_move) { sc <- sc * exp(acc_c / 50 - 0.44); acc_c <- 0 }
    }
    if (it > n_burn) {
      out[it - n_burn, ] <- c(Sig[1, 1], Sig[1, 2], Sig[2, 2], mean(lam), a, lam)
    }
  }
  attr(out, "accept") <- c(scale = acc_u / (n_iter - n_burn), shear = acc_c / (n_iter - n_burn))
  out
}

# Functions of the draws that are invariant to both the scale and the shear
identified_draws <- function(S) {
  p <- (ncol(S) - 4) / 2
  s11 <- S[, "S11"]; s12 <- S[, "S12"]; s22 <- S[, "S22"]; lb <- S[, "lbar"]
  A <- S[, paste0("a", 1:p)]; L <- S[, paste0("lam", 1:p)]
  kappa <- (s12 - lb * s11) / sqrt(s11)                 # cov(theta*, tau - lbar theta)
  vstar <- s22 - 2 * lb * s12 + lb^2 * s11              # var(tau - lbar theta)
  list(summary = cbind(kappa = kappa, vstar = vstar, rho_star = kappa / sqrt(vstar)),
       a_id = A * sqrt(s11), lamdev_id = (L - lb) * sqrt(s11),
       raw = cbind(S11 = s11, S12 = s12, lbar = lb))
}
