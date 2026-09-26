# smi_rtmb.R — stage 1 of SMI RT-IRT with RTMB (Laplace approximation over theta).
# Same model and priors as smi/smi_gibbs.R with temper = "marginal":
#   Y_ij = 1{a_j theta_i - d_j + e > 0},  theta_i ~ N(X_i b_th, 1)
#   log T_i ~ N(xi - X_i b_tau + gamma theta_i, Omega),  Omega = diag(s^2) + v 11'   (tau integrated out)
# The objective is  -log p(Y | .) - eta * log L_RT(.) - log p(theta) - log priors:  the RT marginal
# likelihood is raised to eta directly, so no augmentation or |Omega| correction is needed.
# Fixed effects are estimated at the posterior mode (priors as penalties, log-scale variances with
# their Jacobians); theta is integrated by the Laplace approximation. eta = 0 maps the RT parameters
# out (the cut). Targets are linear functionals W' theta, reported with sdreport().
suppressMessages(library(RTMB))
smi_rtmb <- function(Y, logT, eta, X = NULL, W, silent = TRUE) {
  n <- nrow(Y); p <- ncol(Y)
  X <- if (is.null(X)) matrix(0, n, 1) else as.matrix(X)
  q <- ncol(X)
  S <- 2 * Y - 1
  par <- list(theta = rep(0, n), a = rep(1, p), d = rep(0, p), b_th = rep(0, q),
              xi = colMeans(logT), gam = rep(0, p), ls2 = rep(log(.25), p), lv = log(.3),
              b_tau = rep(0, q))
  f <- function(par) {
    getAll(par)
    lin <- matrix(theta, n, p) * matrix(a, n, p, byrow = TRUE) - matrix(d, n, p, byrow = TRUE)
    nll <- -sum(log(pnorm(S * lin)))
    nll <- nll - sum(dnorm(theta, drop(X %*% b_th), 1, log = TRUE))
    nll <- nll - sum(dnorm(a, 1, 1, log = TRUE)) - sum(dnorm(d, 0, 2, log = TRUE)) -
      sum(dnorm(b_th, 0, sqrt(10), log = TRUE))
    if (eta > 0) {
      w <- exp(-ls2); v <- exp(lv); Sw <- sum(w)
      m <- drop(X %*% b_tau)
      R <- logT - matrix(xi, n, p, byrow = TRUE) + matrix(m, n, p) -
        matrix(theta, n, p) * matrix(gam, n, p, byrow = TRUE)
      rw <- drop(R %*% w)
      quad <- drop((R * R) %*% w) - v * rw^2 / (1 + v * Sw)
      ll <- -0.5 * (p * log(2 * pi) + sum(ls2) + log(1 + v * Sw) + quad)
      nll <- nll - eta * sum(ll)
      # priors: xi ~ N(4, 10^2), gamma ~ N(0, 1), s^2 ~ IG(1, 1), v ~ IG(1, 1) (on the log scale)
      nll <- nll - sum(dnorm(xi, 4, 10, log = TRUE)) - sum(dnorm(gam, 0, 1, log = TRUE)) -
        sum(-ls2 - exp(-ls2)) - (-lv - exp(-lv)) - sum(dnorm(b_tau, 0, sqrt(10), log = TRUE))
    }
    tgt <- drop(t(W) %*% theta)
    ADREPORT(tgt)
    nll
  }
  map <- if (eta > 0) list() else
    list(xi = factor(rep(NA, p)), gam = factor(rep(NA, p)), ls2 = factor(rep(NA, p)),
         lv = factor(NA), b_tau = factor(rep(NA, q)))
  obj <- MakeADFun(f, par, random = "theta", map = map, silent = silent)
  opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(iter.max = 500, eval.max = 1000))
  sdr <- sdreport(obj)
  rep <- summary(sdr, "report")
  list(est = rep[, 1], se = rep[, 2], theta = summary(sdr, "random")[, 1], opt = opt,
       conv = opt$convergence, eta = eta)
}

# weights for a gap (binary b) and an OLS slope (continuous z) as linear functionals of theta
w_gap <- function(b) ifelse(b == 1, 1 / sum(b == 1), -1 / sum(b == 0))
w_slope <- function(z) { zc <- z - mean(z); zc / sum(zc^2) }
