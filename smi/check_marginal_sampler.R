# check_marginal_sampler.R — the DA sampler for marginal tempering (temper = "marginal") vs an
# independent sampler of the same stage-1 target that never introduces tau~: theta is drawn from its
# Gaussian conditional under eta * log N(T_i; xi + gamma theta_i, Omega), and (xi, gamma, log sigma^2,
# log v) by random-walk MH on the exact tempered marginal likelihood. Quantiles must agree.
if (file.exists("../smi/smi_gibbs.R")) setwd("..")
source("pxda/probit_da_gibbs.R"); source("smi/smi_gibbs.R"); source("smi/sim_smi.R")
dat <- sim_smi(80, 4, seed = 5, gam_A = .5, gam_B = .5, v = .15)
Y <- dat$Y; logT <- dat$logT; n <- nrow(Y); p <- ncol(Y); eta <- 0.3
n_iter <- 60000; n_burn <- 5000

direct_mh <- function(seed) {
  set.seed(seed)
  a <- rep(1, p); d <- rep(0, p); th <- rep(0, n)
  xi <- colMeans(logT); gam <- rep(0, p); ls2 <- rep(log(.25), p); lv <- 0
  P0ad <- diag(c(1, 1 / 4)); m0ad <- c(1, 0)
  lp_psi <- function(xi, gam, ls2, lv, th) {
    s2 <- exp(ls2); v <- exp(lv); Om <- diag(s2) + v
    Oi <- solve(Om); ld <- as.numeric(determinant(Om)$modulus)
    R <- logT - matrix(xi, n, p, byrow = TRUE) - outer(th, gam)
    ll <- -n / 2 * ld - 0.5 * sum((R %*% Oi) * R)
    lpri <- sum(dnorm(xi, 4, 10, log = TRUE)) + sum(dnorm(gam, 0, 1, log = TRUE)) +
      sum(-2 * ls2 - exp(-ls2) + ls2) + (-2 * lv - exp(-lv) + lv)     # IG(1,1) on s2, v with log Jacobian
    eta * ll + lpri
  }
  out <- matrix(NA, n_iter - n_burn, 6, dimnames = list(NULL, c("v", "gam1", "s2_1", "th1", "th2", "a1")))
  cur <- lp_psi(xi, gam, ls2, lv, th); sd_rw <- c(.05, .05, .15, .15)
  for (it in 1:n_iter) {
    Z <- matrix(rtnorm_side(outer(th, a) - matrix(d, n, p, byrow = TRUE), Y), n, p)
    s2 <- exp(ls2); v <- exp(lv); Oi <- solve(diag(s2) + v)
    prec <- 1 + sum(a^2) + eta * drop(t(gam) %*% Oi %*% gam)
    b <- drop((Z + matrix(d, n, p, byrow = TRUE)) %*% a) +
      eta * drop((logT - matrix(xi, n, p, byrow = TRUE)) %*% Oi %*% gam)
    th <- b / prec + rnorm(n) / sqrt(prec)
    X <- cbind(th, -1); Vad <- solve(P0ad + crossprod(X))
    AD <- Vad %*% (drop(P0ad %*% m0ad) + crossprod(X, Z)) + t(chol(Vad)) %*% matrix(rnorm(2 * p), 2)
    a <- AD[1, ]; d <- AD[2, ]
    cur <- lp_psi(xi, gam, ls2, lv, th)
    for (blk in 1:4) {                                  # xi, gamma, log sigma^2, log v
      xi1 <- xi; gam1 <- gam; ls21 <- ls2; lv1 <- lv
      if (blk == 1) xi1 <- xi + rnorm(p, 0, sd_rw[1])
      if (blk == 2) gam1 <- gam + rnorm(p, 0, sd_rw[2])
      if (blk == 3) ls21 <- ls2 + rnorm(p, 0, sd_rw[3])
      if (blk == 4) lv1 <- lv + rnorm(1, 0, sd_rw[4])
      new <- lp_psi(xi1, gam1, ls21, lv1, th)
      if (log(runif(1)) < new - cur) { xi <- xi1; gam <- gam1; ls2 <- ls21; lv <- lv1; cur <- new }
    }
    if (it > n_burn) out[it - n_burn, ] <- c(exp(lv), gam[1], exp(ls2[1]), th[1], th[2], a[1])
  }
  out
}
t0 <- proc.time()[3]
S_mh <- direct_mh(1)
f <- smi_rtirt(Y, logT, eta = eta, n_iter = n_iter, n_burn = n_burn, temper = "marginal",
               zscale_px = FALSE, seed = 2)
S_da <- cbind(v = f$v_stage1, gam1 = f$gamma_stage1[, 1], th1 = f$theta[, 1], th2 = f$theta[, 2], a1 = f$a[, 1])
qs <- c(.05, .25, .5, .75, .95)
cmp <- function(k) rbind(direct_MH = quantile(S_mh[, k], qs), DA_marginal = quantile(S_da[, k], qs))
for (k in c("v", "gam1", "th1", "th2", "a1")) { cat("\n", k, "\n"); print(round(cmp(k), 3)) }
cat(sprintf("\n%.0f s\n", proc.time()[3] - t0))
