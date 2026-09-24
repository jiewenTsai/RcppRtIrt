# check_marginal_identity.R — the DA representation used for marginal tempering
#   log ∫ prod_j N(T_j; xi_j - tau + gamma_j theta, sigma_j^2 / eta) N(tau; 0, v / eta) dtau
#     + (1 - eta)/2 log|Omega|  -  eta log N(T; xi + gamma theta, Omega)
# must not depend on the parameters (Omega = diag(sigma^2) + v 11').
if (file.exists("../smi/smi_gibbs.R")) setwd("..")
set.seed(1); p <- 6; eta <- 0.3
T <- rnorm(p, 4, .5); theta <- 0.7
lhs_minus_rhs <- function(xi, gam, s2, v) {
  f <- function(tau) sapply(tau, function(t) exp(sum(dnorm(T, xi - t + gam * theta, sqrt(s2 / eta), log = TRUE)) +
                                               dnorm(t, 0, sqrt(v / eta), log = TRUE) + 30))
  lint <- log(integrate(f, -Inf, Inf, rel.tol = 1e-12)$value) - 30
  Om <- diag(s2) + v
  mu <- xi + gam * theta
  lmarg <- -p / 2 * log(2 * pi) - 0.5 * as.numeric(determinant(Om)$modulus) -
    0.5 * drop(t(T - mu) %*% solve(Om, T - mu))
  logdetOm <- sum(log(s2)) + log(1 + v * sum(1 / s2))            # matrix determinant lemma
  c(diff = lint + (1 - eta) / 2 * logdetOm - eta * lmarg,
    lemma_err = logdetOm - as.numeric(determinant(Om)$modulus))
}
res <- t(replicate(5, lhs_minus_rhs(rnorm(p, 4, .3), rnorm(p, .5, .2), runif(p, .1, .5), runif(1, .05, .5))))
print(signif(res, 8))
cat("spread of diff across parameter draws:", signif(diff(range(res[, "diff"])), 3), "\n")
