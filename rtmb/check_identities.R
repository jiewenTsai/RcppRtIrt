# check_identities.R — numerical checks of the identities used in smi_rtmb.R (tutorial exercises).
set.seed(1); p <- 6
s2 <- runif(p, .1, .4); v <- .3; w <- 1 / s2; r <- rnorm(p)
Om <- diag(s2) + v * matrix(1, p, p)
# 1. Sherman-Morrison: r' Omega^-1 r
quad_direct <- drop(t(r) %*% solve(Om) %*% r)
quad_fast <- sum(w * r^2) - v * sum(w * r)^2 / (1 + v * sum(w))
# 2. matrix determinant lemma: log|Omega|
ld_direct <- as.numeric(determinant(Om)$modulus)
ld_fast <- sum(log(s2)) + log(1 + v * sum(w))
cat("quadratic form:", quad_direct, quad_fast, "\n")
cat("log det      :", ld_direct, ld_fast, "\n")
# 3. integrating tau out: N(T; mu - tau 1, D) N(tau; 0, v) over tau = N(T; mu, D + v 11')
mu <- rnorm(p); Tt <- mu + rnorm(p)
f <- function(tau) sapply(tau, function(t) exp(sum(dnorm(Tt, mu - t, sqrt(s2), log = TRUE)) + dnorm(t, 0, sqrt(v), log = TRUE)))
lhs <- log(integrate(f, -Inf, Inf)$value)
rhs <- -0.5 * (p * log(2 * pi) + ld_fast + sum(w * (Tt - mu)^2) - v * sum(w * (Tt - mu))^2 / (1 + v * sum(w)))
cat("marginal over tau:", lhs, rhs, "\n")
# 4. log-scale IG(1,1) prior with Jacobian: density of u = log s2 integrates to 1
cat("IG(1,1) on log scale integrates to:", integrate(function(u) exp(-u - exp(-u)), -Inf, Inf)$value, "\n")
# 5. tempering the tau prior is improper: integral over tau of N(T|tau)^eta N(tau;0,v)^eta
#    carries a factor v^((1-eta)/2), so the implied density of v does not integrate
eta <- .5
g <- function(v) sapply(v, function(vv) integrate(function(t) dnorm(t, 0, sqrt(vv))^eta, -Inf, Inf)$value)
cat("int N(tau;0,v)^.5 dtau at v = 1, 4, 16:", round(g(c(1, 4, 16)), 3), " (grows like v^(1/4))\n")
