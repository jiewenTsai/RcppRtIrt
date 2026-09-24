suppressMessages(library(aghq)); source("aghq_score.R")
set.seed(3); n <- 5; p <- 12
a <- runif(p, .7, 1.5); b <- rnorm(p, 0, .5); xi <- rnorm(p, 4, .2); sig <- runif(p, .3, .6)
S12 <- .45; S22 <- .4
th <- rnorm(n); ta <- S12 * th + rnorm(n, 0, sqrt(S22 - S12^2))
Y <- matrix(rbinom(n * p, 1, plogis(outer(th, a) - rep(1, n) %o% (a * b))), n)
logT <- matrix(xi, n, p, byrow = TRUE) - ta + matrix(rnorm(n * p), n) %*% diag(sig)
psi <- c(a, b, xi, log(sig), S12, log(S22))
mine <- marg_ll(psi, Y, logT, gh_nodes(9))
# aghq on the full 2-D (theta, tau) integrand, no analytic tau step
Sig <- matrix(c(1, S12, S12, S22), 2); Si <- solve(Sig)
pkg <- sapply(1:n, function(i) {
  fn <- function(u) sum(dbinom(Y[i, ], 1, plogis(a * (u[1] - b)), log = TRUE)) +
    sum(dnorm(logT[i, ], xi - u[2], sig, log = TRUE)) +
    (-0.5 * drop(u %*% Si %*% u) - log(2 * pi) - 0.5 * log(det(Sig)))
  gr <- function(u) numDeriv::grad(fn, u)
  he <- function(u) numDeriv::hessian(fn, u)
  q <- aghq(list(fn = fn, gr = gr, he = he), k = 9, startingvalue = c(0, 0))
  get_log_normconst(q)
})
print(rbind(mine = mine, aghq_pkg_2D = pkg, diff = mine - pkg))
