if (file.exists("../aghq_score.R")) setwd(".."); source("aghq_score.R")
set.seed(1); n <- 5; p <- 10
a <- runif(p, .7, 1.5); b <- rnorm(p, 0, .5); xi <- rnorm(p, 4, .2); sig <- runif(p, .3, .6)
S12 <- 0.5; S22 <- 0.4
th <- rnorm(n); ta <- S12 * th + rnorm(n, 0, sqrt(S22 - S12^2))
Y <- matrix(rbinom(n * p, 1, plogis(outer(th, a) - rep(1, n) %o% (a * b))), n)
logT <- matrix(xi, n, p, byrow = TRUE) - ta + matrix(rnorm(n * p), n) %*% diag(sig)
psi <- c(a, b, xi, log(sig), S12, log(S22))

# brute force: 2-D grid over (theta, tau), includes full RT density
brute <- function(i) {
  g <- seq(-7, 7, length.out = 1401); h <- diff(g)[1]
  Sig <- matrix(c(1, S12, S12, S22), 2)
  f <- outer(g, g, function(t1, t2) {
    ly <- sapply(seq_along(t1), function(m) sum(dbinom(Y[i, ], 1, plogis(a * (t1[m] - b)), log = TRUE)))
    lt <- sapply(seq_along(t2), function(m) sum(dnorm(logT[i, ], xi - t2[m], sig, log = TRUE)))
    q <- cbind(t1, t2) %*% solve(Sig)
    lp <- -0.5 * rowSums(q * cbind(t1, t2)) - log(2 * pi) - 0.5 * log(det(Sig))
    ly + lt + lp
  })
  m <- max(f); m + log(sum(exp(f - m)) * h^2)
}
bf <- sapply(1:n, brute)
res <- sapply(c(1, 3, 5, 9, 15), function(K) marg_ll(psi, Y, logT, gh_nodes(K)) - bf)
colnames(res) <- paste0("K=", c(1, 3, 5, 9, 15))
cat("AGHQ minus brute-force 2-D grid (per person):\n"); print(signif(res, 3))
