if (file.exists("../pxda/probit_da_gibbs.R")) setwd("..")
suppressMessages(library(coda)); source("pxda/probit_da_gibbs.R")
# z-scale PX step for one probit regression block: z -> g z, with beta integrated out
#   log p(g | z) = (n - 1) log g - (A g^2 - 2 B g) / 2  (Haar dg/g, Gaussian prior on beta)
set.seed(2); th <- MASS::mvrnorm(500, c(0, 0), matrix(c(1, .3, .3, .4), 2))[, 1]
run <- function(a0, px) {
  set.seed(3); y <- rbinom(500, 1, pnorm(a0 * th - 0.3))
  X <- cbind(th, -1); P0 <- diag(c(1, 1/4)); m0 <- c(1, 0)
  V <- solve(P0 + crossprod(X)); Rv <- t(chol(V))
  b <- c(1, 0); out <- matrix(NA, 5000, 2)
  for (it in 1:6000) {
    z <- rtnorm_side(drop(X %*% b), y)
    if (px) z <- z * zscale_move(z, X, P0, m0, V)
    b <- drop(V %*% (P0 %*% m0 + crossprod(X, z))) + drop(Rv %*% rnorm(2))
    if (it > 1000) out[it - 1000, ] <- b
  }
  c(ESS_a = unname(effectiveSize(out[, 1])), mean_a = mean(out[, 1]), sd_a = sd(out[, 1]))
}
for (a0 in c(0.8, 1.2, 2.0)) {
  r0 <- run(a0, FALSE); r1 <- run(a0, TRUE)
  cat(sprintf("a = %.1f | DA: ESS %5.0f mean %.3f sd %.3f | DA + z-scale PX: ESS %5.0f mean %.3f sd %.3f\n",
              a0, r0[1], r0[2], r0[3], r1[1], r1[2], r1[3]))
}
