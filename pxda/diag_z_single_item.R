if (file.exists("../pxda/probit_da_gibbs.R")) setwd("..")
# theta 固定為真值，只更新一題的 (a, d)：Albert–Chib DA 的 ESS 是多少？
suppressMessages(library(coda)); source("pxda/probit_da_gibbs.R"); source("pxda/sim_data.R")
dat <- sim_pxda(500, 15, seed = 1)
set.seed(2); th <- MASS::mvrnorm(500, c(0, 0), matrix(c(1, .3, .3, .4), 2))[, 1]
# regenerate one item from the known theta so the truth is exact
res <- sapply(c(0.8, 1.2, 2.0), function(a0) {
  y <- rbinom(500, 1, pnorm(a0 * th - 0.3))
  X <- cbind(th, -1); P0 <- diag(c(1, 1/4)); m0 <- c(1, 0)
  V <- solve(P0 + crossprod(X)); Rv <- t(chol(V))
  b <- c(1, 0); out <- matrix(NA, 5000, 2)
  for (it in 1:6000) {
    z <- rtnorm_side(drop(X %*% b), y)
    b <- drop(V %*% (P0 %*% m0 + crossprod(X, z))) + drop(Rv %*% rnorm(2))
    if (it > 1000) out[it - 1000, ] <- b
  }
  c(a_true = a0, ESS_a = unname(effectiveSize(out[, 1])), lag1 = acf(out[, 1], plot = FALSE)$acf[2])
})
print(round(res, 3))
