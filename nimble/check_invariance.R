# 用未編譯的 nimble 模型檢查：
#  (1) 資料節點的 logProb 在尺度變換下不變（skipInvariant 的前提）
#  (2) 接受率中 theta 的先驗變化與 Jacobian 的 +n 恰好抵消 → 接受率只剩 Sigma, a, b 的先驗
#  (3) expandNodeNames 對矩陣元素的順序
suppressMessages(library(nimble)); nimbleOptions(verbose = FALSE)
if (file.exists("../nimble/px_sampler.R")) setwd("..")
set.seed(1); n <- 20; p <- 5
a <- runif(p, .7, 1.5); b <- rnorm(p, 0, .5); xi <- rnorm(p, 4, .2); sig <- runif(p, .3, .6)
pers <- MASS::mvrnorm(n, c(0, 0), matrix(c(1, .5, .5, .4), 2))
Y <- matrix(rbinom(n * p, 1, plogis(outer(pers[, 1], a) - rep(1, n) %o% (a * b))), n)
logT <- matrix(xi, n, p, byrow = TRUE) - pers[, 2] + matrix(rnorm(n * p), n) %*% diag(sig)
code <- nimbleCode({
  Sigma[1:2, 1:2] ~ dinvwish(S = S0[1:2, 1:2], df = 4)
  for (j in 1:p) { a[j] ~ T(dnorm(1, sd = 1), 0, ); b[j] ~ dnorm(0, sd = 2)
                   xi[j] ~ dnorm(4, sd = 10); sig[j] ~ dinvgamma(1, 1) }
  for (i in 1:n) { pers[i, 1:2] ~ dmnorm(mu0[1:2], cov = Sigma[1:2, 1:2])
    for (j in 1:p) { Y[i, j] ~ dbern(expit(a[j] * (pers[i, 1] - b[j])))
                     logT[i, j] ~ dnorm(xi[j] - pers[i, 2], sd = sig[j]) } }
})
m <- nimbleModel(code, constants = list(n = n, p = p, S0 = diag(2), mu0 = c(0, 0)),
                 data = list(Y = Y, logT = logT),
                 inits = list(Sigma = matrix(c(1.3, .4, .4, .5), 2), a = a, b = b, xi = xi,
                              sig = sig, pers = pers))
m$calculate()
cat("expandNodeNames order:", m$expandNodeNames("Sigma[1:2, 1:2]", returnScalarComponents = TRUE), "\n")
u <- 0.37; alpha <- exp(u)
dataN <- m$getNodeNames(dataOnly = TRUE)
lp_data0 <- m$getLogProb(dataN); lp_pers0 <- m$getLogProb("pers"); lp_S0 <- m$getLogProb("Sigma")
lp_a0 <- m$getLogProb("a"); lp_b0 <- m$getLogProb("b")
m$pers[, 1] <- m$pers[, 1] * alpha; m$b <- m$b * alpha; m$a <- m$a / alpha
m$Sigma <- diag(c(alpha, 1)) %*% m$Sigma %*% diag(c(alpha, 1))
m$calculate()
cat(sprintf("(1) data logProb change: %.3e  (should be ~0)\n", m$getLogProb(dataN) - lp_data0))
cat(sprintf("(2) pers prior change + n*u: %.3e  (should be ~0)\n", m$getLogProb("pers") - lp_pers0 + n * u))
cat(sprintf("    remaining acceptance terms: Sigma %.4f + 3u %.4f | a %.4f - p*u %.4f | b %.4f + p*u %.4f\n",
            m$getLogProb("Sigma") - lp_S0, 3 * u, m$getLogProb("a") - lp_a0, -p * u, m$getLogProb("b") - lp_b0, p * u))
