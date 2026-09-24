# test_px_shear.R — sampler_pxShear 的兩個檢查
#  (A) 有資料（含 cross-loading lambda）的模型：手動做一次剪切，資料 logProb 與人參數先驗都應不變
#  (B) 無資料模型（後驗 = 先驗）：MCMC 應重現精確先驗分位數；故意破壞接受率時應偏離
suppressMessages(library(nimble)); nimbleOptions(verbose = FALSE, MCMCprogressBar = FALSE)
if (file.exists("../nimble/px_sampler.R")) setwd("..")
source("nimble/px_sampler.R")

## (A) ----------------------------------------------------------------------
set.seed(1); n <- 20; p <- 5
a <- runif(p, .7, 1.5); b <- rnorm(p, 0, .5); xi <- rnorm(p, 4, .2); sig <- runif(p, .3, .6)
lam <- rnorm(p, .2, .1)
pers <- MASS::mvrnorm(n, c(0, 0), matrix(c(1, .5, .5, .4), 2))
Y <- matrix(rbinom(n * p, 1, plogis(outer(pers[, 1], a) - rep(1, n) %o% (a * b))), n)
logT <- matrix(xi, n, p, byrow = TRUE) - pers[, 2] + outer(pers[, 1], lam) +
  matrix(rnorm(n * p), n) %*% diag(sig)
codeA <- nimbleCode({
  Sigma[1:2, 1:2] ~ dinvwish(S = S0[1:2, 1:2], df = 4)
  for (j in 1:p) { a[j] ~ T(dnorm(1, sd = 1), 0, ); b[j] ~ dnorm(0, sd = 2)
                   xi[j] ~ dnorm(4, sd = 10); sig[j] ~ dinvgamma(1, 1); lam[j] ~ dnorm(0, sd = 1) }
  for (i in 1:n) { pers[i, 1:2] ~ dmnorm(mu0[1:2], cov = Sigma[1:2, 1:2])
    for (j in 1:p) { Y[i, j] ~ dbern(expit(a[j] * (pers[i, 1] - b[j])))
                     logT[i, j] ~ dnorm(xi[j] - pers[i, 2] + lam[j] * pers[i, 1], sd = sig[j]) } }
})
mA <- nimbleModel(codeA, constants = list(n = n, p = p, S0 = diag(2), mu0 = c(0, 0)),
                  data = list(Y = Y, logT = logT),
                  inits = list(Sigma = matrix(c(1.3, .4, .4, .5), 2), a = a, b = b, xi = xi,
                               sig = sig, lam = lam, pers = pers))
invisible(mA$calculate())
cc <- 0.8
dataN <- mA$getNodeNames(dataOnly = TRUE)
l0 <- c(data = mA$getLogProb(dataN), pers = mA$getLogProb("pers"),
        Sigma = mA$getLogProb("Sigma"), lam = mA$getLogProb("lam"))
L <- matrix(c(1, cc, 0, 1), 2)
mA$pers[, 2] <- mA$pers[, 2] + cc * mA$pers[, 1]
mA$lam <- mA$lam + cc
mA$Sigma <- L %*% mA$Sigma %*% t(L)
invisible(mA$calculate())
l1 <- c(data = mA$getLogProb(dataN), pers = mA$getLogProb("pers"),
        Sigma = mA$getLogProb("Sigma"), lam = mA$getLogProb("lam"))
cat("(A) logProb change after shear c = 0.8:\n"); print(signif(l1 - l0, 4))

## (B) ----------------------------------------------------------------------
nB <- 4; pB <- 3
codeB <- nimbleCode({
  Sigma[1:2, 1:2] ~ dinvwish(S = S0[1:2, 1:2], df = 4)
  for (j in 1:p) lam[j] ~ dnorm(0, sd = 1)
  for (i in 1:n) pers[i, 1:2] ~ dmnorm(mu0[1:2], cov = Sigma[1:2, 1:2])
})
mB <- nimbleModel(codeB, constants = list(n = nB, p = pB, S0 = diag(2), mu0 = c(0, 0)),
                  inits = list(Sigma = diag(2), lam = rep(0, pB), pers = matrix(0, nB, 2)))
cmB <- compileNimble(mB)
make_mcmc <- function(slope = 0) {
  conf <- configureMCMC(mB, nodes = NULL, monitors = c("Sigma", "lam", "pers"))
  for (j in 1:pB) conf$addSampler(paste0("lam[", j, "]"), "RW")
  for (i in 1:nB) conf$addSampler(paste0("pers[", i, ", 1:2]"), "RW_block")
  conf$addSampler("Sigma[1:2, 1:2]", "RW_wishart")
  conf$addSampler("Sigma[1:2, 1:2]", "pxShear",
                  control = list(sourceNodes = "pers[1:4, 1]", targetNodes = "pers[1:4, 2]",
                                 shiftNodes = "lam[1:3]", covNode = "Sigma[1:2, 1:2]",
                                 scale = 0.5, logRatioSlopeOverride = slope))
  buildMCMC(conf)
}
comp <- compileNimble(make_mcmc(0), make_mcmc(2), project = mB)
keep <- c("Sigma[1, 2]", "Sigma[2, 2]", "lam[1]", "pers[1, 2]")
run <- function(cmc) { cmc$run(400000, thin = 20, reset = TRUE)
  as.matrix(cmc$mvSamples)[-(1:2000), keep] }
set.seed(1); s_ok <- run(comp[[1]]); set.seed(2); s_bad <- run(comp[[2]])
set.seed(3); R <- 200000
Sig <- replicate(R, solve(rWishart(1, 4, diag(2))[, , 1]))
tau1 <- sapply(1:R, function(r) MASS::mvrnorm(1, c(0, 0), Sig[, , r])[2])
exact <- cbind(Sig[1, 2, ], Sig[2, 2, ], rnorm(R), tau1); colnames(exact) <- keep
qs <- c(.1, .25, .5, .75, .9)
tab <- function(s) sapply(colnames(s), function(k) quantile(s[, k], qs))
cat("\n(B) exact prior quantiles\n");            print(round(tab(exact), 3))
cat("(B) MCMC with pxShear (correct)\n");        print(round(tab(s_ok), 3))
cat("(B) MCMC with pxShear (log r + 2c, wrong)\n"); print(round(tab(s_bad), 3))
