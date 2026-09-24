# bench_px_shear.R — 含 cross-loading lambda 的擴張模型，預設 sampler vs + pxShear（+ pxScale）
suppressMessages({library(nimble); library(coda)})
nimbleOptions(verbose = FALSE, MCMCprogressBar = FALSE)
if (file.exists("../nimble/px_sampler.R")) setwd("..")
source("nimble/px_sampler.R")
n <- 300; p <- 10; niter <- 6000; nburn <- 1000
set.seed(11)
a <- runif(p, .7, 1.5); b <- rnorm(p, 0, .5); xi <- rnorm(p, 4, .2); sig <- runif(p, .3, .6)
lam <- rnorm(p, .15, .1)
pers <- MASS::mvrnorm(n, c(0, 0), matrix(c(1, .3, .3, .4), 2))
Y <- matrix(rbinom(n * p, 1, plogis(outer(pers[, 1], a) - rep(1, n) %o% (a * b))), n)
logT <- matrix(xi, n, p, byrow = TRUE) - pers[, 2] + outer(pers[, 1], lam) +
  matrix(rnorm(n * p), n) %*% diag(sig)
code <- nimbleCode({
  Sigma[1:2, 1:2] ~ dinvwish(S = S0[1:2, 1:2], df = 4)
  for (j in 1:p) { a[j] ~ T(dnorm(1, sd = 1), 0, ); b[j] ~ dnorm(0, sd = 2)
                   xi[j] ~ dnorm(4, sd = 10); sig[j] ~ dinvgamma(1, 1); lam[j] ~ dnorm(0, sd = 1) }
  for (i in 1:n) { pers[i, 1:2] ~ dmnorm(mu0[1:2], cov = Sigma[1:2, 1:2])
    for (j in 1:p) { Y[i, j] ~ dbern(expit(a[j] * (pers[i, 1] - b[j])))
                     logT[i, j] ~ dnorm(xi[j] - pers[i, 2] + lam[j] * pers[i, 1], sd = sig[j]) } }
})
m <- nimbleModel(code, constants = list(n = n, p = p, S0 = diag(2), mu0 = c(0, 0)),
                 data = list(Y = Y, logT = logT),
                 inits = list(Sigma = diag(2), a = rep(1, p), b = rep(0, p), xi = colMeans(logT),
                              sig = rep(.5, p), lam = rep(0, p), pers = matrix(0, n, 2)))
cm <- compileNimble(m)
mon <- c("Sigma", "lam", "a")
shear_ctl <- list(sourceNodes = sprintf("pers[1:%d, 1]", n), targetNodes = sprintf("pers[1:%d, 2]", n),
                  shiftNodes = sprintf("lam[1:%d]", p), covNode = "Sigma[1:2, 1:2]",
                  skipInvariant = TRUE, scale = 0.1)
scale_ctl <- list(scaleUp = c(sprintf("pers[1:%d, 1]", n), sprintf("b[1:%d]", p)),
                  scaleDown = c(sprintf("a[1:%d]", p), sprintf("lam[1:%d]", p)),
                  covNode = "Sigma[1:2, 1:2]", skipInvariant = TRUE, scale = 0.2)
conf0 <- configureMCMC(m, monitors = mon)
conf1 <- configureMCMC(m, monitors = mon)
conf1$addSampler("Sigma[1:2, 1:2]", "pxShear", control = shear_ctl)
conf2 <- configureMCMC(m, monitors = mon)
conf2$addSampler("Sigma[1:2, 1:2]", "pxShear", control = shear_ctl)
conf2$addSampler("Sigma[1:2, 1:2]", "pxScale", control = scale_ctl)
comp <- compileNimble(buildMCMC(conf0), buildMCMC(conf1), buildMCMC(conf2), project = m)

summ <- function(cmc, label) {
  t0 <- proc.time()[3]; cmc$run(niter); secs <- proc.time()[3] - t0
  s <- as.matrix(cmc$mvSamples)[-(1:nburn), ]
  s11 <- s[, "Sigma[1, 1]"]; s12 <- s[, "Sigma[1, 2]"]
  L <- s[, grep("^lam\\[", colnames(s))]; lbar <- rowMeans(L)
  kappa <- (s12 - lbar * s11) / sqrt(s11)
  dev <- (L - lbar) * sqrt(s11)
  q <- cbind(lbar_raw = lbar, S12_raw = s12, S11_raw = s11, kappa_id = kappa)
  ess <- effectiveSize(mcmc(q)); ess_dev <- median(effectiveSize(mcmc(dev)))
  aid <- s[, grep("^a\\[", colnames(s))] * sqrt(s11)
  ess_a <- median(effectiveSize(mcmc(aid)))
  cat(sprintf("\n== %s: %.1f s\n", label, secs))
  print(round(rbind(post_mean = colMeans(q), ESS = ess, ESS_per_s = ess / secs), 3))
  cat(sprintf("median ESS (lambda_j - lbar) * sqrt(S11): %.0f | a_j * sqrt(S11): %.0f\n", ess_dev, ess_a))
}
cat(sprintf("truth: kappa = (S12 - lbar S11)/sqrt(S11) = %.3f\n", .3 - mean(lam)))
set.seed(7); summ(comp[[1]], "default")
set.seed(7); summ(comp[[2]], "default + pxShear")
set.seed(7); summ(comp[[3]], "default + pxShear + pxScale")
