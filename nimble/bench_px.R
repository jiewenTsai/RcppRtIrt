# bench_px.R — 同一個擴張（未識別）模型、同一組預設 sampler，比較有無 pxScale
suppressMessages({library(nimble); library(coda)})
nimbleOptions(verbose = FALSE, MCMCprogressBar = FALSE)
if (file.exists("../nimble/px_sampler.R")) setwd("..")
source("nimble/px_sampler.R")

args <- commandArgs(TRUE)
n <- as.integer(if (length(args) >= 1) args[1] else 300)
p <- as.integer(if (length(args) >= 2) args[2] else 15)
niter <- 6000; nburn <- 1000
set.seed(42)
a <- runif(p, .7, 1.5); b <- rnorm(p, 0, .5); xi <- rnorm(p, 4, .2); sig <- runif(p, .3, .6)
S <- matrix(c(1, .5, .5, .4), 2)
pers <- MASS::mvrnorm(n, c(0, 0), S)
Y <- matrix(rbinom(n * p, 1, plogis(outer(pers[, 1], a) - rep(1, n) %o% (a * b))), n)
logT <- matrix(xi, n, p, byrow = TRUE) - pers[, 2] + matrix(rnorm(n * p), n) %*% diag(sig)

code <- nimbleCode({
  Sigma[1:2, 1:2] ~ dinvwish(S = S0[1:2, 1:2], df = 4)
  for (j in 1:p) {
    a[j] ~ T(dnorm(1, sd = 1), 0, )
    b[j] ~ dnorm(0, sd = 2)
    xi[j] ~ dnorm(4, sd = 10)
    sig[j] ~ dinvgamma(1, 1)
  }
  for (i in 1:n) {
    pers[i, 1:2] ~ dmnorm(mu0[1:2], cov = Sigma[1:2, 1:2])
    for (j in 1:p) {
      Y[i, j] ~ dbern(expit(a[j] * (pers[i, 1] - b[j])))
      logT[i, j] ~ dnorm(xi[j] - pers[i, 2], sd = sig[j])
    }
  }
})
m <- nimbleModel(code, constants = list(n = n, p = p, S0 = diag(2), mu0 = c(0, 0)),
                 data = list(Y = Y, logT = logT),
                 inits = list(Sigma = diag(2), a = rep(1, p), b = rep(0, p), xi = colMeans(logT),
                              sig = rep(.5, p), pers = matrix(0, n, 2)))
cm <- compileNimble(m)
mon <- c("Sigma", "a", "b")
conf0 <- configureMCMC(m, monitors = mon)
cat("default samplers:\n"); print(table(sapply(conf0$getSamplers(), function(s) s$name)))
conf1 <- configureMCMC(m, monitors = mon)
conf1$addSampler("Sigma[1:2, 1:2]", "pxScale",
                 control = list(scaleUp = c(sprintf("pers[1:%d, 1]", n), sprintf("b[1:%d]", p)),
                                scaleDown = sprintf("a[1:%d]", p), covNode = "Sigma[1:2, 1:2]",
                                skipInvariant = TRUE, scale = 0.2))
comp <- compileNimble(buildMCMC(conf0), buildMCMC(conf1), project = m)

summ <- function(cmc, label) {
  t0 <- proc.time()[3]; cmc$run(niter); secs <- proc.time()[3] - t0
  s <- as.matrix(cmc$mvSamples)[-(1:nburn), ]
  s11 <- s[, "Sigma[1, 1]"]; s12 <- s[, "Sigma[1, 2]"]; s22 <- s[, "Sigma[2, 2]"]
  aid <- s[, grep("^a\\[", colnames(s))] * sqrt(s11)
  q <- cbind(S11_raw = s11, S12_id = s12 / sqrt(s11), rho = s12 / sqrt(s11 * s22),
             a_id_med = apply(aid, 1, median))
  ess <- effectiveSize(mcmc(q)); ess_a <- median(effectiveSize(mcmc(aid)))
  cat(sprintf("\n== %s: %.1f s for %d iter\n", label, secs, niter))
  print(round(rbind(post_mean = colMeans(q), ESS = ess, ESS_per_s = ess / secs), 3))
  cat(sprintf("median ESS over a_id[j]: %.0f  (%.1f per s)\n", ess_a, ess_a / secs))
  cat(sprintf("truth: S12_id = .5, rho = %.3f, median a = %.3f\n", .5 / sqrt(.4), median(a)))
}
set.seed(7); summ(comp[[1]], "default")
set.seed(7); summ(comp[[2]], "default + pxScale")
