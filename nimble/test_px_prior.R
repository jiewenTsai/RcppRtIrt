# test_px_prior.R — 無資料模型（後驗 = 先驗）下檢查 sampler_pxScale 保持目標分布不變
#   正確 Jacobian 的鏈應重現精確的先驗分位數；故意設錯 Jacobian (jacExp = 0) 的鏈應偏離。
suppressMessages(library(nimble)); nimbleOptions(verbose = FALSE, MCMCprogressBar = FALSE)
if (file.exists("../nimble/px_sampler.R")) setwd("..")
source("nimble/px_sampler.R")

n <- 4; p <- 3
code <- nimbleCode({
  Sigma[1:2, 1:2] ~ dinvwish(S = S0[1:2, 1:2], df = 4)
  for (j in 1:p) { a[j] ~ T(dnorm(1, sd = 1), 0, ); b[j] ~ dnorm(0, sd = 2) }
  for (i in 1:n) pers[i, 1:2] ~ dmnorm(mu0[1:2], cov = Sigma[1:2, 1:2])
})
m <- nimbleModel(code, constants = list(n = n, p = p, S0 = diag(2), mu0 = c(0, 0)),
                 inits = list(Sigma = diag(2), a = rep(1, p), b = rep(0, p),
                              pers = matrix(0, n, 2)))
cm <- compileNimble(m)

make_mcmc <- function(jac = NULL) {
  conf <- configureMCMC(m, nodes = NULL, monitors = c("Sigma", "a", "b", "pers"))
  for (j in 1:p) { conf$addSampler(paste0("a[", j, "]"), "RW"); conf$addSampler(paste0("b[", j, "]"), "RW") }
  for (i in 1:n) conf$addSampler(paste0("pers[", i, ", 1:2]"), "RW_block")
  conf$addSampler("Sigma[1:2, 1:2]", "RW_wishart")
  ctl <- list(scaleUp = c("pers[1:4, 1]", "b[1:3]"), scaleDown = "a[1:3]",
              covNode = "Sigma[1:2, 1:2]", covIndex = 1, scale = 0.5)
  if (!is.null(jac)) ctl$jacExpOverride <- jac
  conf$addSampler("Sigma[1:2, 1:2]", "pxScale", control = ctl)
  buildMCMC(conf)
}
mc_ok  <- make_mcmc()
mc_bad <- make_mcmc(jac = 0)
comp <- compileNimble(mc_ok, mc_bad, project = m)

run <- function(cmc) {
  cmc$run(400000, thin = 20, reset = TRUE); s <- as.matrix(cmc$mvSamples)
  s[-(1:2000), c("Sigma[1, 1]", "Sigma[1, 2]", "a[1]", "b[1]", "pers[1, 1]")]
}
set.seed(1); s_ok <- run(comp[[1]])
set.seed(2); s_bad <- run(comp[[2]])

# 精確先驗抽樣
set.seed(3); R <- 200000
Sig <- replicate(R, solve(rWishart(1, 4, diag(2))[, , 1]))
a1 <- qnorm(runif(R, pnorm(0, 1, 1), 1), 1, 1)
exact <- cbind(Sig[1, 1, ], Sig[1, 2, ], a1, rnorm(R, 0, 2), rnorm(R) * sqrt(Sig[1, 1, ]))
colnames(exact) <- colnames(s_ok)

qs <- c(.1, .25, .5, .75, .9)
tab <- function(s) sapply(colnames(s), function(k) quantile(s[, k], qs))
cat("== exact prior quantiles\n");            print(round(tab(exact), 3))
cat("== MCMC with pxScale (correct Jacobian)\n"); print(round(tab(s_ok), 3))
cat("== MCMC with pxScale (jacExp = 0, wrong)\n"); print(round(tab(s_bad), 3))
