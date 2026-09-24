suppressMessages({library(nimble); library(nimbleQuad)})
nimbleOptions(verbose = FALSE)
if (file.exists("../aghq_score.R")) setwd(".."); source("aghq_score.R")
set.seed(5); n <- 100; p <- 6
a <- runif(p, .7, 1.5); b <- rnorm(p, 0, .5); xi <- rnorm(p, 4, .2); sig <- runif(p, .3, .6)
S12 <- .45; S22 <- .4
th <- rnorm(n); ta <- S12 * th + rnorm(n, 0, sqrt(S22 - S12^2))
Y <- matrix(rbinom(n * p, 1, plogis(outer(th, a) - rep(1, n) %o% (a * b))), n)
logT <- matrix(xi, n, p, byrow = TRUE) - ta + matrix(rnorm(n * p), n) %*% diag(sig)

code <- nimbleCode({
  for (j in 1:p) {
    a[j] ~ T(dnorm(1, sd = 1), 0, )
    b[j] ~ dnorm(0, sd = 2)
    xi[j] ~ dnorm(4, sd = 10)
    sig[j] ~ dinvgamma(1, 1)
  }
  s22 ~ dinvgamma(1, 1)
  s12 ~ dunif(-1, 1)
  for (i in 1:n) {
    th[i] ~ dnorm(0, sd = 1)
    tau[i] ~ dnorm(s12 * th[i], var = s22 - s12^2)
    for (j in 1:p) {
      Y[i, j] ~ dbern(expit(a[j] * (th[i] - b[j])))
      logT[i, j] ~ dnorm(xi[j] - tau[i], sd = sig[j])
    }
  }
})
t0 <- proc.time()[3]
m <- nimbleModel(code, constants = list(n = n, p = p), data = list(Y = Y, logT = logT),
                 inits = list(a = a, b = b, xi = xi, sig = sig, s12 = S12, s22 = S22,
                              th = rep(0, n), tau = rep(0, n)), buildDerivs = TRUE)
cm <- compileNimble(m)
pn <- c("a", "b", "xi", "sig", "s12", "s22")
# casewise: wrapper nimbleFunction over the per-person AGHQ sets
casewise_nf <- nimbleFunction(
  setup = function(agh) { nfl <- agh$AGHQuad_nfl; ns <- length(nfl) },
  run = function(pv = double(1)) {
    ans <- matrix(0, ns, length(pv) + 1)
    for (i in 1:ns) {
      ans[i, 1] <- nfl[[i]]$calcLogLik2(pv)
      ans[i, 2:(length(pv) + 1)] <- nfl[[i]]$gr_logLik2(pv)
    }
    returnType(double(2)); return(ans)
  },
  methods = list(adDummy = function(x = double(1)) { returnType(double()); return(sum(x)) }),
  buildDerivs = "adDummy")
agh <- buildAGHQ(m, nQuad = 9, paramNodes = pn, randomEffectsNodes = c("th", "tau"))
cw <- casewise_nf(agh)
comp <- compileNimble(agh, cw, project = m)
cagh <- comp[[1]]; ccw <- comp[[2]]
cat(sprintf("model + AGHQ build/compile: %.1f s\n", proc.time()[3] - t0))

pvec <- c(a, b, xi, sig, S12, S22)
cat("param order:", head(cagh$getNodeNamesVec(TRUE), 3), "...", tail(cagh$getNodeNamesVec(TRUE), 2), "\n")
ll_nim <- cagh$calcLogLik(pvec)
psi <- c(a, b, xi, log(sig), S12, log(S22))
ll_mine <- sum(marg_ll(psi, Y, logT, gh_nodes(9)))
cat(sprintf("total marginal loglik: nimble %.8f  mine %.8f  diff %.2e\n", ll_nim, ll_mine, ll_nim - ll_mine))

t4 <- proc.time()[3]
t5 <- proc.time()[3]; out <- ccw$run(pvec); cat(sprintf("casewise run: %.3f s\n", proc.time()[3] - t5))
ll_i <- out[, 1]; G <- out[, -1]
mine_i <- marg_ll(psi, Y, logT, gh_nodes(9))
cat("max |casewise loglik nimble - mine|:", signif(max(abs(ll_i - mine_i)), 3), "\n")
S_mine <- casewise_scores(psi, Y, logT, gh_nodes(9))
jac <- c(rep(1, 3 * p), 1 / sig, 1, 1 / S22)   # d log(x)/dx for sig, s22
cat("max |AD casewise score - my numeric score|:", signif(max(abs(G - sweep(S_mine, 2, jac, "*"))), 3), "\n")


# --- diagnose: nimble AD casewise score vs finite differences of nimble's own casewise loglik
h <- 1e-5
FD <- sapply(seq_along(pvec), function(k) {
  e <- replace(numeric(length(pvec)), k, h)
  (ccw$run(pvec + e)[, 1] - ccw$run(pvec - e)[, 1]) / (2 * h)
})
nm <- cagh$getNodeNamesVec(TRUE)
blk <- sub("\\[.*", "", nm)
mine_nat <- sweep(S_mine, 2, jac, "*")
cat("\nmax |diff| by parameter block:\n")
print(round(rbind(
  AD_vs_nimbleFD = tapply(apply(abs(G - FD), 2, max), blk, max),
  mine_vs_nimbleFD = tapply(apply(abs(mine_nat - FD), 2, max), blk, max)), 5))
cat("sum of AD casewise vs top-level gr_logLik:", signif(max(abs(colSums(G) - cagh$gr_logLik(pvec))), 3), "\n")
