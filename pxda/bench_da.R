# bench_da.R — one dataset, seven samplers; writes an RDS with ESS / ESS-per-second / means
# Usage: Rscript pxda/bench_da.R <seed> <outdir>
suppressMessages({library(nimble); library(coda)})
nimbleOptions(verbose = FALSE, MCMCprogressBar = FALSE)
if (file.exists("../pxda/probit_da_gibbs.R")) setwd("..")
source("pxda/probit_da_gibbs.R"); source("pxda/sim_data.R"); source("nimble/px_sampler.R")
args <- commandArgs(TRUE)
seed <- as.integer(if (length(args) >= 1) args[1] else 1)
outdir <- if (length(args) >= 2) args[2] else "."
n <- 500; p <- 15; n_iter <- 6000; n_burn <- 1000
dat <- sim_pxda(n, p, seed = seed)

summarise <- function(S, secs, label) {
  idd <- identified_draws(S)
  ess_s <- effectiveSize(mcmc(idd$summary))
  ess_raw <- effectiveSize(mcmc(idd$raw))
  c(method = label, secs = secs,
    setNames(ess_s, paste0("ESS_", names(ess_s))),
    ESS_a_med = median(effectiveSize(mcmc(idd$a_id))),
    ESS_lamdev_med = median(effectiveSize(mcmc(idd$lamdev_id))),
    setNames(ess_raw, paste0("ESS_raw_", names(ess_raw))),
    setNames(colMeans(idd$summary), paste0("mean_", colnames(idd$summary))))
}
res <- list()
rpart_file <- file.path(outdir, sprintf("bench_da_seed%d_Rpart.rds", seed))

## R probit DA Gibbs ---------------------------------------------------------
if (file.exists(rpart_file)) res <- readRDS(rpart_file) else {
cfgs <- list(M1_identified = list("identified", FALSE, FALSE),
             M2_expanded   = list("expanded", FALSE, FALSE),
             M3_scale      = list("expanded", TRUE, FALSE),
             M4_scale_shear = list("expanded", TRUE, TRUE))
for (nm in names(cfgs)) {
  cf <- cfgs[[nm]]
  t0 <- proc.time()[3]
  S <- probit_da_gibbs(dat$Y, dat$logT, n_iter, n_burn, method = cf[[1]],
                       scale_move = cf[[2]], shear_move = cf[[3]], seed = seed)
  res[[nm]] <- summarise(S, proc.time()[3] - t0, nm)
}
saveRDS(res, rpart_file)
}

## nimble (logit) -------------------------------------------------------------
code <- nimbleCode({
  Sigma[1:2, 1:2] ~ dinvwish(S = S0[1:2, 1:2], df = 4)
  for (j in 1:p) {
    a[j] ~ dnorm(1, sd = 1); d[j] ~ dnorm(0, sd = 2); xi[j] ~ dnorm(4, sd = 10)
    s2[j] ~ dinvgamma(1, 1); lam[j] ~ dnorm(0, sd = 1)
  }
  for (i in 1:n) {
    pers[i, 1:2] ~ dmnorm(mu0[1:2], cov = Sigma[1:2, 1:2])
    for (j in 1:p) {
      logit(pr[i, j]) <- a[j] * pers[i, 1] - d[j]
      Y[i, j] ~ dbern(pr[i, j])
      logT[i, j] ~ dnorm(xi[j] - pers[i, 2] + lam[j] * pers[i, 1], var = s2[j])
    }
  }
})
m <- nimbleModel(code, constants = list(n = n, p = p, S0 = diag(2), mu0 = c(0, 0)),
                 data = list(Y = dat$Y, logT = dat$logT),
                 inits = list(Sigma = diag(2), a = rep(1, p), d = rep(0, p), xi = colMeans(dat$logT),
                              s2 = rep(.25, p), lam = rep(0, p), pers = matrix(0, n, 2)))
cm <- compileNimble(m)
mon <- c("Sigma", "a", "lam")
add_pg <- function(conf) {
  for (j in 1:p) {
    tg <- c(sprintf("a[%d]", j), sprintf("d[%d]", j))
    conf$removeSamplers(tg)
    conf$addSampler(target = tg, type = "polyagamma",
                    control = list(fixedDesignColumns = c(FALSE, TRUE),
                                   nonTargetNodes = sprintf("pers[1:%d, 1:2]", n)))
  }
  conf
}
add_px <- function(conf) {
  conf$addSampler("Sigma[1:2, 1:2]", "pxScale",
    control = list(scaleUp = sprintf("pers[1:%d, 1]", n),
                   scaleDown = c(sprintf("a[1:%d]", p), sprintf("lam[1:%d]", p)),
                   covNode = "Sigma[1:2, 1:2]", skipInvariant = TRUE, scale = 0.05))
  conf$addSampler("Sigma[1:2, 1:2]", "pxShear",
    control = list(sourceNodes = sprintf("pers[1:%d, 1]", n), targetNodes = sprintf("pers[1:%d, 2]", n),
                   shiftNodes = sprintf("lam[1:%d]", p), covNode = "Sigma[1:2, 1:2]",
                   skipInvariant = TRUE, scale = 0.05))
  conf
}
confs <- list(N1_nimble_default = configureMCMC(m, monitors = mon),
              N2_nimble_PG      = add_pg(configureMCMC(m, monitors = mon)),
              N3_nimble_PG_px   = add_px(add_pg(configureMCMC(m, monitors = mon))))
cmcs <- do.call(compileNimble, c(lapply(confs, buildMCMC), list(project = m)))
if (!is.list(cmcs)) cmcs <- list(cmcs)
for (k in seq_along(confs)) {
  set.seed(seed)
  t0 <- proc.time()[3]; cmcs[[k]]$run(n_iter); secs <- proc.time()[3] - t0
  s <- as.matrix(cmcs[[k]]$mvSamples)[-(1:n_burn), ]
  S <- cbind(S11 = s[, "Sigma[1, 1]"], S12 = s[, "Sigma[1, 2]"], S22 = s[, "Sigma[2, 2]"],
             lbar = rowMeans(s[, sprintf("lam[%d]", 1:p)]),
             s[, sprintf("a[%d]", 1:p)], s[, sprintf("lam[%d]", 1:p)])
  colnames(S) <- c("S11", "S12", "S22", "lbar", paste0("a", 1:p), paste0("lam", 1:p))
  res[[names(confs)[k]]] <- summarise(S, secs, names(confs)[k])
}
out <- list(seed = seed, truth = dat$truth, res = do.call(rbind, res))
saveRDS(out, file.path(outdir, sprintf("bench_da_seed%d.rds", seed)))
print(out$res[, c("method", "secs", "ESS_kappa", "ESS_rho_star", "ESS_a_med")])
