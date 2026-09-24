# Gibbs + AGHQ score-test prototypes

Run from the repository root (`Rscript proto/<file>.R`); they source `aghq_score.R`.

| Script | What it checks | Result (2026-09-24) |
|---|---|---|
| `check_aghq_bruteforce.R` | `marg_ll()` (analytic tau + 1-D AGHQ) vs a 2-D brute-force grid | K = 9: ~1e-8, K = 15: ~1e-11 per person; Laplace (K = 1): ~7e-3 |
| `check_vs_aghq_pkg.R` | `marg_ll()` vs `aghq::aghq()` on the full 2-D integrand | agree to ~1e-10 |
| `sim_drift.R` | DGP with Sigma_12 depending on a covariate, plus `run_one()` | helper for the simulation |
| `sim_size_power.R` | DM and median-split LM rejection rates, n = 1000, p = 15 | H0: DM .02 / LM .06; linear drift: .98 / .94; step at z = .7: 1.00 / .96 |
| `nimble_aghq_scores.R` | nimbleQuad `buildAGHQ` casewise loglik and AD gradient | loglik matches (1e-5); AD gradient for 2-D random effects with nQuad = 9 does **not** match finite differences of nimble's own loglik (off by up to 14.8 on Sigma_12) |
