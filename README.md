# Joint RT-IRT Modeling (Rcpp)

Rcpp implementation of joint response-time and item-response modeling with:

- PG-based logistic IRT sampling
- Mean RT model
- Quantile RT model (BQR-style augmentation)
- Cross-quantile RT model
- Optional RT slope and hierarchical RT slope

This project is designed for reproducible comparison of `theta` recovery under different RT assumptions, including Quantile Invariance (QI) checks.

---

## Why This Repo

Traditional joint RT-IRT models often assume no direct ability effect on RT (or fixed effect across conditions).  
This codebase provides a practical workflow to test when that assumption is stable vs. quantile-dependent.

Key goals:

- Compare mean vs quantile vs cross-quantile RT formulations
- Evaluate effect on `theta` estimation quality (RMSE, bias, correlation)
- Support simulation and real-data analysis in one consistent pipeline

---

## Current Model Coverage

- `structural` (mean RT + latent regression)
- `quantile` (partial quantile RT)
- `cross_quantile` (quantile RT with cross loading)
- `ml_irt`, `null` variants

Notes:

- `cross_quantile` is implemented as a practical first version.
- Bifactor extension is not included yet.

---

## Installation

In R:

```r
install.packages(c("Rcpp", "RcppArmadillo", "MASS", "coda"))
# pg package is required (GitHub package tmsalab/pg)
```

Then compile/load:

```r
library(Rcpp)
sourceCpp("RtIrtGibbs.cpp")
source("RtIrtGibbs_helpers.R")
.rtirt_cpp_path <- normalizePath("RtIrtGibbs.cpp")
```

---

## Quick Start

### 1) Simulate data

```r
sim <- sim_joint_rt_irt(
  n_subj = 300, n_item = 15, n_feat = 2,
  seed = 42, rt_error = "t", df_t = 4
)
```

### 2) Fit mean vs quantile and compare `theta`

```r
out <- compare_ability_mean_vs_quantile(
  sim,
  n_iter = 2000,
  n_burnin = 1000,
  n_chains = 2,
  n_cores = 2,
  parallel_backend = "psock",  # recommended in RStudio
  q_rt = 0.75,
  quantile_estimate_rt_slope = TRUE,
  quantile_hierarchical_rt_slope = TRUE,
  quantile_adaptive_slope_mh = TRUE
)

out$theta_compare
```

### 3) Run q-grid benchmark

```r
simgen <- function(i) sim_joint_rt_irt(
  n_subj = 200, n_item = 12, n_feat = 2,
  seed = 100 + i, rt_error = "t", df_t = 4
)

out_grid <- run_q_grid_ability_compare(
  sim_generator = simgen,
  q_grid = c(0.25, 0.5, 0.75, 0.85),
  n_rep = 10,
  n_iter = 2000,
  n_burnin = 1000,
  n_chains = 2,
  n_cores = 2,
  parallel_backend = "psock",
  quantile_estimate_rt_slope = TRUE,
  quantile_hierarchical_rt_slope = TRUE,
  quantile_adaptive_slope_mh = TRUE
)

out_grid$summary
plot_q_grid_ability(out_grid)
```

---

## Parallel Strategy (Recommended)

For maintainability/stability (R2jags-style chain parallelism):

- In RStudio: use `parallel_backend = "psock"`
- In terminal/Linux: `parallel_backend = "auto"` is usually fine
- Start with `n_chains = 2`, `n_cores = 2`

---

## Naming Convention

Descriptive output aliases are provided in addition to legacy names:

- `person_ability` (`theta`)
- `person_speed` (`tau`)
- `item_time_intensity` (`xi`)
- `rt_resid_sd` (`sigma_t`)
- `rt_speed_loading` (`slope_rt`)
- `person_cov` (`Sigma_p`)
- `log_likelihood` (`loglik`)

See `naming_mapping.md` for full mapping.

---

## Known Limitations

- Cross-quantile model is a practical first implementation.
- Bifactor cross-quantile is not yet implemented.
- Extreme quantiles may need stronger identification/regularization depending on data condition.

---

## License / Citation

If you use this code in research, please cite your dissertation and related model papers (PG logistic augmentation, Bayesian quantile regression, and joint RT-IRT references).
