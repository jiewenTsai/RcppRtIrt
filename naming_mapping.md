# Naming Mapping (Old -> New)

This note standardizes variable and function naming for the joint RT-IRT project.
It keeps backward compatibility while providing descriptive names for research use.

## Parameter Names

| Old name | New name | Meaning |
|---|---|---|
| `theta` | `person_ability` | Person latent ability |
| `tau` | `person_speed` | Person latent speed |
| `a` | `item_discrimination` | IRT item discrimination |
| `b` | `item_difficulty` | IRT item difficulty |
| `xi` | `item_time_intensity` | RT item time intensity (intercept) |
| `sigma_t` | `rt_resid_sd` | RT residual SD by item |
| `slope_rt` | `rt_speed_loading` | Item loading of person speed in RT model |
| `beta` | `person_regression` | Regression coefficients for person layer |
| `Sigma_p` | `person_cov` | Person-layer covariance matrix |
| `loglik` | `log_likelihood` | Log-likelihood by saved draw |

## RT Model Equation Naming

Current RT mean structure:

`log_rt_ij = item_time_intensity_j - rt_speed_loading_j * person_speed_i + error_ij`

Where:
- `item_time_intensity_j` corresponds to old `xi_j`
- `rt_speed_loading_j` corresponds to old `slope_rt_j`
- `person_speed_i` corresponds to old `tau_i`

## C++ Function Alias Mapping

| New function | Backward target |
|---|---|
| `draw_item_time_intensity()` | `draw_xi()` |
| `draw_item_time_resid_sd()` | `draw_sigma_t()` |
| `draw_person_speed()` | `draw_tau()` |
| `loglik_response_time()` | `loglik_rt()` |
| `gibbs_joint_rt_irt_null()` | `gibbs_rtirt_null()` |

## R Helper Alias Mapping

| New function | Backward target |
|---|---|
| `sim_joint_rt_irt_null()` | `sim_rtirt_null()` |
| `sim_joint_rt_irt()` | `sim_rtirt()` |
| `sample_joint_rt_irt()` | `sample_rtirt()` |
| `to_coda_joint_rt_irt()` | `to_coda()` |
| `diagnose_joint_rt_irt()` | `diagnose()` |
| `posterior_means_joint_rt_irt()` | `posterior_means()` |
| `compare_ability_mean_vs_quantile()` | `compare_theta_mean_vs_quantile()` |
| `run_q_grid_ability_compare()` | `run_q_grid_theta_compare()` |
| `plot_q_grid_ability()` | `plot_q_grid_theta()` |

## Recommended Naming Style Going Forward

- Prefer `snake_case` with explicit role:
  - `person_*` for person-level latent/state parameters
  - `item_*` for item-level parameters
  - `rt_*` for response-time-specific terms
- Keep old names only for compatibility with legacy scripts.
- In new analysis scripts and manuscripts, use only the new names.

## Recommended Parallel Strategy (Maintenance Mode)

For stable and easy maintenance (R2jags-style workflow):

- Prefer chain-level parallelism only.
- In RStudio, use PSOCK backend:
  - `parallel_backend = "psock"`
- In terminal/Linux scripts, `parallel_backend = "auto"` is usually fine.
- Start from:
  - `n_chains = 2`
  - `n_cores = 2`
  - then scale up after confirming stability.

Example:

```r
fit <- sample_joint_rt_irt(
  model = "quantile",
  Y = sim$Y, logT = sim$logT, X_cov = sim$X,
  n_chains = 2, n_cores = 2,
  parallel_backend = "psock",
  n_iter = 2000, n_burnin = 1000,
  q_rt = 0.75,
  estimate_rt_slope = TRUE,
  hierarchical_rt_slope = TRUE,
  adaptive_slope_mh = TRUE
)
```

