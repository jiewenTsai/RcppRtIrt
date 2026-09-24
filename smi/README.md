# Semi-modular inference (SMI) for probit RT-IRT

`smi_gibbs.R` — stage 1 of SMI (RT likelihood raised to eta; tau prior untempered) plus optional
nested stage 2 (`K_inner` sweeps of the RT module given theta). Identified parameterisation:
tau ⟂ theta, the ability-to-RT relation is carried by gamma_j. eta = 0 is the cut posterior for
(theta, a, d); eta = 1 is the full posterior. Uses `rtnorm_side()` and `zscale_move()` from
`pxda/probit_da_gibbs.R`; run scripts from the repository root.

## Experiments (2026-09-24; n = 500, 3 datasets per cell, 3000 iterations, 500 burn-in)

`exp_eta.R` (p = 15, gamma ~ .3, v = .3): RT carries little information about theta (posterior
SD .361 at eta = 0 vs .341 at eta = 1), and the "misspecified" scenario (gamma .4 vs -.2) pooled to
gamma ≈ .1, so nothing moved. Kept as a negative control.

`exp_speed_shift.R` and `exp_temper.R` (p = 10, gamma ~ .5, v = .15): group B equally able but
slower by `shift` on the log-time scale; the model assumes one speed distribution. `gap_bias` is
the error in the estimated B - A ability gap.

| shift | eta | gap_bias | RMSE theta | mean post SD |
|---|---|---|---|---|
| 0 | 0 | .048 | .447 | .439 |
| 0 | .05 | .046 | .431 | .429 |
| 0 | .10 | .041 | .414 | .418 |
| 0 | .25 | .032 | .386 | .391 |
| 0 | 1 | .027 | .382 | .376 |
| .5 | 0 | .048 | .447 | .439 |
| .5 | .05 | .076 | .434 | .430 |
| .5 | .10 | .110 | .419 | .420 |
| .5 | .25 | .188 | .401 | .397 |
| .5 | 1 | .190 | .398 | .386 |

(The .048 at eta = 0 is sampling noise shared by all rows: same 3 datasets.)

- Using RT lowers theta RMSE by ~15% when the RT module is right.
- A speed shift in one group distorts the estimated group gap by ~0.14 SD at eta ≥ .25, yet
  overall RMSE still favours eta = 1 (.398 vs .447): an RMSE criterion would pick the model that
  biases group comparisons.
- The effect of eta saturates by ~.25. With the tau prior untempered, the information about
  theta that passes through tau is capped by the speed variance v, so eta is a nonlinear dial.
- Tempering the whole module including the tau prior is improper for eta < 1 (see the header of
  `smi_gibbs.R`); a graded version has to temper the marginal RT likelihood.

## Marginal tempering (`temper = "marginal"`)

Raises the marginal RT likelihood N(T_i; xi + gamma theta_i, Omega), Omega = diag(sigma^2) + v 11',
to the power eta, via tau~ ~ N(0, v/eta), T | tau~ ~ N(., diag(sigma^2)/eta) plus the correction
|Omega|^(n(1-eta)/2); sigma^2 and v by log-scale slice sampling.

- `check_marginal_identity.R`: the tau~-integral differs from eta * log L_i by a constant
  (spread 3.6e-15 over random parameter draws); determinant lemma exact.
- `check_marginal_sampler.R`: n = 80, p = 4, eta = .3, 55k draws each, against an independent
  sampler that never introduces tau~ (Gaussian theta conditional under the tempered marginal,
  random-walk MH on xi, gamma, log sigma^2, log v). Quantiles of theta_1, theta_2, a_1 agree to
  <= .03; v and gamma_1 medians differ by .005 and .015 (≈ .07 posterior SD).
- `exp_temper_compare.R`: same design as `exp_temper.R`. info_frac = share of the cut-to-full gain
  in mean posterior precision of theta.

| eta | info_frac likelihood | info_frac marginal | gap_bias likelihood (shift .5) | gap_bias marginal (shift .5) |
|---|---|---|---|---|
| 0 | 0 | 0 | .048 | .048 |
| .05 | .13 | .06 | .076 | .052 |
| .10 | .29 | .08 | .110 | .058 |
| .25 | .71 | .19 | .188 | .073 |
| .50 | .90 | .38 | .196 | .105 |
| .75 | .97 | .63 | .192 | .145 |
| 1 | 1 | 1 | .190 | .191 |

Marginal tempering turns eta into a graded dial (information borrowed rises roughly in step with
eta, slightly convex), whereas likelihood tempering passes ~70% of the information by eta = .25.
At matched information the two give similar bias (e.g. info ≈ .3–.4: gap bias ≈ .11 for both), so
the choice mainly changes how eta maps to borrowing, not the bias–information frontier.

## Choosing eta (`exp_choose_eta.R`)

Target: the group gap Delta = mean(theta_B) - mean(theta_A). The cut posterior (eta = 0) is the
unbiased reference: D(eta) = Delta_hat(eta) - Delta_hat(0); under a correct RT module
Var(D) ≈ V(0) - V(eta) (Hausman). Rules on the grid eta = 0, .1, .25, .4, .55, .75, 1:
- **risk**: argmin R(eta) = D^2 - V(0) + 2 V(eta) (unbiased estimate of MSE of Delta_hat(eta));
- **hausman**: largest eta before the first H(eta) = D^2 / (V(0) - V(eta)) ≥ 3.84.

Marginal tempering; n = 500, p = 10, 6 datasets per shift (seeds 301–306), 3000 iterations.

| shift | rule | mean eta | mean gap error | RMSE gap | RMSE theta |
|---|---|---|---|---|---|
| 0 | cut | 0 | .003 | .018 | .421 |
| 0 | full | 1 | .023 | .029 | .368 |
| 0 | risk | .72 | .013 | .020 | .376 |
| 0 | hausman | .90 | .019 | .029 | .371 |
| .25 | cut | 0 | .003 | .018 | .421 |
| .25 | full | 1 | .106 | .108 | .373 |
| .25 | risk | .10 | .009 | .019 | .416 |
| .25 | hausman | .15 | .015 | .025 | .412 |
| .5 | cut | 0 | .003 | .018 | .421 |
| .5 | full | 1 | .159 | .161 | .381 |
| .5 | risk | .07 | .010 | .021 | .418 |
| .5 | hausman | .08 | .012 | .020 | .417 |

Both rules keep most of the precision gain when the RT module is (nearly) right and fall back to
near-cut when one group is slower: gap RMSE .019–.025 instead of .108–.161. The small full-model gap
error at shift 0 (.023) comes from the item-level gamma differences between groups in `sim_smi()`.

## Stage 2: inference on the RT module (`rt_sweeps.cpp`, `exp_stage2.R`, `exp_stage2_eta.R`)

- `rt_sweeps_cpp()` reproduces the R sweep in distribution (theta fixed, 4500 sweeps: gamma-bar and
  v quantiles agree to ~.002) at 0.076 ms vs 0.85 ms per sweep (~11x).
- Nested MCMC with warm-started inner chains: K = 1 already matches a reference that runs a fresh
  300-sweep chain for each of 500 thinned theta draws (gamma-bar median .368 vs .372, v .232 vs .230).
  The naive-cut bias is negligible here because the RT-module Gibbs mixes fast and theta moves slowly.
- **Attenuation under the cut.** The cut draws theta from p(theta | Y), not conditioned on RT, so the
  stage-2 regression of log T on theta draws is an errors-in-variables regression:

| eta | gamma-bar hat / truth | v hat (truth .15) | reliability of theta (1 - mean post var) |
|---|---|---|---|
| 0 | .80 | .235 | .817 |
| .25 | .84 | .222 | .829 |
| .5 | .89 | .206 | .839 |
| .75 | .93 | .188 | .850 |
| 1 | .99 | .167 | .862 |

  (3 datasets, no speed shift.) The cut shrinks gamma by the reliability of theta, as plausible
  values drawn without the secondary variable in the conditioning model would. SMI therefore trades
  two errors: high eta biases group comparisons of theta when one group is slower; low eta
  attenuates the speed–ability relation.

## Cut vs hierarchical (conditioning) model — `exp_hier_vs_cut.R`

n = 500, p = 10, 20 datasets. G: reporting group with a true ability gap of 0.3 SD. Z: correlated
with G (P(Z=1|G) = .3/.7), no ability effect, Z = 1 slower by `shift` in tau. "+G" puts G in the
theta mean (and, for full, the tau mean); "+GZ" puts both. Bias and RMSE of the posterior-mean
theta gap by G and by Z; cov = 95% interval coverage.

| shift | method | bias G | bias Z | cov G | cov Z | RMSE G | RMSE Z | theta RMSE |
|---|---|---|---|---|---|---|---|---|
| 0 | full | -0.043 | 0.004 | .75 | .95 | 0.055 | 0.031 | 0.37 |
| 0 | full+G | -0.004 | 0.019 | .95 | .95 | 0.040 | 0.038 | 0.37 |
| 0 | cut | -0.055 | -0.001 | .70 | .90 | 0.068 | 0.038 | 0.44 |
| 0 | cut+G | 0.003 | 0.021 | .95 | .85 | 0.039 | 0.043 | 0.44 |
| 0.5 | full | 0.011 | 0.150 | .85 | .00 | 0.036 | 0.154 | 0.38 |
| 0.5 | full+G | -0.003 | 0.147 | .95 | .00 | 0.039 | 0.150 | 0.38 |
| 0.5 | full+GZ | -0.004 | 0.019 | .95 | .80 | 0.040 | 0.048 | 0.37 |
| 0.5 | cut | -0.055 | -0.001 | .70 | .90 | 0.068 | 0.038 | 0.44 |
| 0.5 | cut+G | 0.003 | 0.021 | .95 | .85 | 0.039 | 0.043 | 0.44 |

(cut rows do not depend on the shift: the cut never sees RT.)

- Without a conditioning model both full and cut shrink a real group gap (cut by about 1 − reliability).
  So the cut is an unbiased reference for a gap only when theta's prior conditions on that grouping;
  the eta rules in `exp_choose_eta.R` were run with a true gap of 0 and must be redone with `Xth`.
- For a variable in the model (G), the hierarchical full model is unbiased and matches the cut in
  gap RMSE. The RT precision gain for individual theta (0.37 vs 0.44) does not carry over to group
  means (posterior SD of the G gap 0.042 vs 0.044).
- For a speed-related variable left out of the RT model (Z), full+G is biased by 0.15 SD with 0%
  coverage; cut+G is not. Only the oracle full+GZ fixes it.
