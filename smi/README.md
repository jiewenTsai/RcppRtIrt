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
