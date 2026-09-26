# RTMB implementation of SMI stage 1

`smi_rtmb.R` fits the marginally tempered stage-1 target of `smi/smi_gibbs.R` with RTMB: the RT
marginal likelihood (tau integrated analytically, Omega = diag(s^2) + v 11') is raised to eta
directly, theta is integrated by the Laplace approximation, fixed effects are at the posterior mode
(same priors as the Gibbs sampler), and targets (group gaps, slopes) are linear functionals of theta
reported with `sdreport()`. eta = 0 maps the RT parameters out (the cut).

Requires TMB >= 1.9.25 and RTMB 2.0 built from GitHub (Ubuntu's r-cran-tmb 1.9.10 is too old); the
session-start hook installs both.

`compare_gibbs_rtmb.R` — 4 datasets (`sim_gz`, Z slower by 0 or 0.5), eta in {0, .1, .25, .5, .75, 1},
targets: gap by G (in the model) and by Z (not in the model).

- Estimates agree with the Gibbs posterior means to within 0.011 (max over 48 comparisons); RTMB
  SEs / Gibbs posterior SDs: median ratio 1.02.
- Risk-rule eta for the Z gap under the speed shift: Gibbs 0.1 / RTMB 0.0 (both datasets); without
  the shift the two agree (0.25 and 1). For the G gap the risk curve is flat across eta, so the
  chosen eta (0.1-0.5) varies between methods without changing the estimate.
- Time per fit: RTMB 2.1 s, Gibbs 10.1 s (3000 iterations).
