# nimble PX-DA scale sampler

- `px_sampler.R` — `sampler_pxScale`: generalized-Gibbs scale move (Liu & Sabatti 2000).
- `test_px_prior.R` — data-free model (posterior = prior). Correct Jacobian reproduces the exact
  prior quantiles (Sigma_11 median .425 vs .426); jacExp forced to 0 does not (.162).
- `bench_px.R` — expanded (unidentified) RT-IRT model, same default samplers with and without
  `pxScale`. One dataset, n = 300, p = 15, 6000 iterations (1000 burn-in), 2026-09-24:

| | default | + pxScale |
|---|---|---|
| seconds | 51.4 | 51.9 |
| ESS Sigma_11 (raw, unidentified) | 7 | 950 |
| ESS Sigma_12 / sqrt(Sigma_11) | 1059 | 567 |
| ESS rho | 330 | 205 |
| median ESS of a_j * sqrt(Sigma_11) | 436 | 402 |

The move mixes the unidentified scale direction as intended, but on this single run it does not
improve the identified quantities. Replications and an identified-parameterization baseline are
still needed before drawing conclusions.

## pxShear (cross-loading direction)

- `test_px_shear.R` — (A) with data and cross-loadings, a manual shear changes the data
  log-likelihood by 7e-14 and the person prior by 0; (B) data-free model reproduces the exact
  prior quantiles (Sigma_22 median .428 vs .424); adding 2c to log r breaks it (.116).
- `bench_px_shear.R` — n = 300, p = 10, cross-loadings, one dataset, 6000 iterations
  (1000 burn-in), 2026-09-24. kappa = (Sigma_12 - lbar Sigma_11) / sqrt(Sigma_11), truth .179.

| | default | + pxShear | + pxShear + pxScale |
|---|---|---|---|
| seconds | 58.6 | 61.1 | 63.7 |
| ESS lbar (raw, unidentified) | 1.7 | 1068 | 1067 |
| ESS Sigma_12 (raw, unidentified) | 3.5 | 1074 | 1065 |
| ESS Sigma_11 (raw, unidentified) | 9.8 | 10.8 | 693 |
| ESS kappa (identified) | 461 | 612 | 437 |
| median ESS (lambda_j - lbar) sqrt(Sigma_11) | 581 | 568 | 544 |
| median ESS a_j sqrt(Sigma_11) | 309 | 261 | 289 |
| posterior mean kappa | .173 | .167 | .170 |

Same pattern as pxScale: each move fixes mixing in its own unidentified direction at ~4% extra
time, while ESS for identified quantities stays within single-run noise. With random-walk
samplers on theta, a and lambda there is no data-augmentation coupling for the moves to break.
