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
