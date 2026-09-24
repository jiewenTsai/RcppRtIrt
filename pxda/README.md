# Probit DA Gibbs with PX group moves vs nimble

`probit_da_gibbs.R` — Albert–Chib DA Gibbs for RT-IRT with cross-loadings (R, vectorised);
`sim_data.R` — data generator; `bench_da.R <seed> <outdir>` — one dataset, seven samplers.

Setting (2026-09-24): n = 500, p = 15, probit DGP, Sigma_12 = .3, Sigma_22 = .4,
lambda_j ~ N(.15, .1); 6000 iterations, 1000 burn-in (5000 saved); 4 datasets.
Identified quantities: kappa = (S12 - lbar S11)/sqrt(S11), v* = var(tau - lbar theta),
rho* = kappa/sqrt(v*), a_j sqrt(S11), (lambda_j - lbar) sqrt(S11).

Mean over the 4 datasets:

| method | sec | ESS kappa | ESS rho* | med ESS a_id | med ESS (lam - lbar) | ESS raw S11 | ESS raw lbar | ESS/s rho* |
|---|---|---|---|---|---|---|---|---|
| M1 identified DA (S11 = 1) | 24.6 | 2346 | 2302 | 347 | 2588 | 0 | 5 | 93 |
| M2 expanded DA | 24.3 | 2981 | 2891 | 456 | 3093 | 36 | 7 | 119 |
| M3 = M2 + scale move | 25.6 | 2794 | 2618 | 443 | 3034 | 936 | 8 | 102 |
| M4 = M3 + shear move | 26.7 | 3191 | 2960 | 447 | 3090 | 1027 | 1113 | 111 |
| N1 nimble default (logit, RW) | 89.4 | 1109 | 1008 | 519 | 1299 | 14 | 3 | 11 |
| N2 nimble + PG on (a_j, d_j) | 104.0 | 1558 | 1334 | 625 | 1634 | 8 | 4 | 13 |
| N3 = N2 + pxScale + pxShear | 116.5 | 1122 | 1022 | 502 | 1322 | 934 | 1190 | 9 |

- Expanded vs identified DA (M2 vs M1): ~25–30% more ESS for every identified quantity at the
  same cost (marginal augmentation).
- The group moves (M3, M4) fix mixing in the unidentified raw directions but do not improve the
  identified quantities; per-dataset differences go both ways.
- In nimble the moves reduce identified ESS (N3 < N2), plausibly because they keep changing the
  raw scale that the adaptive random-walk proposals are tuned on.
- nimble's PG sampler on the item blocks raises ESS for a_id by ~20% over the default, but costs
  ~16% more time. nimble cannot use PG on persons (their children include log T).
- Reflection: in 2 of 4 datasets all nimble chains sat in the mirror mode (theta, a, lambda ->
  negative), giving kappa of the wrong sign with the right magnitude. a ~ N(1, 1) cannot be
  truncated because nimble's PG sampler requires a dnorm prior; relabel post hoc.
- Caveats: one condition, 4 datasets; nimble uses a logit link on probit data; M1 has a
  different prior on the identified quantities than M2–M4.
