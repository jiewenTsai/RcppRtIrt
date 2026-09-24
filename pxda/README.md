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

## Where the discrimination bottleneck is (2026-09-24)

`a_j sqrt(S11)` is the slowest quantity in every DA variant (median ESS ~450 of 5000).

1. **Collapsing theta does not help.** `collapse_items = TRUE` draws each (a_j, d_j) with the persons
   integrated out (leave-one-item-out Gaussian posterior, slice on a_j, d_j exact), then the
   persons. `check_collapse_posterior.R` confirms it targets the same posterior (n = 60, p = 5,
   55k draws each: all quantiles agree). ESS for a is unchanged (456 → 477) at twice the R cost.
2. **The DA latent z is the bottleneck.** `diag_z_single_item.R`, theta fixed at the truth, one
   item updated by Albert–Chib: ESS for a = 1161 / 520 / 210 when a = 0.8 / 1.2 / 2.0
   (lag-1 autocorrelation .65 / .81 / .92).
3. **Liu–Wu PX-DA on the latent-response scale fixes it.** Move (z_j, a_j, d_j) → g (z_j, a_j, d_j)
   with (a_j, d_j) integrated out (`zscale_move()`; Gamma proposal plus MH weight exp(B(g − 1))).
   Single item (`diag_zscale_px_single_item.R`): ESS 818 → 1425, 437 → 1027, 176 → 662 for
   a = 0.8, 1.2, 2.0; posterior means and SDs unchanged.

Full RT-IRT, 4 datasets (`bench_collapse.R`, n = 500, p = 15, 5000 saved draws):

| variant | sec | ESS kappa | ESS rho* | med ESS a_id | min ESS a_id | ESS/s a_id (med) |
|---|---|---|---|---|---|---|
| standard expanded DA | 24.1 | 2981 | 2891 | 456 | 200 | 19.5 |
| + collapsed theta for items | 43.2 | 2927 | 2851 | 477 | 266 | 11.0 |
| + z-scale PX-DA per item | 27.1 | 2798 | 2619 | 712 | 495 | 26.3 |

The z-scale move raises median item ESS by ~56% and the worst item by ~2.5x at ~12% more time.
The theta-scale and shear moves (earlier table) act on unidentified directions only; the working
parameter that matters here is the latent-response scale of each item.
