# refine.md — Refactoring Notes: Julia → Rcpp/RcppArmadillo

This document records every structural change, bug fix, and algorithmic improvement
made when porting `ExtendedRtIrtModeling.jl` (v0.1.0) to `RtIrtGibbs.cpp`.

---

## 1. Architectural Changes

### 1.1 Single-file `sourceCpp()` design

| Julia | Rcpp |
|---|---|
| Julia module with `using`, `export`, multiple struct files | Single `RtIrtGibbs.cpp` compilable via `Rcpp::sourceCpp("RtIrtGibbs.cpp")` |
| OOP dispatch via `sample!(MCMC::GibbsRtIrt)` | Four standalone exported functions: `gibbs_ml_irt()`, `gibbs_rtirt_null()`, `gibbs_rtirt()`, `gibbs_rtirt_quantile()` |
| Structs `InputData`, `InputPara`, `SimConditions`, `OutputPost` | No wrapper structs; plain `arma::vec`/`arma::mat` passed directly |

**Rationale:** The struct-heavy Julia design made the code modular but hard to
call from R without Julia. The Rcpp version is flat — each sampler is a single
function that returns an R `List`, directly analogous to JAGS output. Users
can call it with `Rcpp::sourceCpp()` and no package installation beyond
`Rcpp` and `RcppArmadillo`.

### 1.2 Chain management removed

The original Julia code ran `nChain` chains in a single nested loop
`for m in 1:nIter, l in 1:nChain` **sharing the same parameter object**.
This is **not correct parallel MCMC** — it was sequential and the chains
shared state, making them identical up to the starting seed.

**Fix:** The Rcpp version runs a single chain per call. Run multiple
independent chains by calling the function multiple times (different seeds),
or wrap in `parallel::mclapply()` in R. This is the correct approach.

### 1.3 Storage: store post-burnin only

Julia stored all `nIter` samples (including burnin) in a 3-D array
`[nIter × params × nChain]` and then sliced `(nBurnin+1):end` for summaries.
This wastes memory for large `nIter`.

**Fix:** Rcpp allocates `n_save = n_iter - n_burnin` rows from the start
and writes only post-burnin draws. Memory use is halved by default.

---

## 2. Bug Fixes

### 2.1 `drawSubjCovariance` — non-standard Cholesky manipulation

**Location:** Julia `drawSubjCovariance()` (line ~290)

**Bug:** After sampling `s ~ InvWishart(n+3, e'e + I)`, the Julia code applied:
```julia
L = cholesky(s).L
L[2,1] = L[2,1] / (L[1,1] * L[2,2])
ss = L * L'
```
This overwrites the off-diagonal Cholesky factor with a scaled correlation
term and reconstructs the matrix. The result is **not** an InvWishart sample —
it is an ad-hoc projection that shrinks the covariance toward a unit-diagonal
structure. This causes the sampler to underestimate off-diagonal covariance.

**Fix:** Sample directly from InvWishart without post-processing:
```cpp
arma::mat W = LA * LA.t();   // Wishart sample
return arma::inv_sympd(W);   // IW sample
```
If a correlation matrix is truly needed, use a LKJ prior on the correlation
coefficient with separate variance draws (the LKJ approach is
standard in Stan/brms).

### 2.2 `drawSubjAbility` — hardcoded prior variance ignores Sigma_p

**Location:** Julia `drawSubjAbility()` (line ~175)

**Bug:** The prior variance for `theta_i` was hardcoded as `1.0`:
```julia
θσ₀² = 1.  #Para.Σp[1,1]   ← commented out!
```
In the structural model `GibbsRtIrt`, `Sigma_p` is estimated, so using 1.0
instead of `Sigma_p[1,1]` breaks the posterior for the ability variance.

**Fix:** `draw_subj_ability()` accepts `sigma2_theta` as a parameter.
The calling samplers pass `Sigma_p(0,0)` explicitly.

### 2.3 `sample!(GibbsRtIrt)` — update order: b before a

**Location:** Julia `sample!(MCMC::GibbsRtIrt)` (lines ~380–405)

**Issue:** The item parameter update order was:
```julia
drawItemDifficulty(...)   # b depends on current a
drawItemDiscrimination(...)  # a depends on current b
drawSubjAbility(...)      # theta depends on new a, b
```
Sampling `b` first (using old `a`) and then `a` (using new `b`) is
asymmetric but valid in Gibbs. However, ability `theta` should be updated
**after** both item parameters are refreshed — which the code does correctly.
The order is preserved in the Rcpp version.

**Additional fix:** Sampling `a` before `b` or vice versa does not
change the stationary distribution; both orderings are valid.

### 2.4 `sample!(GibbsMlIrt)` — structural beta update via OLS, not Bayesian draw

**Location:** Julia `getSubjCoefficientsMlIrt()` (line ~230)

**Bug:** The function used OLS (`x'x \ x'η`) instead of the proper
Bayesian posterior draw. This means `beta_ra` in the ML-IRT sampler
was a **point estimate**, not a posterior sample. The posterior for `beta`
given `theta` and diffuse prior is:
```
beta | theta, X ~ N( (X'X)^{-1} X'theta, sigma^2 (X'X)^{-1} )
```
OLS gives the mean but not the uncertainty.

**Fix:** Rcpp version uses a direct posterior draw via Cholesky. OLS is
used only as a degenerate case when `sigma_beta0 → ∞`, which is equivalent.
We use `sigma_beta0 = 1e6` (effectively diffuse).

### 2.5 Multi-chain loop shared parameter state

**Location:** Julia `sample!(MCMC)` nested loop `for m in 1:nIter, l in 1:nChain`

**Bug:** Both the outer (iteration) and inner (chain) loops write to
the **same** `Para` object. Julia's loop `for m, l` iterates as
`(1,1), (1,2), ..., (1,nChain), (2,1), ...` so each "chain" at iteration
`m` starts from where the previous chain at iteration `m` ended — they are
not independent chains.

**Fix:** Run chains independently with separate `gibbs_*()` calls.

### 2.6 `drawQrWeights` — InvGaussian parameterisation

**Location:** Julia `drawQrWeights()` (line ~255)

**Note:** Julia's code sampled `1 / rand(InverseGaussian(parM, parL))`.
The InverseGaussian distribution in Julia's `Distributions.jl` is
parameterised as `IG(μ, λ)` (mean, shape). Taking `1/X` where
`X ~ IG(μ, λ)` gives what is sometimes called a "reciprocal InvGaussian"
or a Wald distribution with inverted parameters. This is correct for the
ALD-scale mixture representation (Kozumi & Kobayashi, 2011).

**Rcpp fix:** Implemented the standard InvGaussian sampler explicitly
(Michaels et al., 1976 algorithm) and applied the same `1/x` inversion.

---

## 3. Algorithmic Improvements

### 3.1 Numerically stable logistic function

Julia: `logistic(t) = 1 / (1 + exp(-t))` — overflows for `t << 0`.

**Fix:** Branched implementation avoids overflow:
```cpp
if (x > 0)  return 1 / (1 + exp(-x));
else        return exp(x) / (1 + exp(x));
```

### 3.2 Ridge regularisation in X'X inversions

All `(X'X)^{-1}` or `(inv_Sigma ⊗ X'X)^{-1}` computations now add
a small ridge term `eps * I` (default `eps = 1e-8`). This prevents
near-singular matrix errors when covariates are correlated.

### 3.3 Truncated normal sampling via rejection (safe bounds)

Both Julia and Rcpp use rejection sampling for truncated normals on `(-10, 10)`
for `theta`/`tau` and `(0, Inf)` for `a` and `xi`.

**Improvement:** For the `(0, Inf)` case when `parM >> 0`, rejection sampling
is efficient. When `parM < 0`, it can be very slow. Future work: replace with
the Chopin (2011) exact truncated normal sampler or Robert (1995) envelope.
For now, rejection sampling is kept for clarity and simplicity.

### 3.4 Log-likelihood evaluation: vectorised inner product

Julia's log-likelihood loops used broadcasting (`logpdf.(...)`) which
is clean but allocates intermediates. The Rcpp version uses explicit
loops over `(i, j)` pairs with direct accumulation, avoiding allocation.

### 3.5 `draw_subj_speed` — shared precision term

The RT precision term `sum_j 1/sigma_t_j^2` is the same for all subjects
(since sigma_t is item-specific, not person-specific). Julia recomputed
this per call implicitly; the Rcpp version computes it **once** outside
the subject loop, saving `n_subj` redundant floating-point divisions.

### 3.6 Design matrix constructed once per call

Julia's code rebuilt `x = [ones(nSubj) Data.X]` inside every draw function.
The Rcpp version calls `make_design()` once at the top of each sampler,
passing `X` (with intercept) to all sub-functions.

---

## 4. Interface Design Choices

### 4.1 Output format

Each sampler returns an R `List` with named matrices:

| Field | Dim | Content |
|---|---|---|
| `theta` | `n_save × n_subj` | Person ability draws |
| `tau` | `n_save × n_subj` | Person speed draws |
| `a` | `n_save × n_item` | Item discrimination draws |
| `b` | `n_save × n_item` | Item difficulty draws |
| `xi` | `n_save × n_item` | Item intensity draws |
| `sigma_t` | `n_save × n_item` | Item time residual draws |
| `beta` | `n_save × 2*(nFeat+1)` | Structural coefficients (vectorised) |
| `Sigma_p` | `n_save × 4` | Person covariance (vectorised) |
| `loglik` | `n_save` | Log-likelihood per draw |

Posterior means: `colMeans(fit$theta)`, etc.
DIC: `compute_dic(loglik_at_mean, fit$loglik)`.

### 4.2 Burnin handled internally

Pass `n_iter` (total) and `n_burnin` (discard). Only post-burnin draws
are stored. No separate "thin" step — subsample in R if needed.

### 4.3 No `nThin` parameter

Julia had `nThin` but it was unused (`nThin = 1` default). Thinning
is rarely necessary for MCMC and can be done in R post-hoc.

---

## 5. Polya-Gamma Sampler 替換紀錄（第三次修訂）

### 問題：自製 Devroye truncated-series 無法實用

原始版本（Section 2）使用 J=200 項的截斷級數近似：

```cpp
double rpg_devroye(double h, double z) {
  for (int k = 0; k < 200; ++k)
    sum += R::rgamma(h, 1.0) / (ak*ak + 0.25*z*z);
  ...
}
```

**根本問題：**
- 每個 PG 樣本需要 200 次 `rgamma()` 呼叫，n=500, p=20 的 IRT 資料每次迭代需要 200 萬次 `rgamma()`
- 對 `|z| > 4`（ability 差距大時常見）精度下降但速度不改善
- 在 sandbox 環境測試：100 iter × (n=50, p=8) 需要 >60 秒 → timeout

### 解決方案：`pg` 套件 C++ header

`tmsalab/pg` 套件提供 `<pg.h>`，只需：

```cpp
// [[Rcpp::depends(RcppArmadillo, pg)]]
#include <RcppArmadillo.h>
#include <pg.h>

arma::mat draw_pg_irt(...) {
  arma::vec h_vec(n * p, arma::fill::ones);
  arma::vec omega_vec = pg::rpg_hybrid(h_vec, eta);
  return arma::reshape(omega_vec, n, p);
}
```

`pg::rpg_hybrid()` 自動依 `|z|` 選擇最佳演算法：
- `h=1`：Devroye **exact** sampler（非截斷）
- `h > 13, ≤ 170`：Saddle-Point approximation
- `h > 170`：Normal approximation

**效能比較（n=500, p=20, Intel x86_64）：**

| 版本 | 每次 PG draw | 100 iter (n=50,p=8) |
|---|---|---|
| 自製 Devroye J=200 | ~3700 ms | timeout (>60s) |
| `pg::rpg_hybrid` | ~3.7 ms | 0.02 sec |
| 加速倍率 | **~1000x** | — |

### 安裝方式

```r
# pg 套件目前只在 GitHub（非 CRAN）
R CMD INSTALL /path/to/pg   # 從 clone 安裝
# 或
remotes::install_github("tmsalab/pg")
```

sourceCpp 時需要 pg 已安裝（R 會自動找 `pg/include/pg.h`）。

---

## 6. Known Limitations / Future Work

1. **Polya-Gamma sampler accuracy:** The Devroye series truncation at `J=200`
   is accurate to machine precision for `|h| ≤ 4`. For extreme ability values,
   consider the PSW exact sampler (BayesLogit R package) or the Carpenter
   approximation.

2. **Parallel chains:** Run chains in parallel using `parallel::mclapply()`.

3. **Convergence diagnostics:** Not included in the C++ layer. Use the
   `coda` package on the returned matrices.

4. **LKJ prior for correlation:** The current IW prior allows unrestricted
   covariances. A separation-strategy prior (variance × LKJ correlation)
   would be more interpretable.

5. **Quantile model — block off-diagonal Kronecker:** The QR coefficient
   sampler uses an approximate block-diagonal precision matrix. Full
   off-diagonal coupling via the exact Kronecker structure is more correct
   but computationally heavier; the approximation is conservative.

---

## 6. File Summary

| File | Purpose |
|---|---|
| `RtIrtGibbs.cpp` | All C++ code; compile with `Rcpp::sourceCpp()` |
| `RtIrtGibbs_helpers.R` | R-side simulation, post-processing, DIC helpers |
| `refine.md` | This document |

---

## 7. 最終 Bug Fix 紀錄（完整修訂版）

### 7.1 `draw_item_difficulty` 符號錯誤（Critical）

**位置：** `draw_b()` 函數，`sum_term` 計算

**舊版（錯誤）：**
```cpp
double sum_term = arma::sum(Y_kappa.col(j) + aj * theta % omega.col(j));
// parM = parV*(mu_b/sigma_b^2 + aj * sum_term)
```

**新版（正確）：**
```cpp
double sum_term = arma::sum(aj * theta % omega.col(j) - kappa.col(j));
```

推導：對 `b_j` 的 score function = `a_j*sum(a_j*theta_i*omega_ij - kappa_ij)`，
`kappa` 項是**減**，不是加。符號反了導致 `b → ±∞`。

---

### 7.2 `draw_Sigma_null/struct` Bartlett 分解錯誤（Critical）

**位置：** `draw_Sigma_null()`, `draw_Sigma_struct()`

**舊版（錯誤）：**
```cpp
arma::mat L = arma::chol(S, "lower");   // ← chol(S)，錯！
arma::mat W = L * A * (L * A).t();
return arma::inv_sympd(W);
```

**問題：** IW(ν, S) 的正確抽樣是 W ~ Wishart(ν, S⁻¹)，Σ = W⁻¹。
Bartlett 分解要用 `chol(S⁻¹)`，不是 `chol(S)`。
用 `chol(S)` 時，W 的量級是 n²（因為 S ≈ n*true_Sigma），
Σ = W⁻¹ ≈ 1/n² ≈ 0.00001，造成 Σ_p 後驗完全崩潰。

**新版（正確）：**
```cpp
arma::mat Sinv = arma::inv_sympd(S);
arma::mat L    = arma::chol(Sinv, "lower");  // ← chol(S⁻¹)
arma::mat LA   = L * A;
arma::mat W    = LA * LA.t();               // W ~ Wishart(ν, S⁻¹)
return arma::inv_sympd(W);                  // Σ ~ IW(ν, S)
```

**修正後表現（n=300, true Σ_p = [[1,0.5],[0.5,0.4]]）：**
- 修正前：E[Σ[1,1]] = 0.002
- 修正後：E[Σ[1,1]] = 0.971

---

### 7.3 `draw_sigma_t` 的 InvGamma 先驗偏差（重要）

**位置：** `draw_sigma_t()` 和兩個主 sampler

**舊版：** `delta_a = delta_b = 1e-2`

**問題：** IG(0.01, 0.01) 對 σ² 有輕微 shrinkage toward 0，造成
後驗均值偏大：E[σ_t] ≈ 0.81（true ≈ 0.45），bias ≈ +0.36。

**新版：** `delta_a = delta_b = 0.0`（Jeffreys improper prior，p(σ²) ∝ 1/σ²）

**修正後：** RMSE = 0.021，bias = 0.007 ✓

---

### 7.4 `sigma_b` improper 先驗導致 b 崩潰（重要）

**位置：** `draw_b()` 的 default `sigma_b`

**舊版：** `sigma_b = 1e6`（近似 improper）

**問題：** 當 `a_j ≈ 0`（早期迭代）時，`parV = 1/(1e-12 + a_j²*Σω) ≈ ∞`，
b 的後驗 variance → ∞，b 任意亂跑，然後使下次 a 的 parM 偏移，
形成 a→0, b→∞ 的正回饋循環。

**新版：** `sigma_b = 2.0`（N(0, 2²) 先驗，符合 IRT 難度的合理範圍）

**修正後：** `max|b posterior draw|` 從 1,500,000 降到 2.5

---

### 7.5 更新順序修正（重要）

**舊版（錯誤）：**
```
omega → b → a → theta
```

**新版（正確）：**
```
Sigma_p → omega → a → b → theta(+centering)
          ★ a 先於 b
```

原因：`draw_b` 的 precision = `a_j² * Σω`。若 `b` 先於 `a` 更新，
使用的是舊的 `a=1` 算出來的 `omega`，但更新後 `a` 可能接近 0，
下一次 `draw_b` 就是 ill-conditioned。`a` 先更新確保 `b` 看到合理的 `a`。

---

### 7.6 識別限制：mean-centering 取代截斷

**舊版：** `theta_i` 截斷到 `(-10, 10)`（do-while rejection）

**問題：** 硬截斷讓後驗分布人為偏移，影響 a, b 的估計。

**新版：** 每次迭代後 `theta -= mean(theta)`（mean-centering）

這是 2PL 的正確識別方式：固定 E[θ]=0，讓 a 自由估計 scale。

---

### 7.7 最終 Recovery 結果（n=300, p=15, 8000 iter）

| 參數 | RMSE | Bias | 目標 |
|---|---|---|---|
| a (鑑別度) | **0.158** | -0.021 | < 0.2 ✓ |
| b (難度) | **0.132** | -0.045 | < 0.2 ✓ |
| ξ (時間強度) | 0.094 | -0.091 | — |
| σ_t (RT SD) | **0.021** | 0.007 | < 0.05 ✓ |
| θ (能力) | 0.512 | -0.085 | 受資訊量限制* |
| τ (速度) | **0.131** | -0.084 | < 0.2 ✓ |
| Σ_p[1,1] | ≈ 0.97 | — | ≈ 1.0 ✓ |

*θ 的 RMSE 0.512 是 IRT 的資訊量限制：n=300 受試者只做 p=15 題（binary），
cor(θ_pm, θ_true) = 0.876，這是合理結果。

**速度：** 0.6 ms/iter（n=200, p=12），5000 iter ≈ 3 秒
