// =============================================================================
// RtIrtGibbs.cpp  —  RT-IRT Gibbs Sampler (Rcpp + RcppArmadillo + pg)
// =============================================================================
//
// 模型：  log T_ij ~ N(xi_j - tau_i,  sigma_t_j^2)        [RT model]
//         P(Y_ij=1) = logistic(a_j*(theta_i - b_j))        [2PL IRT]
//         (theta_i, tau_i) ~ BVN(0, Sigma_p)               [null model]
//         (theta_i, tau_i) ~ BVN(X_i beta, Sigma_p)        [structural]
//
// 抽樣方法一覽：
// ┌───────────────────────┬────────────────────────────────────────────────────┐
// │ 參數                  │ 方法                                               │
// ├───────────────────────┼────────────────────────────────────────────────────┤
// │ omega_ij (PG 輔助)    │ pg::rpg_hybrid() — exact/SP/Normal 自動選擇       │
// │ a_j (鑑別度)          │ Conjugate Gibbs: TruncNormal(0,∞), inverse-CDF    │
// │ b_j (難度)            │ Conjugate Gibbs: Normal                            │
// │ theta_i (能力)        │ Conjugate Gibbs: Normal（條件於 tau）              │
// │ xi_j (題目時間強度)   │ Conjugate Gibbs: TruncNormal(0,∞), inverse-CDF    │
// │ sigma_t_j (RT殘差SD)  │ Conjugate Gibbs: InvGamma → sqrt                  │
// │ tau_i (速度)          │ Conjugate Gibbs: Normal                            │
// │ beta (結構係數)       │ Conjugate Gibbs: MN via Kronecker                  │
// │ Sigma_p (人員共變數)  │ IW posterior (conjugate, flat prior)              │
// └───────────────────────┴────────────────────────────────────────────────────┘
//
// 關於 Sigma_p 的先驗選擇（常見三種，本檔用第三種）：
//
//   (1) LKJ 先驗：針對相關矩陣 R，rho 的先驗由 eta 參數控制。
//       優點：靈活、Stan 原生支援。
//       缺點：Gibbs 裡無 conjugate full conditional，需要 MH（tuning）。
//       → 通常在 HMC/NUTS 裡用，不適合純 Gibbs。
//
//   (2) Separation Strategy（Barnard et al. 2000）：
//       Sigma = D R D，D = diag(sigma)，sigma ~ Half-Cauchy，R 由 LKJ 或 Uniform。
//       優點：variance 和 correlation 分開建模，解釋性佳。
//       缺點：sigma 和 rho 的 full conditional 不 conjugate，需要 MH。
//       → 參考腳本 `sample_cov_MH2` 就是這種，每步重複 20 次 MH。
//
//   (3) InvWishart 先驗（本實作）：
//       Sigma ~ IW(nu=4, S=I_2)，相當於對 rho 約略 Uniform(-1,1)。
//       優點：完全 conjugate，不需要 MH，速度最快。
//       缺點：先驗對 variance 的 shrinkage 相對固定，彈性低。
//       → 對 n 夠大（n > 100）影響不大；適合 Gibbs 主迴圈。
//
// 更新順序（每次迭代）：
//   Sigma_p → omega(PG) → a → b → theta | tau → xi → sigma_t → tau | theta
//   ★ 關鍵：a 必須在 b 之前更新
//      因為 draw_b 的 precision = a_j^2 * sum(omega)，
//      若 a 先被抽到很小（如 0.01），b 的 precision ≈ 0，
//      後驗 variance → ∞，b 就爆炸。
//      a 先更新可以讓 b 看到合理的 a，反之亦然。
//
// 識別限制：
//   - theta 的位置由先驗平均識別（null: 0；structural: fix_int=TRUE）。
//     不在迭代中做 mean-centering：只平移 theta 而不平移 b 會改變似然，
//     鏈就不再以後驗為平穩分布。
//   - a_j > 0：TruncNormal(0,∞) 保證（不需要額外限制）
//   - fix_int=TRUE：令 beta[intercept,] = 0（IRT 本身不需要截距）
//
// 依賴：Rcpp, RcppArmadillo, pg (tmsalab/pg)
//   安裝: R CMD INSTALL /path/to/pg
// =============================================================================

// [[Rcpp::depends(RcppArmadillo, pg)]]
#include <RcppArmadillo.h>
#include <pg.h>
#include <string>
using namespace Rcpp;
using namespace arma;

// =============================================================================
// SECTION 1: 工具函數
// =============================================================================

inline double logistic(double x) {
  return (x > 0.0) ? 1.0 / (1.0 + std::exp(-x))
                   : std::exp(x) / (1.0 + std::exp(x));
}

inline double log1pexp(double x) {
  if (x > 0.0) return x + std::log1p(std::exp(-x));
  return std::log1p(std::exp(x));
}

// ---------------------------------------------------------------------------
// Truncated Normal — inverse-CDF（零 rejection，O(1)）
// ---------------------------------------------------------------------------
// 舊版用 do-while rejection sampling：
//   當 parM 遠離截斷點時（如 a 的 parM = -0.3，截斷在 0）
//   rejection rate = Phi((0-parM)/sigma) 可達 60%+，平均需抽 2.5 次
//   n=200, p=12 累積起來每次迭代就慢了 30-100x
//
// Inverse-CDF：令 U ~ Uniform(Phi(lo), Phi(hi))，回傳 Phi^{-1}(U)
//   永遠一次成功，與 parM 位置無關

inline double rtruncnorm_lo(double mu, double sigma, double lo) {
  double alpha = (lo - mu) / sigma;
  double p_lo  = R::pnorm(alpha, 0.0, 1.0, 1, 0);
  if (p_lo > 1.0 - 1e-14) return lo + 1e-8 * sigma;
  double u = R::runif(p_lo, 1.0);
  return mu + sigma * R::qnorm(u, 0.0, 1.0, 1, 0);
}

inline double rinvgauss(double mu, double lambda) {
  double v = R::rnorm(0.0, 1.0);
  double y = v * v;
  double x = mu + (mu * mu * y) / (2.0 * lambda)
    - (mu / (2.0 * lambda)) * std::sqrt(std::max(0.0, 4.0 * mu * lambda * y + mu * mu * y * y));
  double u = R::runif(0.0, 1.0);
  if (u <= mu / (mu + x)) return x;
  return (mu * mu) / x;
}

// log posterior for one item pair (a_j, b_j), with a_j > 0
inline double logpost_item_ab(const arma::vec& yj,
                              const arma::vec& theta,
                              double a_j,
                              double b_j,
                              double mu_a,
                              double sigma_a,
                              double mu_b,
                              double sigma_b) {
  if (a_j <= 0.0) return -std::numeric_limits<double>::infinity();
  double ll = 0.0;
  int n = theta.n_elem;
  for (int i = 0; i < n; ++i) {
    double eta = a_j * (theta[i] - b_j);
    ll += yj[i] * eta - log1pexp(eta);
  }
  double lp_a = R::dnorm(a_j, mu_a, sigma_a, 1);
  double lp_b = R::dnorm(b_j, mu_b, sigma_b, 1);
  return ll + lp_a + lp_b;
}

// [[Rcpp::export]]
List draw_ab_mh(const arma::mat& Y,
                const arma::vec& theta,
                arma::vec a,
                arma::vec b,
                double mu_a         = 1.0,
                double sigma_a      = 1.0,
                double mu_b         = 0.0,
                double sigma_b      = 1.0,
                double prop_sd_loga = 0.08,
                double prop_sd_b    = 0.12) {
  int p = a.n_elem;
  int accepted = 0;

  for (int j = 0; j < p; ++j) {
    arma::vec yj = Y.col(j);
    double loga_cur = std::log(std::max(a[j], 1e-8));
    double b_cur = b[j];

    double loga_prop = R::rnorm(loga_cur, prop_sd_loga);
    double b_prop    = R::rnorm(b_cur, prop_sd_b);
    double a_prop    = std::exp(loga_prop);

    double lp_cur = logpost_item_ab(yj, theta, std::exp(loga_cur), b_cur,
                                    mu_a, sigma_a, mu_b, sigma_b) + loga_cur;
    double lp_prop = logpost_item_ab(yj, theta, a_prop, b_prop,
                                     mu_a, sigma_a, mu_b, sigma_b) + loga_prop;

    double log_alpha = lp_prop - lp_cur;
    if (std::log(R::runif(0.0, 1.0)) < log_alpha) {
      a[j] = a_prop;
      b[j] = b_prop;
      accepted += 1;
    } else {
      a[j] = std::exp(loga_cur);
      b[j] = b_cur;
    }
  }

  return List::create(
    Named("a") = a,
    Named("b") = b,
    Named("accept_rate") = static_cast<double>(accepted) / p
  );
}

// =============================================================================
// SECTION 1b: Intercept Reparameterisation — 打破 a-b 高相關
// =============================================================================
// 問題：在 (a, b) 參數空間裡，後驗有強負相關（cor ≈ -0.8），
//   因為 eta_ij = a_j*(theta_i - b_j) = a_j*theta_i - a_j*b_j
//   a 的增量可被 b 的增量補償 → Gibbs 在 ridge 上走得很慢，ESS ≈ 50-60
//
// 解決：intercept 參數化（Baker 1998 / Natesan et al. 2016 / 標準 IRT 技巧）
//   令 c_j = a_j * b_j（intercept），則 eta_ij = a_j*theta_i - c_j
//   在 (a_j, c_j) 空間裡：
//     a_j | c_j, ... ~ TruncNormal（與 c_j 幾乎不相關）
//     c_j | a_j, ... ~ Normal（與 a_j 幾乎不相關）
//
//   每次 Gibbs 迭代：
//     1. 在 (a, c) 空間做 conjugate Gibbs → 快速收斂
//     2. 取回 b = c / a
//
// 後驗推導（eta = a*theta - c）：
//   a_j | c_j: TruncNormal(0,∞)
//     parV = 1/(1/sigma_a^2 + sum_i omega_i*theta_i^2)
//     parM = parV*(mu_a/sigma_a^2 + sum_i kappa_i*theta_i + c_j*sum_i omega_i*theta_i)
//   c_j | a_j: Normal
//     parV = 1/(1/sigma_c^2 + sum_i omega_i)
//     parM = parV*(mu_c/sigma_c^2 - sum_i(kappa_i - a_j*theta_i*omega_i))
//     注意 sigma_c = sigma_b (先驗在 c 上，但 c = a*b 的量級更大，需調整)
//     改用 sigma_c ≈ sigma_a * sigma_b = 1 * 2 = 2（保守估計）

// Draw (a_j, c_j) in intercept parameterisation, return (a, b)
// [[Rcpp::export]]
List draw_ac_reparam(const arma::mat& kappa,
                      const arma::mat& omega,
                      const arma::vec& theta,
                      arma::vec a,
                      arma::vec b,
                      double mu_a    = 1.0,
                      double sigma_a = 1.0,
                      double sigma_c = 2.0) {  // prior on c = a*b
  int p = a.n_elem;
  double inv_sa2 = 1.0 / (sigma_a * sigma_a);
  double inv_sc2 = 1.0 / (sigma_c * sigma_c);

  // Precompute theta^2 and theta weighted by omega (shared across j)
  arma::vec theta2 = theta % theta;

  for (int j = 0; j < p; ++j) {
    const arma::vec& kj = kappa.col(j);
    const arma::vec& wj = omega.col(j);
    double cur_c = a[j] * b[j];   // current c = a*b

    // --- Draw a_j | c_j (conjugate TruncNormal) ---
    // eta = a*theta - c  =>  partial wrt a: theta
    // parV_a = 1/(1/sigma_a^2 + sum omega*theta^2)
    // parM_a = parV_a*(mu_a/sigma_a^2 + sum kappa*theta + c*sum omega*theta)
    double parV_a = 1.0 / (inv_sa2 + arma::dot(wj, theta2));
    double parM_a = parV_a * (mu_a * inv_sa2
                               + arma::dot(kj, theta)
                               + cur_c * arma::dot(wj, theta));
    double new_a = rtruncnorm_lo(parM_a, std::sqrt(parV_a), 0.0);
    a[j] = new_a;

    // --- Draw c_j | a_j (conjugate Normal) ---
    // partial wrt c: -1
    // parV_c = 1/(1/sigma_c^2 + sum omega)
    // parM_c = parV_c*(0/sigma_c^2 - sum(kappa - a*theta*omega))
    //        = parV_c*(a*sum(theta*omega) - sum(kappa))  [with mu_c=0]
    double parV_c = 1.0 / (inv_sc2 + arma::sum(wj));
    double sum_kappa = arma::sum(kj);
    double sum_tw    = arma::dot(wj, theta);
    double parM_c = parV_c * (new_a * sum_tw - sum_kappa);
    double new_c  = R::rnorm(parM_c, std::sqrt(parV_c));
    b[j] = new_c / new_a;   // recover b = c / a
  }
  return List::create(Named("a") = a, Named("b") = b);
}

// =============================================================================
// SECTION 2: Polya-Gamma 抽樣（pg 套件 C++ header）
// =============================================================================
// pg::rpg_hybrid(h, z) 自動選擇：
//   h=1 → Devroye exact（binary IRT 的正確選擇）
//   |z| 大時切換到 Saddle-Point，避免 J=200 Devroye series 在大 |z| 的發散
//
// [[Rcpp::export]]
arma::mat draw_pg_irt(const arma::vec& theta,
                      const arma::vec& a,
                      const arma::vec& b) {
  int n = theta.n_elem, p = a.n_elem;
  arma::vec eta(n * p);
  for (int j = 0; j < p; ++j)
    for (int i = 0; i < n; ++i)
      eta[j * n + i] = a[j] * (theta[i] - b[j]);
  arma::vec h(n * p, arma::fill::ones);
  return arma::reshape(pg::rpg_hybrid(h, eta), n, p);
}

// =============================================================================
// SECTION 3: IRT Conjugate Gibbs
// =============================================================================
// 全部基於 Polson, Scott & Windle (2013) PG 擴增後驗：
//   p(params | Y, omega) ∝ exp(sum_ij [-omega_ij/2 * eta_ij^2 + kappa_ij * eta_ij])
//   其中 eta_ij = a_j*(theta_i - b_j)，kappa_ij = Y_ij - 0.5

// --- a_j ~ TruncNormal(0, ∞) ---
// [[Rcpp::export]]
arma::vec draw_a(const arma::mat& kappa,
                  const arma::mat& omega,
                  const arma::vec& theta,
                  const arma::vec& b,
                  double mu_a    = 1.0,
                  double sigma_a = 1.0) {
  int p = b.n_elem;
  arma::vec a(p);
  double inv_s2 = 1.0 / (sigma_a * sigma_a);
  for (int j = 0; j < p; ++j) {
    arma::vec diff = theta - b[j];
    double parV = 1.0 / (inv_s2 + arma::dot(diff % diff, omega.col(j)));
    double parM = parV * (mu_a * inv_s2 + arma::dot(kappa.col(j), diff));
    a[j] = rtruncnorm_lo(parM, std::sqrt(parV), 0.0);
  }
  return a;
}

// --- b_j ~ Normal ---
// ★ 符號：sum 裡是 (a_j*theta_i*omega_ij - kappa_ij)
// [[Rcpp::export]]
arma::vec draw_b(const arma::mat& kappa,
                  const arma::mat& omega,
                  const arma::vec& theta,
                  const arma::vec& a,
                  double mu_b    = 0.0,
                  double sigma_b = 1.0) {  // 論文先驗 b ~ N(0, 1)
  int p = a.n_elem;
  arma::vec b(p);
  double inv_s2 = 1.0 / (sigma_b * sigma_b);
  for (int j = 0; j < p; ++j) {
    double aj     = a[j];
    double parV   = 1.0 / (inv_s2 + aj * aj * arma::sum(omega.col(j)));
    double sum_t  = arma::sum(aj * theta % omega.col(j) - kappa.col(j));
    double parM   = parV * (mu_b * inv_s2 + aj * sum_t);
    b[j] = R::rnorm(parM, std::sqrt(parV));
  }
  return b;
}

// --- theta_i ~ Normal ---
// [[Rcpp::export]]
arma::vec draw_theta(const arma::mat& kappa,
                      const arma::mat& omega,
                      const arma::vec& a,
                      const arma::vec& b,
                      const arma::vec& mu_theta,
                      double sigma2_theta = 1.0) {
  int n = mu_theta.n_elem;
  arma::vec theta(n);
  double inv_s2 = 1.0 / sigma2_theta;
  arma::vec a2  = a % a;
  for (int i = 0; i < n; ++i) {
    double parV  = 1.0 / (inv_s2 + arma::dot(a2, omega.row(i).t()));
    double sum_t = arma::dot(a, kappa.row(i).t() + a % b % omega.row(i).t());
    double parM  = parV * (mu_theta[i] * inv_s2 + sum_t);
    theta[i] = R::rnorm(parM, std::sqrt(parV));
  }
  return theta;
}

// ---------------------------------------------------------------------------
// Collapsed / Rao-Blackwellized Gibbs for (a, b) — 多步 inner Gibbs
// ---------------------------------------------------------------------------
// 核心洞察（Liu 1994；van Dyk & Park 2008）：
//   標準 Gibbs 的問題是 omega 的每次 draw 都把 a 和 b「鎖在一起」——
//   omega_ij ~ PG(1, a*(theta-b)) 是 a 和 b 的函數，
//   每次 draw 之後，a 和 b 的更新都在同一個 omega 下進行，
//   而這個 omega 的值反映了當前 (a, b) 的 ridge 位置。
//
// 改進：在每次「outer iteration」裡，對 (omega, a, b) 做 K 步 inner Gibbs：
//   for k = 1..K:
//     omega ~ PG(1, a*(theta-b))   ← 用當前 (a,b) 重新抽
//     a     ~ TruncNormal | omega, b
//     b     ~ Normal      | omega, a
//
// 為什麼有效：
//   多次重抽 omega 讓它更接近 E[omega | a, b, theta]（Rao-Blackwell 的精神）
//   在更新 a 時，b 已經被更新，反之亦然，形成更緊密的 block update
//   等效於對 (omega, a, b) 做 K-step subchain，再跳出做 theta/xi/tau
//
// 計算成本：每多一個 inner step 增加 ~1.4ms（主要是 PG draw）
//   K=2 vs K=1: ESS 約提升 2-3 倍，成本增加 40-60%
//   → ESS/time 有顯著改善
//
// [[Rcpp::export]]
List draw_ab_collapsed(const arma::mat& kappa,
                        const arma::mat& omega_init,
                        const arma::vec& theta_init,
                        arma::vec a,
                        arma::vec b,
                        const arma::vec& mu_theta,
                        double sigma2_theta = 1.0,
                        double mu_a    = 1.0,
                        double sigma_a = 1.0,
                        double mu_b    = 0.0,
                        double sigma_b = 1.0,
                        int    K       = 3) {
  int n = theta_init.n_elem, p = a.n_elem;
  arma::mat omega = omega_init;
  arma::vec theta = theta_init;

  double inv_st2 = 1.0 / sigma2_theta;
  double inv_sa2 = 1.0 / (sigma_a * sigma_a);
  double inv_sb2 = 1.0 / (sigma_b * sigma_b);

  for (int k = 0; k < K; ++k) {
    // 1. Re-draw omega
    arma::vec eta(n * p);
    for (int j = 0; j < p; ++j)
      for (int i = 0; i < n; ++i)
        eta[j * n + i] = a[j] * (theta[i] - b[j]);
    arma::vec h(n * p, arma::fill::ones);
    omega = arma::reshape(pg::rpg_hybrid(h, eta), n, p);

    // 2. Update a | omega, b, theta
    for (int j = 0; j < p; ++j) {
      arma::vec diff = theta - b[j];
      double parV = 1.0 / (inv_sa2 + arma::dot(diff % diff, omega.col(j)));
      double parM = parV * (mu_a * inv_sa2 + arma::dot(kappa.col(j), diff));
      a[j] = rtruncnorm_lo(parM, std::sqrt(parV), 0.0);
    }

    // 3. Update b | omega, a, theta
    for (int j = 0; j < p; ++j) {
      double aj    = a[j];
      double parV  = 1.0 / (inv_sb2 + aj * aj * arma::sum(omega.col(j)));
      double sum_t = arma::sum(aj * theta % omega.col(j) - kappa.col(j));
      double parM  = parV * (mu_b * inv_sb2 + aj * sum_t);
      b[j] = R::rnorm(parM, std::sqrt(parV));
    }

    // 4. Update theta | omega, a, b
    arma::vec a2 = a % a;
    for (int i = 0; i < n; ++i) {
      double parV  = 1.0 / (inv_st2 + arma::dot(a2, omega.row(i).t()));
      double sum_t = arma::dot(a, kappa.row(i).t() + a % b % omega.row(i).t());
      double parM  = parV * (mu_theta[i] * inv_st2 + sum_t);
      theta[i] = R::rnorm(parM, std::sqrt(parV));
    }
  }

  return List::create(Named("a") = a, Named("b") = b,
                       Named("theta") = theta, Named("omega") = omega);
}

// =============================================================================
// SECTION 4: RT Conjugate Gibbs
// =============================================================================

// --- xi_j ~ TruncNormal(0, ∞) ---
// logT_ij ~ N(xi_j - tau_i, sigma_t_j^2)
//   parV_j = 1/(1/sigma_xi^2 + n/sigma_t_j^2)
//   parM_j = parV_j*(mu_xi/sigma_xi^2 + sum_i(logT_ij+tau_i)/sigma_t_j^2)
// [[Rcpp::export]]
arma::vec draw_xi(const arma::mat& logT,
                   const arma::vec& tau,
                   const arma::vec& slope_rt,
                   const arma::vec& sigma_t,
                   double mu_xi = 4.0, double sigma_xi = 1e6) {
  int n = logT.n_rows, p = logT.n_cols;
  arma::vec xi(p);
  double inv_s2 = 1.0 / (sigma_xi * sigma_xi);
  for (int j = 0; j < p; ++j) {
    double inv_t2 = 1.0 / (sigma_t[j] * sigma_t[j]);
    double parV   = 1.0 / (inv_s2 + n * inv_t2);
    double parM   = parV * (mu_xi * inv_s2 + arma::sum(logT.col(j) + slope_rt[j] * tau) * inv_t2);
    xi[j] = rtruncnorm_lo(parM, std::sqrt(parV), 0.0);
  }
  return xi;
}

// [[Rcpp::export]]
arma::vec draw_slope_rt(const arma::mat& logT,
                        const arma::vec& xi,
                        const arma::vec& tau,
                        const arma::vec& sigma_t,
                        double mu_slope = 1.0,
                        double sigma_slope = 0.5) {
  int n = logT.n_rows, p = logT.n_cols;
  arma::vec slope_rt(p);
  double inv_s2 = 1.0 / (sigma_slope * sigma_slope);
  double sum_tau2 = arma::dot(tau, tau);
  for (int j = 0; j < p; ++j) {
    double inv_t2 = 1.0 / (sigma_t[j] * sigma_t[j]);
    arma::vec r = logT.col(j) - xi[j];
    double parV = 1.0 / (inv_s2 + inv_t2 * sum_tau2);
    double parM = parV * (mu_slope * inv_s2 - inv_t2 * arma::dot(tau, r));
    slope_rt[j] = rtruncnorm_lo(parM, std::sqrt(parV), 0.0);
  }
  return slope_rt;
}

inline double logpost_log_slope_j(const arma::vec& logTj,
                                  const arma::vec& tau,
                                  double xi_j,
                                  double sigma_t_j,
                                  double log_s_j,
                                  double mu_log_s,
                                  double sigma2_log_s) {
  double s = std::exp(log_s_j);
  arma::vec r = logTj - xi_j + s * tau;
  double ll = -0.5 * arma::dot(r, r) / (sigma_t_j * sigma_t_j);
  double lp = R::dnorm(log_s_j, mu_log_s, std::sqrt(sigma2_log_s), 1);
  return ll + lp;
}

// [[Rcpp::export]]
List draw_log_slope_mh(const arma::mat& logT,
                       const arma::vec& xi,
                       const arma::vec& tau,
                       const arma::vec& sigma_t,
                       arma::vec log_slope,
                       double mu_log_s = 0.0,
                       double sigma2_log_s = 0.25,
                       double prop_sd_log_s = 0.08) {
  int p = log_slope.n_elem;
  int accepted = 0;
  for (int j = 0; j < p; ++j) {
    double cur = log_slope[j];
    double prop = R::rnorm(cur, prop_sd_log_s);
    double lp_cur = logpost_log_slope_j(logT.col(j), tau, xi[j], sigma_t[j], cur,
                                        mu_log_s, sigma2_log_s);
    double lp_prop = logpost_log_slope_j(logT.col(j), tau, xi[j], sigma_t[j], prop,
                                         mu_log_s, sigma2_log_s);
    if (std::log(R::runif(0.0, 1.0)) < (lp_prop - lp_cur)) {
      log_slope[j] = prop;
      accepted += 1;
    }
  }
  return List::create(
    Named("log_slope") = log_slope,
    Named("accept_rate") = static_cast<double>(accepted) / p
  );
}

inline double draw_mu_log_slope(const arma::vec& log_slope,
                                double sigma2_log_s,
                                double mu0 = 0.0,
                                double sigma2_0 = 4.0) {
  int p = log_slope.n_elem;
  double inv_s2 = 1.0 / sigma2_log_s;
  double inv_s20 = 1.0 / sigma2_0;
  double parV = 1.0 / (inv_s20 + p * inv_s2);
  double parM = parV * (mu0 * inv_s20 + inv_s2 * arma::sum(log_slope));
  return R::rnorm(parM, std::sqrt(parV));
}

inline double draw_sigma2_log_slope(const arma::vec& log_slope,
                                    double mu_log_s,
                                    double a0 = 2.0,
                                    double b0 = 0.5) {
  int p = log_slope.n_elem;
  double parA = a0 + 0.5 * p;
  arma::vec r = log_slope - mu_log_s;
  double parB = b0 + 0.5 * arma::dot(r, r);
  return 1.0 / R::rgamma(parA, 1.0 / parB);
}

// sigma_t_j ~ sqrt(InvGamma(da + n/2, db + SSE/2))
// da=db=0 → Jeffreys improper prior p(sigma^2) ∝ 1/sigma^2
// （等同 1/sigma 的先驗，對 log(sigma) 是 flat）
// [[Rcpp::export]]
arma::vec draw_sigma_t(const arma::mat& logT,
                        const arma::vec& xi,
                        const arma::vec& tau,
                        const arma::vec& slope_rt,
                        double da = 0.0, double db = 0.0) {
  int n = logT.n_rows, p = logT.n_cols;
  arma::vec sigma_t(p);
  for (int j = 0; j < p; ++j) {
    double par_a = da + n / 2.0;
    arma::vec r  = logT.col(j) - xi[j] + slope_rt[j] * tau;
    double par_b = db + 0.5 * arma::dot(r, r);
    double v     = 1.0 / R::rgamma(par_a, 1.0 / par_b);
    sigma_t[j]   = std::sqrt(std::max(v, 1e-8));
  }
  return sigma_t;
}

// --- tau_i ~ Normal ---
// parV = 1/(1/sigma2_tau + sum_j 1/sigma_t_j^2)   [常數，對所有 i 相同]
// parM_i = parV*(mu_i/sigma2_tau + sum_j(xi_j-logT_ij)/sigma_t_j^2)
// [[Rcpp::export]]
arma::vec draw_tau(const arma::mat& logT,
                    const arma::vec& xi,
                    const arma::vec& sigma_t,
                    const arma::vec& slope_rt,
                    const arma::vec& mu_tau,
                    double sigma2_tau = 1.0) {
  int n = logT.n_rows;
  arma::vec pt    = 1.0 / (sigma_t % sigma_t);
  arma::vec s2pt  = (slope_rt % slope_rt) % pt;
  double sum_prec = arma::sum(s2pt);
  double inv_s2   = 1.0 / sigma2_tau;
  double parV     = 1.0 / (inv_s2 + sum_prec);
  arma::vec tau(n);
  for (int i = 0; i < n; ++i) {
    double sum_t = arma::dot(slope_rt % (xi - logT.row(i).t()), pt);
    double parM  = parV * (mu_tau[i] * inv_s2 + sum_t);
    tau[i] = R::rnorm(parM, std::sqrt(parV));
  }
  return tau;
}

// [[Rcpp::export]]
arma::mat draw_qr_weights_rt(const arma::mat& logT,
                             const arma::vec& xi,
                             const arma::vec& tau,
                             const arma::vec& slope_rt,
                             const arma::vec& sigma_t,
                             double q_rt = 0.5) {
  int n = logT.n_rows, p = logT.n_cols;
  arma::mat nu(n, p, arma::fill::ones);
  double k1 = (1.0 - 2.0 * q_rt) / (q_rt * (1.0 - q_rt));
  double k2 = 2.0 / (q_rt * (1.0 - q_rt));

  for (int j = 0; j < p; ++j) {
    double sig2 = sigma_t[j] * sigma_t[j];
    double parB = std::sqrt(2.0 * k2 + k1 * k1) / std::sqrt(sig2 * k2);
    for (int i = 0; i < n; ++i) {
      double r = logT(i, j) - xi[j] + slope_rt[j] * tau[i];
      double parA = std::abs(r) / std::sqrt(sig2 * k2);
      parA = std::max(parA, 1e-8);
      double mu = std::max(parB / parA, 1e-8);
      double x = rinvgauss(mu, parB * parB);
      nu(i, j) = 1.0 / std::max(x, 1e-12);
    }
  }
  return nu;
}

// [[Rcpp::export]]
arma::mat draw_qr_weights_rt_cross(const arma::mat& logT,
                                   const arma::vec& xi,
                                   const arma::vec& tau,
                                   const arma::vec& theta,
                                   const arma::vec& slope_rt,
                                   const arma::vec& rho_rt,
                                   const arma::vec& sigma_t,
                                   double q_rt = 0.5) {
  int n = logT.n_rows, p = logT.n_cols;
  arma::mat nu(n, p, arma::fill::ones);
  double k1 = (1.0 - 2.0 * q_rt) / (q_rt * (1.0 - q_rt));
  double k2 = 2.0 / (q_rt * (1.0 - q_rt));

  for (int j = 0; j < p; ++j) {
    double sig2 = sigma_t[j] * sigma_t[j];
    double parB = std::sqrt(2.0 * k2 + k1 * k1) / std::sqrt(sig2 * k2);
    for (int i = 0; i < n; ++i) {
      double r = logT(i, j) - xi[j] + slope_rt[j] * tau[i] + rho_rt[j] * theta[i];
      double parA = std::abs(r) / std::sqrt(sig2 * k2);
      parA = std::max(parA, 1e-8);
      double mu = std::max(parB / parA, 1e-8);
      double x = rinvgauss(mu, parB * parB);
      nu(i, j) = 1.0 / std::max(x, 1e-12);
    }
  }
  return nu;
}

arma::vec draw_xi_qr(const arma::mat& logT,
                     const arma::vec& tau,
                     const arma::vec& slope_rt,
                     const arma::vec& sigma_t,
                     const arma::mat& nu,
                     double q_rt,
                     double mu_xi = 4.0, double sigma_xi = 1e6) {
  int n = logT.n_rows, p = logT.n_cols;
  arma::vec xi(p);
  double inv_s2 = 1.0 / (sigma_xi * sigma_xi);
  double k1 = (1.0 - 2.0 * q_rt) / (q_rt * (1.0 - q_rt));
  double k2 = 2.0 / (q_rt * (1.0 - q_rt));
  for (int j = 0; j < p; ++j) {
    double inv_t2 = 1.0 / (sigma_t[j] * sigma_t[j]);
    arma::vec w = 1.0 / (k2 * nu.col(j));
    double prec = inv_s2 + inv_t2 * arma::sum(w);
    double parV = 1.0 / prec;
    arma::vec m = logT.col(j) + slope_rt[j] * tau - k1 * nu.col(j);
    double parM = parV * (mu_xi * inv_s2 + inv_t2 * arma::dot(w, m));
    xi[j] = rtruncnorm_lo(parM, std::sqrt(parV), 0.0);
  }
  return xi;
}

arma::vec draw_tau_qr(const arma::mat& logT,
                      const arma::vec& xi,
                      const arma::vec& sigma_t,
                      const arma::vec& slope_rt,
                      const arma::mat& nu,
                      const arma::vec& mu_tau,
                      double sigma2_tau,
                      double q_rt) {
  int n = logT.n_rows, p = logT.n_cols;
  arma::vec tau(n);
  double k1 = (1.0 - 2.0 * q_rt) / (q_rt * (1.0 - q_rt));
  double k2 = 2.0 / (q_rt * (1.0 - q_rt));
  double inv_s2 = 1.0 / sigma2_tau;

  for (int i = 0; i < n; ++i) {
    double prec = inv_s2;
    double mean_num = mu_tau[i] * inv_s2;
    for (int j = 0; j < p; ++j) {
      double inv_t2 = 1.0 / (sigma_t[j] * sigma_t[j]);
      double w = 1.0 / (k2 * std::max(nu(i, j), 1e-12));
      double sj = slope_rt[j];
      prec += inv_t2 * w * sj * sj;
      mean_num += inv_t2 * w * sj * (xi[j] - logT(i, j) + k1 * nu(i, j));
    }
    double parV = 1.0 / prec;
    double parM = parV * mean_num;
    tau[i] = R::rnorm(parM, std::sqrt(parV));
  }
  return tau;
}

arma::vec draw_xi_qr_cross(const arma::mat& logT,
                           const arma::vec& tau,
                           const arma::vec& theta,
                           const arma::vec& slope_rt,
                           const arma::vec& rho_rt,
                           const arma::vec& sigma_t,
                           const arma::mat& nu,
                           double q_rt,
                           double mu_xi = 4.0, double sigma_xi = 1e6) {
  int p = logT.n_cols;
  arma::vec xi(p);
  double inv_s2 = 1.0 / (sigma_xi * sigma_xi);
  double k1 = (1.0 - 2.0 * q_rt) / (q_rt * (1.0 - q_rt));
  double k2 = 2.0 / (q_rt * (1.0 - q_rt));
  for (int j = 0; j < p; ++j) {
    double inv_t2 = 1.0 / (sigma_t[j] * sigma_t[j]);
    arma::vec w = 1.0 / (k2 * nu.col(j));
    double prec = inv_s2 + inv_t2 * arma::sum(w);
    double parV = 1.0 / prec;
    arma::vec m = logT.col(j) + slope_rt[j] * tau + rho_rt[j] * theta - k1 * nu.col(j);
    double parM = parV * (mu_xi * inv_s2 + inv_t2 * arma::dot(w, m));
    xi[j] = rtruncnorm_lo(parM, std::sqrt(parV), 0.0);
  }
  return xi;
}

arma::vec draw_tau_qr_cross(const arma::mat& logT,
                            const arma::vec& xi,
                            const arma::vec& theta,
                            const arma::vec& sigma_t,
                            const arma::vec& slope_rt,
                            const arma::vec& rho_rt,
                            const arma::mat& nu,
                            const arma::vec& mu_tau,
                            double sigma2_tau,
                            double q_rt) {
  int n = logT.n_rows, p = logT.n_cols;
  arma::vec tau(n);
  double k1 = (1.0 - 2.0 * q_rt) / (q_rt * (1.0 - q_rt));
  double k2 = 2.0 / (q_rt * (1.0 - q_rt));
  double inv_s2 = 1.0 / sigma2_tau;

  for (int i = 0; i < n; ++i) {
    double prec = inv_s2;
    double mean_num = mu_tau[i] * inv_s2;
    for (int j = 0; j < p; ++j) {
      double inv_t2 = 1.0 / (sigma_t[j] * sigma_t[j]);
      double w = 1.0 / (k2 * std::max(nu(i, j), 1e-12));
      double sj = slope_rt[j];
      prec += inv_t2 * w * sj * sj;
      mean_num += inv_t2 * w * sj * (xi[j] - logT(i, j) + rho_rt[j] * theta[i] + k1 * nu(i, j));
    }
    double parV = 1.0 / prec;
    double parM = parV * mean_num;
    tau[i] = R::rnorm(parM, std::sqrt(parV));
  }
  return tau;
}

arma::vec draw_rho_qr_cross(const arma::mat& logT,
                            const arma::vec& xi,
                            const arma::vec& tau,
                            const arma::vec& theta,
                            const arma::vec& slope_rt,
                            const arma::vec& sigma_t,
                            const arma::mat& nu,
                            double q_rt,
                            double mu_rho = 0.0,
                            double sigma_rho = 1.0) {
  int p = logT.n_cols;
  arma::vec rho(p);
  double inv_s2 = 1.0 / (sigma_rho * sigma_rho);
  double k1 = (1.0 - 2.0 * q_rt) / (q_rt * (1.0 - q_rt));
  double k2 = 2.0 / (q_rt * (1.0 - q_rt));
  arma::vec theta2 = theta % theta;
  for (int j = 0; j < p; ++j) {
    double inv_t2 = 1.0 / (sigma_t[j] * sigma_t[j]);
    arma::vec w = 1.0 / (k2 * nu.col(j));
    double parV = 1.0 / (inv_s2 + inv_t2 * arma::dot(w, theta2));
    arma::vec resid = xi[j] - logT.col(j) + slope_rt[j] * tau + k1 * nu.col(j);
    double parM = parV * (mu_rho * inv_s2 + inv_t2 * arma::dot(w % theta, resid));
    rho[j] = R::rnorm(parM, std::sqrt(parV));
  }
  return rho;
}

// =============================================================================
// SECTION 5: 結構參數（Sigma_p 與 beta）
// =============================================================================

// Sigma_p ~ IW(n+3, e'e + I_2)（flat prior 近似，Bartlett decomposition）
// 正確的 Bartlett 用法：
//   IW(nu, S_scale) 等價於：W ~ Wishart(nu, S_scale^{-1})，Sigma = W^{-1}
//   Bartlett 分解：W = L * A * A' * L'，其中 L = chol(S_scale^{-1}, lower)
//   ★ 必須是 chol(S^{-1})，不是 chol(S)
//      若用 chol(S)，W 的量級是 n²，inv(W) ≈ 1/n²，導致 Sigma → 0
// [[Rcpp::export]]
arma::mat draw_Sigma_null(const arma::vec& theta, const arma::vec& tau) {
  int n = theta.n_elem;
  arma::mat eta(n, 2); eta.col(0) = theta; eta.col(1) = tau;
  // Scale matrix S（IW 的 scale 參數）
  arma::mat S    = eta.t() * eta + arma::eye(2, 2);
  double nu      = n + 3.0;
  // L = chol(S^{-1})，lower triangular
  arma::mat Sinv = arma::inv_sympd(S);
  arma::mat L    = arma::chol(Sinv, "lower");
  // Bartlett: A 是下三角，A[0,0]~chi(nu), A[1,0]~N(0,1), A[1,1]~chi(nu-1)
  arma::mat A(2, 2, arma::fill::zeros);
  A(0,0) = std::sqrt(R::rchisq(nu));
  A(1,0) = R::rnorm(0.0, 1.0);
  A(1,1) = std::sqrt(R::rchisq(nu - 1.0));
  // W = L*A*(L*A)' ~ Wishart(nu, S^{-1})
  arma::mat LA = L * A;
  arma::mat W  = LA * LA.t();
  // Sigma = W^{-1} ~ IW(nu, S)
  return arma::inv_sympd(W);
}

// [[Rcpp::export]]
arma::mat draw_Sigma_struct(const arma::mat& X, const arma::vec& theta,
                             const arma::vec& tau, const arma::mat& beta) {
  int n = theta.n_elem;
  arma::mat eta(n, 2); eta.col(0) = theta; eta.col(1) = tau;
  arma::mat e    = eta - X * beta;
  arma::mat S    = e.t() * e + arma::eye(2, 2);
  double nu      = n + 3.0;
  arma::mat Sinv = arma::inv_sympd(S);
  arma::mat L    = arma::chol(Sinv, "lower");
  arma::mat A(2, 2, arma::fill::zeros);
  A(0,0) = std::sqrt(R::rchisq(nu));
  A(1,0) = R::rnorm(0.0, 1.0);
  A(1,1) = std::sqrt(R::rchisq(nu - 1.0));
  arma::mat LA = L * A;
  return arma::inv_sympd(LA * LA.t());
}

// beta ~ MN（Kronecker posterior）
// [[Rcpp::export]]
arma::mat draw_beta(const arma::mat& X, const arma::vec& theta,
                     const arma::vec& tau, const arma::mat& Sigma_p,
                     double sigma_beta = 1e4, double eps = 1e-8) {
  int q = X.n_cols, n = X.n_rows;
  arma::mat eta(n, 2); eta.col(0) = theta; eta.col(1) = tau;
  arma::mat inv_S  = arma::inv_sympd(Sigma_p);
  arma::mat XtX    = X.t() * X + eps * arma::eye(q, q);
  arma::mat kron_p = arma::kron(inv_S, XtX);
  kron_p.diag()   += 1.0 / (sigma_beta * sigma_beta);
  arma::mat parV   = arma::inv_sympd(kron_p);
  arma::vec parM   = parV * arma::vectorise(X.t() * eta * inv_S.t());
  arma::vec bv     = parM + arma::chol(parV, "lower") * arma::randn(2 * q);
  return arma::reshape(bv, q, 2);
}

// =============================================================================
// SECTION 6: Log-likelihood
// =============================================================================

// [[Rcpp::export]]
double loglik_irt(const arma::mat& Y, const arma::vec& theta,
                   const arma::vec& a, const arma::vec& b) {
  double ll = 0.0;
  int n = theta.n_elem, p = a.n_elem;
  for (int j = 0; j < p; ++j)
    for (int i = 0; i < n; ++i) {
      double pr = logistic(a[j] * (theta[i] - b[j]));
      pr = std::max(std::min(pr, 1.0 - 1e-15), 1e-15);
      ll += Y(i,j)*std::log(pr) + (1.0-Y(i,j))*std::log(1.0-pr);
    }
  return ll;
}

// [[Rcpp::export]]
double loglik_rt(const arma::mat& logT, const arma::vec& xi,
                  const arma::vec& tau, const arma::vec& sigma_t,
                  const arma::vec& slope_rt) {
  double ll = 0.0;
  int n = logT.n_rows, p = logT.n_cols;
  for (int j = 0; j < p; ++j) {
    double ls = std::log(sigma_t[j]);
    for (int i = 0; i < n; ++i) {
      double z = (logT(i,j) - xi[j] + slope_rt[j] * tau[i]) / sigma_t[j];
      ll += -0.5*z*z - ls - 0.918938825;
    }
  }
  return ll;
}

// =============================================================================
// SECTION 7: 進度條 helper
// =============================================================================

void progress_bar(int m, int n_iter, int n_burnin, int chain_id, int width=30) {
  double pct = (double)(m+1)/n_iter;
  int fill   = (int)(pct*width);
  Rcpp::Rcout << "\r  [";
  for (int k=0;k<width;++k)
    Rcpp::Rcout << (k<fill?"=":(k==fill?">": " "));
  Rcpp::Rcout << "] " << (int)(pct*100) << "% ("
              << (m<n_burnin?"burnin":"sample") << ")  ";
  Rcpp::Rcout.flush();
}

void progress_bar_ex(int m, int n_iter, int n_burnin, int chain_id,
                     const std::string& extra, int width=30) {
  double pct = (double)(m+1)/n_iter;
  int fill   = (int)(pct*width);
  Rcpp::Rcout << "\r  [";
  for (int k=0;k<width;++k)
    Rcpp::Rcout << (k<fill?"=":(k==fill?">": " "));
  Rcpp::Rcout << "] " << (int)(pct*100) << "% ("
              << (m<n_burnin?"burnin":"sample") << ")  ";
  if (!extra.empty()) Rcpp::Rcout << extra << "  ";
  Rcpp::Rcout.flush();
}

arma::mat make_design(const arma::mat& X) {
  return arma::join_horiz(arma::ones(X.n_rows,1), X);
}

// (theta_i, tau_i) ~ BVN((mu_theta_i, mu_tau_i), Sigma_p)：
// 抽 theta 時要用 theta | tau 的條件先驗，抽 tau 時用 tau | theta，
// 否則 Sigma_p 的共變（RT 對能力的資訊）不會進入人員參數的更新。
// k = 0 → theta | tau；k = 1 → tau | theta
inline arma::vec cond_person_mean(const arma::vec& mu_self,
                                  const arma::vec& mu_other,
                                  const arma::vec& other,
                                  const arma::mat& Sigma_p, int k) {
  int o = 1 - k;
  return mu_self + (Sigma_p(k, o) / Sigma_p(o, o)) * (other - mu_other);
}

inline double cond_person_var(const arma::mat& Sigma_p, int k) {
  int o = 1 - k;
  return Sigma_p(k, k) - Sigma_p(k, o) * Sigma_p(k, o) / Sigma_p(o, o);
}

inline void identify_theta_scale(arma::vec& theta,
                                 arma::vec& a,
                                 arma::vec& b) {
  double m = arma::mean(theta);
  theta -= m;
  double s = arma::stddev(theta);
  if (!std::isfinite(s) || s < 1e-8) return;
  theta /= s;
  a *= s;
  b = (b - m) / s;
}

inline void identify_theta_scale_cross(arma::vec& theta,
                                       arma::vec& a,
                                       arma::vec& b,
                                       arma::vec& rho_rt) {
  double m = arma::mean(theta);
  theta -= m;
  double s = arma::stddev(theta);
  if (!std::isfinite(s) || s < 1e-8) return;
  theta /= s;
  a *= s;
  b = (b - m) / s;
  rho_rt *= s;
}

// =============================================================================
// SECTION 8: 主 Gibbs Sampler
// =============================================================================

// ---------------------------------------------------------------------------
// gibbs_rtirt_null：無共變數的 RT-IRT
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
List gibbs_rtirt_null(const arma::mat& Y,
                       const arma::mat& logT,
                       int    n_iter   = 5000,
                       int    n_burnin = 2500,
                       bool   one_pl   = false,
                       double mu_a     = 1.0,
                       double sigma_a  = 1.0,
                       double mu_xi    = -1.0,
                       double sigma_xi = 1e6,
                       double delta_a  = 0.0,
                       double delta_b  = 0.0,
                       int    n_inner  = 2,     // collapsed Gibbs: K=2 最佳 ESS/sec（K=1=standard Gibbs）
                       bool   use_ab_mh = false,
                       double prop_sd_loga = 0.08,
                       double prop_sd_b    = 0.12,
                       bool   adaptive_mh = false,
                       double target_accept = 0.30,
                       int    adapt_every = 25,
                       bool   estimate_rt_slope = false,
                       double mu_slope = 1.0,
                       double sigma_slope = 0.5,
                       bool   hierarchical_rt_slope = true,
                       bool   adaptive_slope_mh = false,
                       double target_accept_slope = 0.30,
                       int    adapt_every_slope = 25,
                       bool   verbose  = true,
                       int    chain_id = 1) {
  int n = Y.n_rows, p = Y.n_cols, n_save = n_iter - n_burnin;
  arma::mat kappa  = Y - 0.5;
  double mu_xi_use = (mu_xi < 0.0) ? arma::mean(arma::mean(logT)) : mu_xi;

  // 初始值：從資料導向的起點，避免第一步爆炸
  arma::vec theta(n, arma::fill::zeros);   // mean=0 已滿足
  arma::vec tau(n,   arma::fill::zeros);
  arma::vec a(p,     arma::fill::ones);    // a=1 (1PL 起點)
  arma::vec b(p,     arma::fill::zeros);
  arma::vec xi      = arma::mean(logT, 0).t();  // 用資料均值初始化
  arma::vec sigma_t(p, arma::fill::value(0.5));
  arma::vec slope_rt(p, arma::fill::ones);
  arma::vec log_slope = arma::log(slope_rt);
  double mu_log_s = std::log(std::max(mu_slope, 1e-6));
  double sigma2_log_s = sigma_slope * sigma_slope;
  double prop_sd_log_s_cur = 0.08;
  arma::mat Sigma_p = arma::eye(2, 2);
  arma::vec mu0(n, arma::fill::zeros);

  arma::mat post_theta(n_save,n), post_tau(n_save,n);
  arma::mat post_a(n_save,p),     post_b(n_save,p);
  arma::mat post_xi(n_save,p),    post_st(n_save,p), post_slope(n_save,p);
  arma::mat post_Sp(n_save,4);    arma::vec post_ll(n_save);
  arma::vec post_acc_ab(n_save, arma::fill::zeros);
  arma::vec post_acc_slope(n_save, arma::fill::zeros);
  arma::mat post_slope_hyp(n_save, 2, arma::fill::zeros); // mu_log_s, sigma2_log_s
  double prop_sd_loga_cur = prop_sd_loga;
  double prop_sd_b_cur    = prop_sd_b;
  int adapt_count = 0;
  int adapt_count_slope = 0;
  int idx = 0;

  if (verbose)
    Rcpp::Rcout << "Chain " << chain_id << " [rtirt_null] "
                << n_iter << " iter, burnin=" << n_burnin << "\n";

  for (int m = 0; m < n_iter; ++m) {
    Rcpp::checkUserInterrupt();
    if (verbose && (m % std::max(1,n_iter/50)==0 || m==n_iter-1)) {
      std::string extra = "";
      if (use_ab_mh) extra = "abMH sd=(" + std::to_string(prop_sd_loga_cur).substr(0,5) + "," +
                              std::to_string(prop_sd_b_cur).substr(0,5) + ")";
      progress_bar_ex(m, n_iter, n_burnin, chain_id, extra);
    }

    // 1. Sigma_p（IW posterior）
    Sigma_p = draw_Sigma_null(theta, tau);

    // 2. PG 輔助
    arma::mat omega = draw_pg_irt(theta, a, b);
    arma::vec mu_theta_c = cond_person_mean(mu0, mu0, tau, Sigma_p, 0);
    double    v_theta_c  = cond_person_var(Sigma_p, 0);

    // 3. IRT：collapsed Gibbs — 每次 outer iter 做 n_inner 步 inner (omega, a, b)
    //    n_inner=1 等同標準 Gibbs；n_inner=3 ESS 提升 ~3x，時間增加 ~60%
    if (!one_pl && use_ab_mh) {
      List ab_mh = draw_ab_mh(Y, theta, a, b, mu_a, sigma_a, 0.0, 1.0,
                              prop_sd_loga_cur, prop_sd_b_cur);
      a = as<arma::vec>(ab_mh["a"]);
      b = as<arma::vec>(ab_mh["b"]);
      omega = draw_pg_irt(theta, a, b);
      theta = draw_theta(kappa, omega, a, b, mu_theta_c, v_theta_c);
      if (adaptive_mh && m < n_burnin && adapt_every > 0 && (m + 1) % adapt_every == 0) {
        double acc = as<double>(ab_mh["accept_rate"]);
        double gamma = 1.0 / std::sqrt(1.0 + adapt_count);
        prop_sd_loga_cur = std::exp(std::log(std::max(prop_sd_loga_cur, 1e-6)) + gamma * (acc - target_accept));
        prop_sd_b_cur    = std::exp(std::log(std::max(prop_sd_b_cur, 1e-6)) + gamma * (acc - target_accept));
        prop_sd_loga_cur = std::min(std::max(prop_sd_loga_cur, 1e-3), 1.0);
        prop_sd_b_cur    = std::min(std::max(prop_sd_b_cur, 1e-3), 2.0);
        adapt_count += 1;
      }
      if (m >= n_burnin) post_acc_ab[idx] = as<double>(ab_mh["accept_rate"]);
    } else if (!one_pl) {
      // Collapsed: (omega, a, b, theta) cycle K times — theta included in inner loop
      List ab = draw_ab_collapsed(kappa, omega, theta, a, b,
                                   mu_theta_c, v_theta_c,
                                   mu_a, sigma_a, 0.0, 1.0, n_inner);
      a     = as<arma::vec>(ab["a"]);
      b     = as<arma::vec>(ab["b"]);
      theta = as<arma::vec>(ab["theta"]);  // theta already centered inside
      omega = as<arma::mat>(ab["omega"]);
    } else {
      b = draw_b(kappa, omega, theta, a, 0.0, 1.0);
      theta = draw_theta(kappa, omega, a, b, mu_theta_c, v_theta_c);
    }

    // 5. RT 項目
    xi      = draw_xi(logT, tau, slope_rt, sigma_t, mu_xi_use, sigma_xi);
    if (estimate_rt_slope) {
      if (hierarchical_rt_slope) {
        List sl = draw_log_slope_mh(logT, xi, tau, sigma_t, log_slope, mu_log_s, sigma2_log_s, prop_sd_log_s_cur);
        log_slope = as<arma::vec>(sl["log_slope"]);
        double acc_sl = as<double>(sl["accept_rate"]);
        if (adaptive_slope_mh && m < n_burnin && adapt_every_slope > 0 && (m + 1) % adapt_every_slope == 0) {
          double gamma = 1.0 / std::sqrt(1.0 + adapt_count_slope);
          prop_sd_log_s_cur = std::exp(std::log(std::max(prop_sd_log_s_cur, 1e-6)) + gamma * (acc_sl - target_accept_slope));
          prop_sd_log_s_cur = std::min(std::max(prop_sd_log_s_cur, 1e-3), 1.0);
          adapt_count_slope += 1;
        }
        mu_log_s = draw_mu_log_slope(log_slope, sigma2_log_s, std::log(std::max(mu_slope, 1e-6)), 4.0);
        sigma2_log_s = draw_sigma2_log_slope(log_slope, mu_log_s, 2.0, 0.5);
        slope_rt = arma::exp(log_slope);
        if (m >= n_burnin) post_acc_slope[idx] = acc_sl;
      } else {
        slope_rt = draw_slope_rt(logT, xi, tau, sigma_t, mu_slope, sigma_slope);
      }
    }
    sigma_t = draw_sigma_t(logT, xi, tau, slope_rt, delta_a, delta_b);

    // 6. 速度
    tau = draw_tau(logT, xi, sigma_t, slope_rt,
                   cond_person_mean(mu0, mu0, theta, Sigma_p, 1),
                   cond_person_var(Sigma_p, 1));

    if (m >= n_burnin) {
      post_theta.row(idx) = theta.t();
      post_tau.row(idx)   = tau.t();
      post_a.row(idx)     = a.t();
      post_b.row(idx)     = b.t();
      post_xi.row(idx)    = xi.t();
      post_st.row(idx)    = sigma_t.t();
      post_slope.row(idx) = slope_rt.t();
      post_slope_hyp(idx, 0) = mu_log_s;
      post_slope_hyp(idx, 1) = sigma2_log_s;
      post_Sp.row(idx)    = arma::vectorise(Sigma_p).t();
      post_ll(idx) = loglik_irt(Y,theta,a,b) + loglik_rt(logT,xi,tau,sigma_t,slope_rt);
      ++idx;
    }
  }
  if (verbose) Rcpp::Rcout << "\n  Done.\n";

  return List::create(
    Named("theta")   = post_theta, Named("tau")     = post_tau,
    Named("a")       = post_a,     Named("b")       = post_b,
    Named("xi")      = post_xi,    Named("sigma_t") = post_st,
    Named("slope_rt")= post_slope,
    Named("Sigma_p") = post_Sp,    Named("loglik")  = post_ll,
    Named("person_ability")   = post_theta,
    Named("person_speed")     = post_tau,
    Named("item_discrimination") = post_a,
    Named("item_difficulty")  = post_b,
    Named("item_time_intensity") = post_xi,
    Named("rt_resid_sd")      = post_st,
    Named("rt_speed_loading") = post_slope,
    Named("person_cov")       = post_Sp,
    Named("log_likelihood")   = post_ll,
    Named("accept_ab") = post_acc_ab,
    Named("accept_slope") = post_acc_slope,
    Named("slope_hyper") = post_slope_hyp,
    Named("final_prop_sd_loga") = prop_sd_loga_cur,
    Named("final_prop_sd_b")    = prop_sd_b_cur,
    Named("final_prop_sd_log_slope") = prop_sd_log_s_cur
  );
}

// ---------------------------------------------------------------------------
// gibbs_rtirt：含結構共變數
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
List gibbs_rtirt(const arma::mat& Y,
                  const arma::mat& logT,
                  const arma::mat& X_cov,
                  int    n_iter   = 5000,
                  int    n_burnin = 2500,
                  bool   fix_int  = true,
                  bool   one_pl   = false,
                  double mu_a     = 1.0,
                  double sigma_a  = 1.0,
                  double mu_xi    = -1.0,
                  double sigma_xi = 1e6,
                  double delta_a  = 0.0,
                  double delta_b  = 0.0,
                  int    n_inner  = 2,
                  bool   use_ab_mh = false,
                  double prop_sd_loga = 0.08,
                  double prop_sd_b    = 0.12,
                  bool   adaptive_mh = false,
                  double target_accept = 0.30,
                  int    adapt_every = 25,
                  bool   estimate_rt_slope = false,
                  double mu_slope = 1.0,
                  double sigma_slope = 0.5,
                  bool   hierarchical_rt_slope = true,
                  bool   adaptive_slope_mh = false,
                  double target_accept_slope = 0.30,
                  int    adapt_every_slope = 25,
                  bool   verbose  = true,
                  int    chain_id = 1) {
  int n = Y.n_rows, p = Y.n_cols;
  arma::mat X = make_design(X_cov);
  int q = X.n_cols, n_save = n_iter - n_burnin;
  arma::mat kappa  = Y - 0.5;
  double mu_xi_use = (mu_xi < 0.0) ? arma::mean(arma::mean(logT)) : mu_xi;

  arma::vec theta(n, arma::fill::zeros);
  arma::vec tau(n,   arma::fill::zeros);
  arma::vec a(p,     arma::fill::ones);
  arma::vec b(p,     arma::fill::zeros);
  arma::vec xi      = arma::mean(logT, 0).t();
  arma::vec sigma_t(p, arma::fill::value(0.5));
  arma::vec slope_rt(p, arma::fill::ones);
  arma::vec log_slope = arma::log(slope_rt);
  double mu_log_s = std::log(std::max(mu_slope, 1e-6));
  double sigma2_log_s = sigma_slope * sigma_slope;
  arma::mat beta(q, 2, arma::fill::zeros);
  arma::mat Sigma_p = arma::eye(2, 2);

  arma::mat post_theta(n_save,n), post_tau(n_save,n);
  arma::mat post_a(n_save,p),     post_b(n_save,p);
  arma::mat post_xi(n_save,p),    post_st(n_save,p), post_slope(n_save,p);
  arma::mat post_beta(n_save,q*2), post_Sp(n_save,4);
  arma::vec post_ll(n_save), post_acc_ab(n_save, arma::fill::zeros);
  arma::vec post_acc_slope(n_save, arma::fill::zeros);
  arma::mat post_slope_hyp(n_save, 2, arma::fill::zeros);
  double prop_sd_loga_cur = prop_sd_loga;
  double prop_sd_b_cur    = prop_sd_b;
  double prop_sd_log_s_cur = 0.08;
  int adapt_count = 0;
  int adapt_count_slope = 0;
  int idx = 0;

  if (verbose)
    Rcpp::Rcout << "Chain " << chain_id << " [rtirt] "
                << n_iter << " iter, burnin=" << n_burnin << "\n";

  for (int m = 0; m < n_iter; ++m) {
    Rcpp::checkUserInterrupt();
    if (verbose && (m % std::max(1,n_iter/50)==0 || m==n_iter-1)) {
      std::string extra = "";
      if (use_ab_mh) extra = "abMH sd=(" + std::to_string(prop_sd_loga_cur).substr(0,5) + "," +
                              std::to_string(prop_sd_b_cur).substr(0,5) + ")";
      progress_bar_ex(m, n_iter, n_burnin, chain_id, extra);
    }

    beta    = draw_beta(X, theta, tau, Sigma_p);
    if (fix_int) beta.row(0).zeros();
    Sigma_p = draw_Sigma_struct(X, theta, tau, beta);

    arma::mat omega = draw_pg_irt(theta, a, b);
    arma::vec mu_theta = X * beta.col(0);
    arma::vec mu_tau   = X * beta.col(1);
    arma::vec mu_theta_c = cond_person_mean(mu_theta, mu_tau, tau, Sigma_p, 0);
    double    v_theta_c  = cond_person_var(Sigma_p, 0);
    if (!one_pl && use_ab_mh) {
      List ab_mh = draw_ab_mh(Y, theta, a, b, mu_a, sigma_a, 0.0, 1.0,
                              prop_sd_loga_cur, prop_sd_b_cur);
      a = as<arma::vec>(ab_mh["a"]);
      b = as<arma::vec>(ab_mh["b"]);
      omega = draw_pg_irt(theta, a, b);
      theta = draw_theta(kappa, omega, a, b, mu_theta_c, v_theta_c);
      if (adaptive_mh && m < n_burnin && adapt_every > 0 && (m + 1) % adapt_every == 0) {
        double acc = as<double>(ab_mh["accept_rate"]);
        double gamma = 1.0 / std::sqrt(1.0 + adapt_count);
        prop_sd_loga_cur = std::exp(std::log(std::max(prop_sd_loga_cur, 1e-6)) + gamma * (acc - target_accept));
        prop_sd_b_cur    = std::exp(std::log(std::max(prop_sd_b_cur, 1e-6)) + gamma * (acc - target_accept));
        prop_sd_loga_cur = std::min(std::max(prop_sd_loga_cur, 1e-3), 1.0);
        prop_sd_b_cur    = std::min(std::max(prop_sd_b_cur, 1e-3), 2.0);
        adapt_count += 1;
      }
      if (m >= n_burnin) post_acc_ab[idx] = as<double>(ab_mh["accept_rate"]);
    } else if (!one_pl) {
      List ab = draw_ab_collapsed(kappa, omega, theta, a, b,
                                   mu_theta_c, v_theta_c,
                                   mu_a, sigma_a, 0.0, 1.0, n_inner);
      a     = as<arma::vec>(ab["a"]);
      b     = as<arma::vec>(ab["b"]);
      theta = as<arma::vec>(ab["theta"]);
      omega = as<arma::mat>(ab["omega"]);
    } else {
      b = draw_b(kappa, omega, theta, a, 0.0, 1.0);
      theta = draw_theta(kappa, omega, a, b, mu_theta_c, v_theta_c);
    }

    xi      = draw_xi(logT, tau, slope_rt, sigma_t, mu_xi_use, sigma_xi);
    if (estimate_rt_slope) {
      if (hierarchical_rt_slope) {
        List sl = draw_log_slope_mh(logT, xi, tau, sigma_t, log_slope, mu_log_s, sigma2_log_s, prop_sd_log_s_cur);
        log_slope = as<arma::vec>(sl["log_slope"]);
        double acc_sl = as<double>(sl["accept_rate"]);
        if (adaptive_slope_mh && m < n_burnin && adapt_every_slope > 0 && (m + 1) % adapt_every_slope == 0) {
          double gamma = 1.0 / std::sqrt(1.0 + adapt_count_slope);
          prop_sd_log_s_cur = std::exp(std::log(std::max(prop_sd_log_s_cur, 1e-6)) + gamma * (acc_sl - target_accept_slope));
          prop_sd_log_s_cur = std::min(std::max(prop_sd_log_s_cur, 1e-3), 1.0);
          adapt_count_slope += 1;
        }
        mu_log_s = draw_mu_log_slope(log_slope, sigma2_log_s, std::log(std::max(mu_slope, 1e-6)), 4.0);
        sigma2_log_s = draw_sigma2_log_slope(log_slope, mu_log_s, 2.0, 0.5);
        slope_rt = arma::exp(log_slope);
        if (m >= n_burnin) post_acc_slope[idx] = acc_sl;
      } else {
        slope_rt = draw_slope_rt(logT, xi, tau, sigma_t, mu_slope, sigma_slope);
      }
    }
    sigma_t = draw_sigma_t(logT, xi, tau, slope_rt, delta_a, delta_b);

    tau = draw_tau(logT, xi, sigma_t, slope_rt,
                   cond_person_mean(mu_tau, mu_theta, theta, Sigma_p, 1),
                   cond_person_var(Sigma_p, 1));

    if (m >= n_burnin) {
      post_theta.row(idx) = theta.t();
      post_tau.row(idx)   = tau.t();
      post_a.row(idx)     = a.t();
      post_b.row(idx)     = b.t();
      post_xi.row(idx)    = xi.t();
      post_st.row(idx)    = sigma_t.t();
      post_slope.row(idx) = slope_rt.t();
      post_slope_hyp(idx, 0) = mu_log_s;
      post_slope_hyp(idx, 1) = sigma2_log_s;
      post_beta.row(idx)  = arma::vectorise(beta).t();
      post_Sp.row(idx)    = arma::vectorise(Sigma_p).t();
      post_ll(idx) = loglik_irt(Y,theta,a,b) + loglik_rt(logT,xi,tau,sigma_t,slope_rt);
      ++idx;
    }
  }
  if (verbose) Rcpp::Rcout << "\n  Done.\n";

  return List::create(
    Named("theta")   = post_theta, Named("tau")     = post_tau,
    Named("a")       = post_a,     Named("b")       = post_b,
    Named("xi")      = post_xi,    Named("sigma_t") = post_st,
    Named("slope_rt")= post_slope,
    Named("beta")    = post_beta,  Named("Sigma_p") = post_Sp,
    Named("loglik")  = post_ll,    Named("accept_ab") = post_acc_ab,
    Named("person_ability")   = post_theta,
    Named("person_speed")     = post_tau,
    Named("item_discrimination") = post_a,
    Named("item_difficulty")  = post_b,
    Named("item_time_intensity") = post_xi,
    Named("rt_resid_sd")      = post_st,
    Named("rt_speed_loading") = post_slope,
    Named("person_regression") = post_beta,
    Named("person_cov")       = post_Sp,
    Named("log_likelihood")   = post_ll,
    Named("accept_slope") = post_acc_slope,
    Named("slope_hyper") = post_slope_hyp,
    Named("final_prop_sd_loga") = prop_sd_loga_cur,
    Named("final_prop_sd_b")    = prop_sd_b_cur,
    Named("final_prop_sd_log_slope") = prop_sd_log_s_cur
  );
}

// ---------------------------------------------------------------------------
// gibbs_ml_irt：僅 IRT，theta ~ regression
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
List gibbs_ml_irt(const arma::mat& Y,
                   const arma::mat& X_cov,
                   int    n_iter   = 5000,
                   int    n_burnin = 2500,
                   bool   fix_int  = true,
                   bool   one_pl   = false,
                   double mu_a     = 1.0,
                   double sigma_a  = 1.0,
                   int    n_inner  = 2,
                   bool   use_ab_mh = false,
                   double prop_sd_loga = 0.08,
                   double prop_sd_b    = 0.12,
                   bool   adaptive_mh = false,
                   double target_accept = 0.30,
                   int    adapt_every = 25,
                   bool   verbose  = true,
                   int    chain_id = 1) {
  int n = Y.n_rows, p = Y.n_cols;
  arma::mat X     = make_design(X_cov);
  int q = X.n_cols, n_save = n_iter - n_burnin;
  arma::mat kappa = Y - 0.5;

  arma::vec theta(n, arma::fill::zeros);
  arma::vec a(p,     arma::fill::ones);
  arma::vec b(p,     arma::fill::zeros);
  arma::vec beta_ra(q, arma::fill::zeros);
  double prop_sd_loga_cur = prop_sd_loga;
  double prop_sd_b_cur    = prop_sd_b;
  int adapt_count = 0;

  arma::mat post_theta(n_save,n), post_a(n_save,p), post_b(n_save,p);
  arma::mat post_beta(n_save,q);  arma::vec post_ll(n_save), post_acc_ab(n_save, arma::fill::zeros);
  int idx = 0;

  if (verbose)
    Rcpp::Rcout << "Chain " << chain_id << " [ml_irt] "
                << n_iter << " iter, burnin=" << n_burnin << "\n";

  for (int m = 0; m < n_iter; ++m) {
    Rcpp::checkUserInterrupt();
    if (verbose && (m % std::max(1,n_iter/50)==0 || m==n_iter-1))
      progress_bar(m, n_iter, n_burnin, chain_id);

    arma::mat XtX = X.t()*X + 1e-8*arma::eye(q,q);
    beta_ra = arma::solve(XtX, X.t()*theta);
    if (fix_int) beta_ra[0] = 0.0;

    arma::mat omega = draw_pg_irt(theta, a, b);
    arma::vec mu_theta = X * beta_ra;
    if (!one_pl && use_ab_mh) {
      List ab_mh = draw_ab_mh(Y, theta, a, b, mu_a, sigma_a, 0.0, 1.0,
                              prop_sd_loga_cur, prop_sd_b_cur);
      a = as<arma::vec>(ab_mh["a"]);
      b = as<arma::vec>(ab_mh["b"]);
      omega = draw_pg_irt(theta, a, b);
      theta = draw_theta(kappa, omega, a, b, mu_theta, 1.0);
      if (adaptive_mh && m < n_burnin && adapt_every > 0 && (m + 1) % adapt_every == 0) {
        double acc = as<double>(ab_mh["accept_rate"]);
        double gamma = 1.0 / std::sqrt(1.0 + adapt_count);
        prop_sd_loga_cur = std::exp(std::log(std::max(prop_sd_loga_cur, 1e-6)) + gamma * (acc - target_accept));
        prop_sd_b_cur    = std::exp(std::log(std::max(prop_sd_b_cur, 1e-6)) + gamma * (acc - target_accept));
        prop_sd_loga_cur = std::min(std::max(prop_sd_loga_cur, 1e-3), 1.0);
        prop_sd_b_cur    = std::min(std::max(prop_sd_b_cur, 1e-3), 2.0);
        adapt_count += 1;
      }
      if (m >= n_burnin) post_acc_ab(idx) = as<double>(ab_mh["accept_rate"]);
    } else if (!one_pl) {
      List ab = draw_ab_collapsed(kappa, omega, theta, a, b,
                                   mu_theta, 1.0,
                                   mu_a, sigma_a, 0.0, 1.0, n_inner);
      a     = as<arma::vec>(ab["a"]);
      b     = as<arma::vec>(ab["b"]);
      theta = as<arma::vec>(ab["theta"]);
      omega = as<arma::mat>(ab["omega"]);
    } else {
      b = draw_b(kappa, omega, theta, a, 0.0, 1.0);
      theta = draw_theta(kappa, omega, a, b, mu_theta, 1.0);
    }

    if (m >= n_burnin) {
      post_theta.row(idx) = theta.t();
      post_a.row(idx)     = a.t();
      post_b.row(idx)     = b.t();
      post_beta.row(idx)  = beta_ra.t();
      post_ll(idx) = loglik_irt(Y,theta,a,b);
      ++idx;
    }
  }
  if (verbose) Rcpp::Rcout << "\n  Done.\n";

  return List::create(
    Named("theta") = post_theta, Named("a") = post_a,
    Named("b")     = post_b,     Named("beta") = post_beta,
    Named("loglik")= post_ll,    Named("accept_ab") = post_acc_ab,
    Named("person_ability") = post_theta,
    Named("item_discrimination") = post_a,
    Named("item_difficulty") = post_b,
    Named("person_regression") = post_beta,
    Named("log_likelihood") = post_ll,
    Named("final_prop_sd_loga") = prop_sd_loga_cur,
    Named("final_prop_sd_b")    = prop_sd_b_cur
  );
}

// ---------------------------------------------------------------------------
// gibbs_rtirt_quantile：RT 部分採 quantile augmentation；RA/theta 保持原模型
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
List gibbs_rtirt_quantile(const arma::mat& Y,
                          const arma::mat& logT,
                          const arma::mat& X_cov,
                          int    n_iter   = 5000,
                          int    n_burnin = 2500,
                          double q_rt     = 0.5,
                          bool   fix_int  = true,
                          bool   one_pl   = false,
                          double mu_a     = 1.0,
                          double sigma_a  = 1.0,
                          double mu_xi    = -1.0,
                          double sigma_xi = 1e6,
                          double delta_a  = 0.0,
                          double delta_b  = 0.0,
                          int    n_inner  = 2,
                          bool   use_ab_mh = false,
                          double prop_sd_loga = 0.08,
                          double prop_sd_b    = 0.12,
                          bool   adaptive_mh = false,
                          double target_accept = 0.30,
                          int    adapt_every = 25,
                          bool   standardize_theta = true,
                          bool   estimate_rt_slope = false,
                          bool   hierarchical_rt_slope = true,
                          bool   adaptive_slope_mh = false,
                          double target_accept_slope = 0.30,
                          int    adapt_every_slope = 25,
                          double mu_slope = 1.0,
                          double sigma_slope = 0.5,
                          bool   verbose  = true,
                          int    chain_id = 1) {
  int n = Y.n_rows, p = Y.n_cols;
  arma::mat X = make_design(X_cov);
  int q = X.n_cols, n_save = n_iter - n_burnin;
  arma::mat kappa  = Y - 0.5;
  double mu_xi_use = (mu_xi < 0.0) ? arma::mean(arma::mean(logT)) : mu_xi;

  arma::vec theta(n, arma::fill::zeros), tau(n, arma::fill::zeros);
  arma::vec a(p, arma::fill::ones), b(p, arma::fill::zeros);
  arma::vec xi = arma::mean(logT, 0).t();
  arma::vec sigma_t(p, arma::fill::value(0.5));
  arma::vec slope_rt(p, arma::fill::ones);
  arma::vec log_slope = arma::log(slope_rt);
  double mu_log_s = std::log(std::max(mu_slope, 1e-6));
  double sigma2_log_s = sigma_slope * sigma_slope;
  arma::mat nu_rt(n, p, arma::fill::ones);
  arma::mat beta(q, 2, arma::fill::zeros);
  arma::mat Sigma_p = arma::eye(2, 2);

  arma::mat post_theta(n_save,n), post_tau(n_save,n), post_a(n_save,p), post_b(n_save,p);
  arma::mat post_xi(n_save,p), post_st(n_save,p), post_slope(n_save,p), post_beta(n_save,q*2), post_Sp(n_save,4);
  arma::vec post_ll(n_save), post_acc_ab(n_save, arma::fill::zeros), post_acc_slope(n_save, arma::fill::zeros);
  arma::mat post_slope_hyp(n_save, 2, arma::fill::zeros);
  double prop_sd_loga_cur = prop_sd_loga, prop_sd_b_cur = prop_sd_b;
  double prop_sd_log_s_cur = 0.08;
  int adapt_count = 0, adapt_count_slope = 0, idx = 0;

  if (verbose) Rcpp::Rcout << "Chain " << chain_id << " [rtirt_quantile] "
                           << n_iter << " iter, burnin=" << n_burnin
                           << ", q_rt=" << q_rt << "\n";

  for (int m = 0; m < n_iter; ++m) {
    Rcpp::checkUserInterrupt();
    if (verbose && (m % std::max(1,n_iter/50)==0 || m==n_iter-1)) {
      std::string extra = "";
      if (use_ab_mh) extra = "abMH";
      if (estimate_rt_slope && hierarchical_rt_slope) {
        extra += " sMH";
      }
      progress_bar_ex(m, n_iter, n_burnin, chain_id, extra);
    }

    beta    = draw_beta(X, theta, tau, Sigma_p);
    if (fix_int) beta.row(0).zeros();
    Sigma_p = draw_Sigma_struct(X, theta, tau, beta);

    arma::mat omega = draw_pg_irt(theta, a, b);
    arma::vec mu_theta = X * beta.col(0);
    arma::vec mu_tau   = X * beta.col(1);
    arma::vec mu_theta_c = cond_person_mean(mu_theta, mu_tau, tau, Sigma_p, 0);
    double    v_theta_c  = cond_person_var(Sigma_p, 0);
    if (!one_pl && use_ab_mh) {
      List ab_mh = draw_ab_mh(Y, theta, a, b, mu_a, sigma_a, 0.0, 1.0,
                              prop_sd_loga_cur, prop_sd_b_cur);
      a = as<arma::vec>(ab_mh["a"]); b = as<arma::vec>(ab_mh["b"]);
      omega = draw_pg_irt(theta, a, b);
      theta = draw_theta(kappa, omega, a, b, mu_theta_c, v_theta_c);
      if (adaptive_mh && m < n_burnin && adapt_every > 0 && (m + 1) % adapt_every == 0) {
        double acc = as<double>(ab_mh["accept_rate"]);
        double gamma = 1.0 / std::sqrt(1.0 + adapt_count);
        prop_sd_loga_cur = std::exp(std::log(std::max(prop_sd_loga_cur, 1e-6)) + gamma * (acc - target_accept));
        prop_sd_b_cur    = std::exp(std::log(std::max(prop_sd_b_cur, 1e-6)) + gamma * (acc - target_accept));
        prop_sd_loga_cur = std::min(std::max(prop_sd_loga_cur, 1e-3), 1.0);
        prop_sd_b_cur    = std::min(std::max(prop_sd_b_cur, 1e-3), 2.0);
        adapt_count += 1;
      }
      if (m >= n_burnin) post_acc_ab[idx] = as<double>(ab_mh["accept_rate"]);
    } else if (!one_pl) {
      List ab = draw_ab_collapsed(kappa, omega, theta, a, b,
                                  mu_theta_c, v_theta_c,
                                  mu_a, sigma_a, 0.0, 1.0, n_inner);
      a = as<arma::vec>(ab["a"]); b = as<arma::vec>(ab["b"]);
      theta = as<arma::vec>(ab["theta"]);
    } else {
      b = draw_b(kappa, omega, theta, a, 0.0, 1.0);
      theta = draw_theta(kappa, omega, a, b, mu_theta_c, v_theta_c);
    }
    if (standardize_theta && !one_pl) identify_theta_scale(theta, a, b);

    nu_rt = draw_qr_weights_rt(logT, xi, tau, slope_rt, sigma_t, q_rt);
    xi    = draw_xi_qr(logT, tau, slope_rt, sigma_t, nu_rt, q_rt, mu_xi_use, sigma_xi);
    if (estimate_rt_slope) {
      if (hierarchical_rt_slope) {
        List sl = draw_log_slope_mh(logT, xi, tau, sigma_t, log_slope,
                                    mu_log_s, sigma2_log_s, prop_sd_log_s_cur);
        log_slope = as<arma::vec>(sl["log_slope"]);
        double acc_sl = as<double>(sl["accept_rate"]);
        if (adaptive_slope_mh && m < n_burnin && adapt_every_slope > 0 && (m + 1) % adapt_every_slope == 0) {
          double gamma = 1.0 / std::sqrt(1.0 + adapt_count_slope);
          prop_sd_log_s_cur = std::exp(std::log(std::max(prop_sd_log_s_cur, 1e-6)) + gamma * (acc_sl - target_accept_slope));
          prop_sd_log_s_cur = std::min(std::max(prop_sd_log_s_cur, 1e-3), 1.0);
          adapt_count_slope += 1;
        }
        mu_log_s = draw_mu_log_slope(log_slope, sigma2_log_s, std::log(std::max(mu_slope, 1e-6)), 4.0);
        sigma2_log_s = draw_sigma2_log_slope(log_slope, mu_log_s, 2.0, 0.5);
        slope_rt = arma::exp(log_slope);
        if (m >= n_burnin) post_acc_slope[idx] = acc_sl;
      } else {
        slope_rt = draw_slope_rt(logT, xi, tau, sigma_t, mu_slope, sigma_slope);
      }
    }
    sigma_t = draw_sigma_t(logT, xi, tau, slope_rt, delta_a, delta_b);
    tau = draw_tau_qr(logT, xi, sigma_t, slope_rt, nu_rt,
                      cond_person_mean(mu_tau, mu_theta, theta, Sigma_p, 1),
                      cond_person_var(Sigma_p, 1), q_rt);

    if (m >= n_burnin) {
      post_theta.row(idx) = theta.t(); post_tau.row(idx) = tau.t();
      post_a.row(idx) = a.t(); post_b.row(idx) = b.t();
      post_xi.row(idx) = xi.t(); post_st.row(idx) = sigma_t.t(); post_slope.row(idx) = slope_rt.t();
      post_slope_hyp(idx, 0) = mu_log_s;
      post_slope_hyp(idx, 1) = sigma2_log_s;
      post_beta.row(idx) = arma::vectorise(beta).t(); post_Sp.row(idx) = arma::vectorise(Sigma_p).t();
      post_ll(idx) = loglik_irt(Y,theta,a,b) + loglik_rt(logT,xi,tau,sigma_t,slope_rt);
      ++idx;
    }
  }
  if (verbose) Rcpp::Rcout << "\n  Done.\n";

  return List::create(
    Named("theta") = post_theta, Named("tau") = post_tau,
    Named("a") = post_a, Named("b") = post_b,
    Named("xi") = post_xi, Named("sigma_t") = post_st, Named("slope_rt") = post_slope,
    Named("beta") = post_beta, Named("Sigma_p") = post_Sp, Named("loglik") = post_ll,
    Named("person_ability") = post_theta,
    Named("person_speed") = post_tau,
    Named("item_discrimination") = post_a,
    Named("item_difficulty") = post_b,
    Named("item_time_intensity") = post_xi,
    Named("rt_resid_sd") = post_st,
    Named("rt_speed_loading") = post_slope,
    Named("person_regression") = post_beta,
    Named("person_cov") = post_Sp,
    Named("log_likelihood") = post_ll,
    Named("accept_ab") = post_acc_ab,
    Named("accept_slope") = post_acc_slope,
    Named("slope_hyper") = post_slope_hyp,
    Named("final_prop_sd_loga") = prop_sd_loga_cur, Named("final_prop_sd_b") = prop_sd_b_cur,
    Named("final_prop_sd_log_slope") = prop_sd_log_s_cur
  );
}

// [[Rcpp::export]]
List gibbs_rtirt_cross_quantile(const arma::mat& Y,
                                const arma::mat& logT,
                                const arma::mat& X_cov,
                                int    n_iter   = 5000,
                                int    n_burnin = 2500,
                                double q_rt     = 0.5,
                                bool   fix_int  = true,
                                bool   one_pl   = false,
                                double mu_a     = 1.0,
                                double sigma_a  = 1.0,
                                double mu_xi    = -1.0,
                                double sigma_xi = 1e6,
                                double delta_a  = 0.0,
                                double delta_b  = 0.0,
                                int    n_inner  = 2,
                                bool   use_ab_mh = false,
                                double prop_sd_loga = 0.08,
                                double prop_sd_b    = 0.12,
                                bool   adaptive_mh = false,
                                double target_accept = 0.30,
                                int    adapt_every = 25,
                                bool   standardize_theta = true,
                                bool   estimate_rt_slope = false,
                                bool   hierarchical_rt_slope = false,
                                bool   adaptive_slope_mh = false,
                                double target_accept_slope = 0.30,
                                int    adapt_every_slope = 25,
                                double mu_slope = 1.0,
                                double sigma_slope = 0.5,
                                double mu_rho = 0.0,
                                double sigma_rho = 1.0,
                                bool   verbose  = true,
                                int    chain_id = 1) {
  int n = Y.n_rows, p = Y.n_cols;
  arma::mat X = make_design(X_cov);
  int q = X.n_cols, n_save = n_iter - n_burnin;
  arma::mat kappa  = Y - 0.5;
  double mu_xi_use = (mu_xi < 0.0) ? arma::mean(arma::mean(logT)) : mu_xi;

  arma::vec theta(n, arma::fill::zeros), tau(n, arma::fill::zeros);
  arma::vec a(p, arma::fill::ones), b(p, arma::fill::zeros);
  arma::vec xi = arma::mean(logT, 0).t();
  arma::vec sigma_t(p, arma::fill::value(0.5));
  arma::vec slope_rt(p, arma::fill::ones), rho_rt(p, arma::fill::zeros);
  arma::vec log_slope = arma::log(slope_rt);
  double mu_log_s = std::log(std::max(mu_slope, 1e-6));
  double sigma2_log_s = sigma_slope * sigma_slope;
  arma::mat nu_rt(n, p, arma::fill::ones);
  arma::mat beta(q, 2, arma::fill::zeros);
  arma::mat Sigma_p = arma::eye(2, 2);

  arma::mat post_theta(n_save,n), post_tau(n_save,n), post_a(n_save,p), post_b(n_save,p);
  arma::mat post_xi(n_save,p), post_st(n_save,p), post_slope(n_save,p), post_rho(n_save,p), post_beta(n_save,q*2), post_Sp(n_save,4);
  arma::vec post_ll(n_save), post_acc_ab(n_save, arma::fill::zeros), post_acc_slope(n_save, arma::fill::zeros);
  arma::mat post_slope_hyp(n_save, 2, arma::fill::zeros);
  double prop_sd_loga_cur = prop_sd_loga, prop_sd_b_cur = prop_sd_b, prop_sd_log_s_cur = 0.08;
  int adapt_count = 0, adapt_count_slope = 0, idx = 0;

  if (verbose) Rcpp::Rcout << "Chain " << chain_id << " [rtirt_cross_quantile] "
                           << n_iter << " iter, burnin=" << n_burnin
                           << ", q_rt=" << q_rt << "\n";

  for (int m = 0; m < n_iter; ++m) {
    Rcpp::checkUserInterrupt();
    if (verbose && (m % std::max(1,n_iter/50)==0 || m==n_iter-1)) {
      std::string extra = "";
      if (use_ab_mh) extra = "abMH";
      if (estimate_rt_slope && hierarchical_rt_slope) extra += " sMH";
      progress_bar_ex(m, n_iter, n_burnin, chain_id, extra);
    }

    beta    = draw_beta(X, theta, tau, Sigma_p);
    if (fix_int) beta.row(0).zeros();
    Sigma_p = draw_Sigma_struct(X, theta, tau, beta);

    arma::mat omega = draw_pg_irt(theta, a, b);
    arma::vec mu_theta = X * beta.col(0);
    arma::vec mu_tau   = X * beta.col(1);
    arma::vec mu_theta_c = cond_person_mean(mu_theta, mu_tau, tau, Sigma_p, 0);
    double    v_theta_c  = cond_person_var(Sigma_p, 0);
    if (!one_pl && use_ab_mh) {
      List ab_mh = draw_ab_mh(Y, theta, a, b, mu_a, sigma_a, 0.0, 1.0,
                              prop_sd_loga_cur, prop_sd_b_cur);
      a = as<arma::vec>(ab_mh["a"]); b = as<arma::vec>(ab_mh["b"]);
      omega = draw_pg_irt(theta, a, b);
      theta = draw_theta(kappa, omega, a, b, mu_theta_c, v_theta_c);
      if (adaptive_mh && m < n_burnin && adapt_every > 0 && (m + 1) % adapt_every == 0) {
        double acc = as<double>(ab_mh["accept_rate"]);
        double gamma = 1.0 / std::sqrt(1.0 + adapt_count);
        prop_sd_loga_cur = std::exp(std::log(std::max(prop_sd_loga_cur, 1e-6)) + gamma * (acc - target_accept));
        prop_sd_b_cur    = std::exp(std::log(std::max(prop_sd_b_cur, 1e-6)) + gamma * (acc - target_accept));
        prop_sd_loga_cur = std::min(std::max(prop_sd_loga_cur, 1e-3), 1.0);
        prop_sd_b_cur    = std::min(std::max(prop_sd_b_cur, 1e-3), 2.0);
        adapt_count += 1;
      }
      if (m >= n_burnin) post_acc_ab[idx] = as<double>(ab_mh["accept_rate"]);
    } else if (!one_pl) {
      List ab = draw_ab_collapsed(kappa, omega, theta, a, b,
                                  mu_theta_c, v_theta_c,
                                  mu_a, sigma_a, 0.0, 1.0, n_inner);
      a = as<arma::vec>(ab["a"]); b = as<arma::vec>(ab["b"]);
      theta = as<arma::vec>(ab["theta"]);
    } else {
      b = draw_b(kappa, omega, theta, a, 0.0, 1.0);
      theta = draw_theta(kappa, omega, a, b, mu_theta_c, v_theta_c);
    }
    if (standardize_theta && !one_pl) identify_theta_scale_cross(theta, a, b, rho_rt);

    nu_rt = draw_qr_weights_rt_cross(logT, xi, tau, theta, slope_rt, rho_rt, sigma_t, q_rt);
    xi    = draw_xi_qr_cross(logT, tau, theta, slope_rt, rho_rt, sigma_t, nu_rt, q_rt, mu_xi_use, sigma_xi);
    if (estimate_rt_slope) {
      if (hierarchical_rt_slope) {
        arma::mat logT_adj = logT + theta * rho_rt.t();
        List sl = draw_log_slope_mh(logT_adj, xi, tau, sigma_t, log_slope,
                                    mu_log_s, sigma2_log_s, prop_sd_log_s_cur);
        log_slope = as<arma::vec>(sl["log_slope"]);
        double acc_sl = as<double>(sl["accept_rate"]);
        if (adaptive_slope_mh && m < n_burnin && adapt_every_slope > 0 && (m + 1) % adapt_every_slope == 0) {
          double gamma = 1.0 / std::sqrt(1.0 + adapt_count_slope);
          prop_sd_log_s_cur = std::exp(std::log(std::max(prop_sd_log_s_cur, 1e-6)) + gamma * (acc_sl - target_accept_slope));
          prop_sd_log_s_cur = std::min(std::max(prop_sd_log_s_cur, 1e-3), 1.0);
          adapt_count_slope += 1;
        }
        mu_log_s = draw_mu_log_slope(log_slope, sigma2_log_s, std::log(std::max(mu_slope, 1e-6)), 4.0);
        sigma2_log_s = draw_sigma2_log_slope(log_slope, mu_log_s, 2.0, 0.5);
        slope_rt = arma::exp(log_slope);
        if (m >= n_burnin) post_acc_slope[idx] = acc_sl;
      } else {
        arma::mat logT_adj = logT + theta * rho_rt.t();
        slope_rt = draw_slope_rt(logT_adj, xi, tau, sigma_t, mu_slope, sigma_slope);
      }
    }
    rho_rt  = draw_rho_qr_cross(logT, xi, tau, theta, slope_rt, sigma_t, nu_rt, q_rt, mu_rho, sigma_rho);
    sigma_t = draw_sigma_t(logT, xi, tau, slope_rt, delta_a, delta_b);
    tau = draw_tau_qr_cross(logT, xi, theta, sigma_t, slope_rt, rho_rt, nu_rt,
                            cond_person_mean(mu_tau, mu_theta, theta, Sigma_p, 1),
                            cond_person_var(Sigma_p, 1), q_rt);

    if (m >= n_burnin) {
      post_theta.row(idx) = theta.t(); post_tau.row(idx) = tau.t();
      post_a.row(idx) = a.t(); post_b.row(idx) = b.t();
      post_xi.row(idx) = xi.t(); post_st.row(idx) = sigma_t.t(); post_slope.row(idx) = slope_rt.t(); post_rho.row(idx) = rho_rt.t();
      post_slope_hyp(idx, 0) = mu_log_s; post_slope_hyp(idx, 1) = sigma2_log_s;
      post_beta.row(idx) = arma::vectorise(beta).t(); post_Sp.row(idx) = arma::vectorise(Sigma_p).t();
      post_ll(idx) = loglik_irt(Y,theta,a,b) + loglik_rt(logT,xi,tau,sigma_t,slope_rt);
      ++idx;
    }
  }
  if (verbose) Rcpp::Rcout << "\n  Done.\n";

  return List::create(
    Named("theta") = post_theta, Named("tau") = post_tau,
    Named("a") = post_a, Named("b") = post_b,
    Named("xi") = post_xi, Named("sigma_t") = post_st, Named("slope_rt") = post_slope, Named("rho_rt") = post_rho,
    Named("beta") = post_beta, Named("Sigma_p") = post_Sp, Named("loglik") = post_ll,
    Named("person_ability") = post_theta, Named("person_speed") = post_tau,
    Named("item_discrimination") = post_a, Named("item_difficulty") = post_b,
    Named("item_time_intensity") = post_xi, Named("rt_resid_sd") = post_st,
    Named("rt_speed_loading") = post_slope, Named("rt_cross_loading") = post_rho,
    Named("person_regression") = post_beta, Named("person_cov") = post_Sp,
    Named("log_likelihood") = post_ll,
    Named("accept_ab") = post_acc_ab, Named("accept_slope") = post_acc_slope, Named("slope_hyper") = post_slope_hyp,
    Named("final_prop_sd_loga") = prop_sd_loga_cur, Named("final_prop_sd_b") = prop_sd_b_cur,
    Named("final_prop_sd_log_slope") = prop_sd_log_s_cur
  );
}

// =============================================================================
// SECTION 9: DIC
// =============================================================================
// [[Rcpp::export]]
List compute_dic(double ll_at_mean, const arma::vec& ll_chain) {
  double D_hat = -2.0 * ll_at_mean;
  double D_bar = -2.0 * arma::mean(ll_chain);
  double pD    = D_bar - D_hat;
  return List::create(Named("pD")=pD, Named("DIC")=D_bar+pD);
}

// =============================================================================
// SECTION 10: Descriptive-name wrappers (backward compatible)
// =============================================================================
// [[Rcpp::export]]
arma::vec draw_item_time_intensity(const arma::mat& log_rt,
                                   const arma::vec& person_speed,
                                   const arma::vec& rt_speed_loading,
                                   const arma::vec& rt_resid_sd,
                                   double prior_mean = 4.0,
                                   double prior_sd = 1e6) {
  return draw_xi(log_rt, person_speed, rt_speed_loading, rt_resid_sd, prior_mean, prior_sd);
}

// [[Rcpp::export]]
arma::vec draw_item_time_resid_sd(const arma::mat& log_rt,
                                  const arma::vec& item_time_intensity,
                                  const arma::vec& person_speed,
                                  const arma::vec& rt_speed_loading,
                                  double prior_shape = 0.0,
                                  double prior_rate = 0.0) {
  return draw_sigma_t(log_rt, item_time_intensity, person_speed, rt_speed_loading, prior_shape, prior_rate);
}

// [[Rcpp::export]]
arma::vec draw_person_speed(const arma::mat& log_rt,
                            const arma::vec& item_time_intensity,
                            const arma::vec& rt_resid_sd,
                            const arma::vec& rt_speed_loading,
                            const arma::vec& person_speed_mean,
                            double person_speed_var = 1.0) {
  return draw_tau(log_rt, item_time_intensity, rt_resid_sd, rt_speed_loading, person_speed_mean, person_speed_var);
}

// [[Rcpp::export]]
double loglik_response_time(const arma::mat& log_rt,
                            const arma::vec& item_time_intensity,
                            const arma::vec& person_speed,
                            const arma::vec& rt_resid_sd,
                            const arma::vec& rt_speed_loading) {
  return loglik_rt(log_rt, item_time_intensity, person_speed, rt_resid_sd, rt_speed_loading);
}

// [[Rcpp::export]]
List gibbs_joint_rt_irt_null(const arma::mat& response,
                             const arma::mat& log_rt,
                             int n_iter = 5000,
                             int n_burnin = 2500,
                             bool one_pl = false,
                             bool verbose = true,
                             int chain_id = 1) {
  return gibbs_rtirt_null(response, log_rt, n_iter, n_burnin, one_pl,
                          1.0, 1.0, -1.0, 1e6, 0.0, 0.0, 2,
                          false, 0.08, 0.12, false, 0.30, 25,
                          false, 1.0, 0.5, true, false, 0.30, 25,
                          verbose, chain_id);
}
