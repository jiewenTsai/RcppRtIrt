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
// │ theta_i (能力)        │ Conjugate Gibbs: Normal + mean-centering          │
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
//   Sigma_p → omega(PG) → a → b → theta(+centering) → xi → sigma_t → tau
//   ★ 關鍵：a 必須在 b 之前更新
//      因為 draw_b 的 precision = a_j^2 * sum(omega)，
//      若 a 先被抽到很小（如 0.01），b 的 precision ≈ 0，
//      後驗 variance → ∞，b 就爆炸。
//      a 先更新可以讓 b 看到合理的 a，反之亦然。
//
// 識別限制：
//   - mean(theta) = 0：每次 iter 對 theta 做 mean-centering（不截斷到 ±10）
//   - a_j > 0：TruncNormal(0,∞) 保證（不需要額外限制）
//   - fix_int=TRUE：令 beta[intercept,] = 0（IRT 本身不需要截距）
//
// 依賴：Rcpp, RcppArmadillo, pg (tmsalab/pg)
//   安裝: R CMD INSTALL /path/to/pg
// =============================================================================

// [[Rcpp::depends(RcppArmadillo, pg)]]
#include <RcppArmadillo.h>
#include <pg.h>
using namespace Rcpp;
using namespace arma;

// =============================================================================
// SECTION 1: 工具函數
// =============================================================================

inline double logistic(double x) {
  return (x > 0.0) ? 1.0 / (1.0 + std::exp(-x))
                   : std::exp(x) / (1.0 + std::exp(x));
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
// 後驗精確推導：
//   -omega_ij/2 * a_j^2*(theta_i-b_j)^2 + kappa_ij*a_j*(theta_i-b_j)
//   完成二次方：
//   parV_j = 1/(1/sigma_a^2 + sum_i omega_ij*(theta_i-b_j)^2)
//   parM_j = parV_j*(mu_a/sigma_a^2 + sum_i kappa_ij*(theta_i-b_j))
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
// 對 b_j 的二次型（eta_ij = a_j*(theta_i - b_j) = a_j*theta_i - a_j*b_j）：
//   線性項：-kappa_ij * a_j * b_j  +  omega_ij/2 * 2*a_j^2*theta_i*b_j
//          = a_j*(-kappa_ij + a_j*theta_i*omega_ij) * b_j
//   二次項：-omega_ij/2 * a_j^2 * b_j^2
//   完成二次方：
//   parV_j = 1/(1/sigma_b^2 + a_j^2 * sum_i omega_ij)
//   parM_j = parV_j*(mu_b/sigma_b^2 + a_j * sum_i(a_j*theta_i*omega_ij - kappa_ij))
//
// ★ 符號陷阱：sum 裡是 (a_j*theta_i*omega_ij - kappa_ij)，不是加 kappa
//   之前寫成 (+kappa_ij + a_j*theta_i*omega_ij) 符號完全反了 → b → ±∞
// [[Rcpp::export]]
arma::vec draw_b(const arma::mat& kappa,
                  const arma::mat& omega,
                  const arma::vec& theta,
                  const arma::vec& a,
                  double mu_b    = 0.0,
                  double sigma_b = 2.0) {   // 合理先驗 N(0, 2^2)，不用 1e6 的 improper 先驗
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

// --- theta_i ~ Normal（無截斷）---
// 對 theta_i 的二次型（eta_ij = a_j*theta_i - a_j*b_j，線性係數 = a_j）：
//   parV_i = 1/(1/sigma2 + sum_j a_j^2*omega_ij)
//   parM_i = parV_i*(mu_i/sigma2 + sum_j a_j*(kappa_ij + a_j*b_j*omega_ij))
//
// ★ 不使用 (-10,10) 截斷：硬截斷偏移後驗，mean-centering 是更好的識別方式
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
                   const arma::vec& sigma_t,
                   double mu_xi = 4.0, double sigma_xi = 1e6) {
  int n = logT.n_rows, p = logT.n_cols;
  arma::vec xi(p);
  double inv_s2 = 1.0 / (sigma_xi * sigma_xi);
  for (int j = 0; j < p; ++j) {
    double inv_t2 = 1.0 / (sigma_t[j] * sigma_t[j]);
    double parV   = 1.0 / (inv_s2 + n * inv_t2);
    double parM   = parV * (mu_xi * inv_s2 + arma::sum(logT.col(j) + tau) * inv_t2);
    xi[j] = rtruncnorm_lo(parM, std::sqrt(parV), 0.0);
  }
  return xi;
}

// sigma_t_j ~ sqrt(InvGamma(da + n/2, db + SSE/2))
// da=db=0 → Jeffreys improper prior p(sigma^2) ∝ 1/sigma^2
// （等同 1/sigma 的先驗，對 log(sigma) 是 flat）
// [[Rcpp::export]]
arma::vec draw_sigma_t(const arma::mat& logT,
                        const arma::vec& xi,
                        const arma::vec& tau,
                        double da = 0.0, double db = 0.0) {
  int n = logT.n_rows, p = logT.n_cols;
  arma::vec sigma_t(p);
  for (int j = 0; j < p; ++j) {
    double par_a = da + n / 2.0;
    arma::vec r  = logT.col(j) - xi[j] + tau;
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
                    const arma::vec& mu_tau,
                    double sigma2_tau = 1.0) {
  int n = logT.n_rows;
  arma::vec pt    = 1.0 / (sigma_t % sigma_t);
  double sum_prec = arma::sum(pt);
  double inv_s2   = 1.0 / sigma2_tau;
  double parV     = 1.0 / (inv_s2 + sum_prec);
  arma::vec tau(n);
  for (int i = 0; i < n; ++i) {
    double sum_t = arma::dot(xi - logT.row(i).t(), pt);
    double parM  = parV * (mu_tau[i] * inv_s2 + sum_t);
    tau[i] = R::rnorm(parM, std::sqrt(parV));
  }
  return tau;
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
                  const arma::vec& tau, const arma::vec& sigma_t) {
  double ll = 0.0;
  int n = logT.n_rows, p = logT.n_cols;
  for (int j = 0; j < p; ++j) {
    double ls = std::log(sigma_t[j]);
    for (int i = 0; i < n; ++i) {
      double z = (logT(i,j) - xi[j] + tau[i]) / sigma_t[j];
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

arma::mat make_design(const arma::mat& X) {
  return arma::join_horiz(arma::ones(X.n_rows,1), X);
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
                       bool   use_mh_ab = false, // FALSE=conjugate Gibbs (recommended); TRUE=intercept reparam (experimental)
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
  arma::mat Sigma_p = arma::eye(2, 2);
  arma::vec mu0(n, arma::fill::zeros);

  arma::mat post_theta(n_save,n), post_tau(n_save,n);
  arma::mat post_a(n_save,p),     post_b(n_save,p);
  arma::mat post_xi(n_save,p),    post_st(n_save,p);
  arma::mat post_Sp(n_save,4);    arma::vec post_ll(n_save);
  int idx = 0;

  if (verbose)
    Rcpp::Rcout << "Chain " << chain_id << " [rtirt_null] "
                << n_iter << " iter, burnin=" << n_burnin << "\n";

  for (int m = 0; m < n_iter; ++m) {
    Rcpp::checkUserInterrupt();
    if (verbose && (m % std::max(1,n_iter/50)==0 || m==n_iter-1))
      progress_bar(m, n_iter, n_burnin, chain_id);

    // 1. Sigma_p（IW posterior）
    Sigma_p = draw_Sigma_null(theta, tau);

    // 2. PG 輔助
    arma::mat omega = draw_pg_irt(theta, a, b);

    // 3. IRT item parameters
    //    use_mh_ab=TRUE:  intercept reparam (a, c=a*b) → 打破 a-b 高相關
    //    use_mh_ab=FALSE: 個別 conjugate Gibbs
    if (use_mh_ab && !one_pl) {
      List ac = draw_ac_reparam(kappa, omega, theta, a, b, mu_a, sigma_a);
      a = as<arma::vec>(ac["a"]);
      b = as<arma::vec>(ac["b"]);
    } else {
      if (!one_pl) a = draw_a(kappa, omega, theta, b, mu_a, sigma_a);
      b = draw_b(kappa, omega, theta, a);
    }

    // 4. 能力 + mean-centering（識別限制）
    theta = draw_theta(kappa, omega, a, b, mu0, Sigma_p(0,0));
    theta -= arma::mean(theta);

    // 5. RT 項目
    xi      = draw_xi(logT, tau, sigma_t, mu_xi_use, sigma_xi);
    sigma_t = draw_sigma_t(logT, xi, tau, delta_a, delta_b);

    // 6. 速度
    tau = draw_tau(logT, xi, sigma_t, mu0, Sigma_p(1,1));

    if (m >= n_burnin) {
      post_theta.row(idx) = theta.t();
      post_tau.row(idx)   = tau.t();
      post_a.row(idx)     = a.t();
      post_b.row(idx)     = b.t();
      post_xi.row(idx)    = xi.t();
      post_st.row(idx)    = sigma_t.t();
      post_Sp.row(idx)    = arma::vectorise(Sigma_p).t();
      post_ll(idx) = loglik_irt(Y,theta,a,b) + loglik_rt(logT,xi,tau,sigma_t);
      ++idx;
    }
  }
  if (verbose) Rcpp::Rcout << "\n  Done.\n";

  return List::create(
    Named("theta")   = post_theta, Named("tau")     = post_tau,
    Named("a")       = post_a,     Named("b")       = post_b,
    Named("xi")      = post_xi,    Named("sigma_t") = post_st,
    Named("Sigma_p") = post_Sp,    Named("loglik")  = post_ll
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
  arma::mat beta(q, 2, arma::fill::zeros);
  arma::mat Sigma_p = arma::eye(2, 2);

  arma::mat post_theta(n_save,n), post_tau(n_save,n);
  arma::mat post_a(n_save,p),     post_b(n_save,p);
  arma::mat post_xi(n_save,p),    post_st(n_save,p);
  arma::mat post_beta(n_save,q*2), post_Sp(n_save,4);
  arma::vec post_ll(n_save);
  int idx = 0;

  if (verbose)
    Rcpp::Rcout << "Chain " << chain_id << " [rtirt] "
                << n_iter << " iter, burnin=" << n_burnin << "\n";

  for (int m = 0; m < n_iter; ++m) {
    Rcpp::checkUserInterrupt();
    if (verbose && (m % std::max(1,n_iter/50)==0 || m==n_iter-1))
      progress_bar(m, n_iter, n_burnin, chain_id);

    beta    = draw_beta(X, theta, tau, Sigma_p);
    if (fix_int) beta.row(0).zeros();
    Sigma_p = draw_Sigma_struct(X, theta, tau, beta);

    arma::mat omega = draw_pg_irt(theta, a, b);
    if (!one_pl) a  = draw_a(kappa, omega, theta, b, mu_a, sigma_a);
    b               = draw_b(kappa, omega, theta, a);

    arma::vec mu_theta = X * beta.col(0);
    theta = draw_theta(kappa, omega, a, b, mu_theta, Sigma_p(0,0));
    theta -= arma::mean(theta);

    xi      = draw_xi(logT, tau, sigma_t, mu_xi_use, sigma_xi);
    sigma_t = draw_sigma_t(logT, xi, tau, delta_a, delta_b);

    arma::vec mu_tau = X * beta.col(1);
    tau = draw_tau(logT, xi, sigma_t, mu_tau, Sigma_p(1,1));

    if (m >= n_burnin) {
      post_theta.row(idx) = theta.t();
      post_tau.row(idx)   = tau.t();
      post_a.row(idx)     = a.t();
      post_b.row(idx)     = b.t();
      post_xi.row(idx)    = xi.t();
      post_st.row(idx)    = sigma_t.t();
      post_beta.row(idx)  = arma::vectorise(beta).t();
      post_Sp.row(idx)    = arma::vectorise(Sigma_p).t();
      post_ll(idx) = loglik_irt(Y,theta,a,b) + loglik_rt(logT,xi,tau,sigma_t);
      ++idx;
    }
  }
  if (verbose) Rcpp::Rcout << "\n  Done.\n";

  return List::create(
    Named("theta")   = post_theta, Named("tau")     = post_tau,
    Named("a")       = post_a,     Named("b")       = post_b,
    Named("xi")      = post_xi,    Named("sigma_t") = post_st,
    Named("beta")    = post_beta,  Named("Sigma_p") = post_Sp,
    Named("loglik")  = post_ll
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

  arma::mat post_theta(n_save,n), post_a(n_save,p), post_b(n_save,p);
  arma::mat post_beta(n_save,q);  arma::vec post_ll(n_save);
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
    if (!one_pl) a  = draw_a(kappa, omega, theta, b, mu_a, sigma_a);
    b               = draw_b(kappa, omega, theta, a);

    arma::vec mu_theta = X * beta_ra;
    theta = draw_theta(kappa, omega, a, b, mu_theta, 1.0);
    theta -= arma::mean(theta);

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
    Named("loglik")= post_ll
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
