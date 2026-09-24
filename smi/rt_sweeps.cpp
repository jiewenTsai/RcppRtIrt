// rt_sweeps.cpp — K Gibbs sweeps of the RT module given theta (SMI stage 2, eta = 1).
// Same model, priors and update order as rt_block() in smi_gibbs.R:
//   log T_ij = xi_j - tau_i + gamma_j theta_i + eps,  eps ~ N(0, sigma_j^2),  tau_i ~ N(0, v)
//   xi ~ N(4, 10^2), gamma ~ N(0, 1), sigma^2 ~ IG(1, 1), v ~ IG(1, 1)
// Uses R's RNG, so set.seed() in R controls it.
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
List rt_sweeps_cpp(const arma::vec& th, const arma::mat& logT,
                   arma::vec tau, arma::vec xi, arma::vec gam, arma::vec s2, double v, int K) {
  const int n = logT.n_rows, p = logT.n_cols;
  arma::mat P0(2, 2, arma::fill::zeros); P0(0, 0) = 1.0 / 100.0; P0(1, 1) = 1.0;
  arma::vec m0 = {4.0, 0.0};
  arma::vec P0m0 = P0 * m0;
  arma::mat W(n, 2); W.col(0).ones(); W.col(1) = th;
  const arma::mat WtW = W.t() * W;
  for (int k = 0; k < K; ++k) {
    // tau_i | rest
    arma::vec w = 1.0 / s2;
    const double prec = 1.0 / v + arma::sum(w);
    const double sdt = 1.0 / std::sqrt(prec);
    for (int i = 0; i < n; ++i) {
      double r = 0.0;
      for (int j = 0; j < p; ++j) r += w[j] * (xi[j] + gam[j] * th[i] - logT(i, j));
      tau[i] = r / prec + sdt * R::rnorm(0.0, 1.0);
    }
    // (xi_j, gamma_j) | rest: log T_j + tau = xi_j + gamma_j theta + eps
    arma::mat Yt = logT; Yt.each_col() += tau;
    const arma::mat Wty = W.t() * Yt;                      // 2 x p
    for (int j = 0; j < p; ++j) {
      arma::mat Vj = arma::inv_sympd(P0 + WtW / s2[j]);
      arma::vec mj = Vj * (P0m0 + Wty.col(j) / s2[j]);
      arma::vec z = {R::rnorm(0.0, 1.0), R::rnorm(0.0, 1.0)};
      arma::vec xg = mj + arma::chol(Vj, "lower") * z;
      xi[j] = xg[0]; gam[j] = xg[1];
    }
    // sigma_j^2 | rest
    for (int j = 0; j < p; ++j) {
      double ss = 0.0;
      for (int i = 0; i < n; ++i) {
        const double e = logT(i, j) - (xi[j] - tau[i] + gam[j] * th[i]);
        ss += e * e;
      }
      s2[j] = 1.0 / R::rgamma(1.0 + n / 2.0, 1.0 / (1.0 + ss / 2.0));
    }
    // v | tau
    v = 1.0 / R::rgamma(1.0 + n / 2.0, 1.0 / (1.0 + arma::dot(tau, tau) / 2.0));
  }
  return List::create(_["tau"] = tau, _["xi"] = xi, _["gam"] = gam, _["s2"] = s2, _["v"] = v);
}
