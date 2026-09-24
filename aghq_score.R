# aghq_score.R — Gibbs + AGHQ 的 score-based 參數不穩定性檢定（原型）
#
# 流程：
#   1. gibbs_rtirt_null() 取得後驗平均 → 轉成識別化參數（Sigma_11 = 1）
#   2. 每人的 marginal log-likelihood：tau 解析積分，theta 用一維 AGHQ
#   3. 以 Gibbs 估計為起點做 BHHH，得到 MML 估計（score test 需要在 MLE 上評估）
#   4. 每人 score（數值微分）→ 對目標參數做 efficient score → 依共變數排序的
#      累積過程 → double-max (DM) 檢定，與「事先分兩組」的 LM 檢定比較
#
# 模型（null RT-IRT，slope_rt = 1）：
#   logit P(Y_ij = 1) = a_j (theta_i - b_j)
#   log T_ij ~ N(xi_j - tau_i, sigma_j^2)
#   (theta_i, tau_i) ~ BVN(0, Sigma),  Sigma_11 = 1
#
# 參數向量 psi = (a[p], b[p], xi[p], log sigma[p], Sigma_12, log Sigma_22)

# --- Gauss–Hermite 節點（physicists'，權重 e^{-x^2}），Golub–Welsch ---
gh_nodes <- function(K) {
  i <- seq_len(K - 1)
  J <- matrix(0, K, K)
  J[cbind(i, i + 1)] <- J[cbind(i + 1, i)] <- sqrt(i / 2)
  e <- eigen(J, symmetric = TRUE)
  list(x = e$values, w = sqrt(pi) * e$vectors[1, ]^2)
}

unpack_psi <- function(psi, p) {
  list(a = psi[1:p], b = psi[p + 1:p], xi = psi[2 * p + 1:p],
       lsig = psi[3 * p + 1:p], s12 = psi[4 * p + 1], ls22 = psi[4 * p + 2])
}

log1pexp <- function(x) ifelse(x > 0, x + log1p(exp(-x)), log1p(exp(x)))

# 每人的 marginal log-likelihood（長度 n 的向量）
#   tau 解析積分：令 r_ij = xi_j - logT_ij = tau_i + e_ij，w_j = 1/sigma_j^2，
#   sum_j w_j (r_ij - tau)^2 = sum_j w_j (r_ij - rbar_i)^2 + W (tau - rbar_i)^2，
#   而 tau | theta ~ N(S12 * theta, S22 - S12^2)，所以
#   p(t_i | theta) = C_i * N(rbar_i; S12 * theta, S22 - S12^2 + 1/W)
marg_ll <- function(psi, Y, logT, gh, return_mode = FALSE) {
  n <- nrow(Y); p <- ncol(Y)
  P <- unpack_psi(psi, p)
  a <- P$a; d <- P$a * P$b
  w <- exp(-2 * P$lsig); W <- sum(w)
  R <- matrix(P$xi, n, p, byrow = TRUE) - logT
  rbar <- drop(R %*% w) / W
  C <- -0.5 * p * log(2 * pi) - sum(P$lsig) -
       0.5 * drop(((R - rbar)^2) %*% w) + 0.5 * log(2 * pi / W)
  beta <- P$s12
  v <- exp(P$ls22) - beta^2
  if (!is.finite(v) || v <= 0) return(rep(-Inf, n))
  vr <- v + 1 / W

  logint <- function(th) {           # log 被積函數（對 theta）
    eta <- outer(th, a) - matrix(d, n, p, byrow = TRUE)
    rowSums(Y * eta - log1pexp(eta)) +
      dnorm(rbar, beta * th, sqrt(vr), log = TRUE) + dnorm(th, 0, 1, log = TRUE)
  }

  # 眾數：Newton（被積函數對 theta 為嚴格凹）
  th <- rep(0, n)
  for (it in 1:50) {
    Pm <- plogis(outer(th, a) - matrix(d, n, p, byrow = TRUE))
    g1 <- drop((Y - Pm) %*% a) + beta * (rbar - beta * th) / vr - th
    g2 <- -drop((Pm * (1 - Pm)) %*% (a^2)) - beta^2 / vr - 1
    step <- g1 / g2
    th <- th - step
    if (max(abs(step)) < 1e-10) break
  }
  sd <- 1 / sqrt(-g2)
  if (return_mode) return(list(mode = th, sd = sd))

  # AGHQ：∫ e^{g} ≈ √2·sd Σ_k w_k exp(g(θ̂ + √2·sd·x_k) + x_k²)
  L <- sapply(seq_along(gh$x), function(k)
    log(gh$w[k]) + gh$x[k]^2 + logint(th + sqrt(2) * sd * gh$x[k]))
  L <- matrix(L, nrow = n)
  mx <- apply(L, 1, max)
  C + log(sqrt(2) * sd) + mx + log(rowSums(exp(L - mx)))
}

# 每人 score：中央差分，回傳 n x q 矩陣
casewise_scores <- function(psi, Y, logT, gh, h = 1e-5) {
  sapply(seq_along(psi), function(k) {
    e <- replace(numeric(length(psi)), k, h)
    (marg_ll(psi + e, Y, logT, gh) - marg_ll(psi - e, Y, logT, gh)) / (2 * h)
  })
}

# BHHH：psi ← psi + (S'S)^{-1} S'1，附步長減半
bhhh_mml <- function(psi, Y, logT, gh, max_iter = 50, tol = 1e-6, verbose = FALSE) {
  ll <- sum(marg_ll(psi, Y, logT, gh))
  for (it in 1:max_iter) {
    S <- casewise_scores(psi, Y, logT, gh)
    g <- colSums(S)
    step <- solve(crossprod(S) + 1e-8 * diag(ncol(S)), g)
    lam <- 1
    repeat {
      cand <- psi + lam * step
      ll_new <- sum(marg_ll(cand, Y, logT, gh))
      if (is.finite(ll_new) && ll_new >= ll - 1e-10) break
      lam <- lam / 2
      if (lam < 1e-6) break
    }
    psi <- cand
    if (verbose) cat(sprintf("  BHHH %2d  ll = %.4f  max|g|/n = %.2e\n",
                             it, ll_new, max(abs(g)) / nrow(Y)))
    if (abs(ll_new - ll) < tol) { ll <- ll_new; break }
    ll <- ll_new
  }
  list(psi = psi, ll = ll, iter = it)
}

# Gibbs 後驗平均 → 識別化的 psi（Sigma_11 = 1）
psi_from_gibbs <- function(fit) {
  a  <- colMeans(fit$a);  b <- colMeans(fit$b)
  xi <- colMeans(fit$xi); st <- colMeans(fit$sigma_t)
  Sp <- colMeans(fit$Sigma_p)          # (S11, S21, S12, S22)
  s  <- sqrt(Sp[1])
  c(a * s, b / s, xi, log(st), Sp[3] / s, log(Sp[4]))
}

# Brownian bridge 的 sup|B| 分布（Kolmogorov）
p_supbb <- function(x, K = 100) {
  k <- 1:K
  min(1, max(0, 2 * sum((-1)^(k + 1) * exp(-2 * k^2 * x^2))))
}

# 對參數 k 做 score-based 不穩定性檢定（沿共變數 z 排序）
#   efficient score：扣掉對 nuisance 的投影後再標準化 → 在 H0 下為 Brownian bridge
score_instability <- function(S, k, z) {
  n <- nrow(S)
  J <- crossprod(S) / n
  s_eff <- S[, k] - S[, -k, drop = FALSE] %*% solve(J[-k, -k], J[-k, k])
  s_eff <- drop(s_eff)
  ord <- order(z)
  B <- cumsum(s_eff[ord]) / sqrt(n * mean(s_eff^2))
  t <- seq_len(n) / n
  dm <- max(abs(B))
  mid <- floor(n / 2)                                 # 事先在中位數切兩組的 LM
  lm2 <- B[mid]^2 / (t[mid] * (1 - t[mid]))
  c(DM = dm, p_DM = p_supbb(dm), LM2 = lm2, p_LM2 = pchisq(lm2, 1, lower.tail = FALSE))
}
