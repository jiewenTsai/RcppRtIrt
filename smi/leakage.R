# leakage.R — leakage of speed differences into ability under the (tempered) joint model.
# Gaussian approximation per person, item parameters at their posterior means. With RT terms
# multiplied by eta (marginal tempering), the joint posterior precision of (theta_i, tau_i) is
#   P = [[I_i + eta sum g^2 w, -eta sum g w], [-eta sum g w, eta / v + eta sum w]],  w = 1 / sigma^2,
# where I_i = 1 + sum_j a_j^2 phi^2 / (Phi (1 - Phi)) is the probit information plus the prior.
# A speed shift delta that the tau mean does not model (person slower by delta) moves the RT
# residuals by +delta, i.e. b by eta delta (sum g w, -sum w), so (theta_hat, tau_hat) move by
# P^{-1} eta delta (sum g w, -sum w). Two summaries:
#   lambda = d theta_hat / d delta    (bias of a theta gap per unit of unmodelled speed gap)
#   c      = d theta_hat / d tau_hat  (leakage ratio: theta shift per unit shift in the speed residual)
# A secondary analyst with draws of theta and u = tau - m_tau can remove leakage on any Z by
#   gap_Z(theta) - c * gap_Z(u)   (computed per draw).
leakage <- function(fit, eta = fit$eta) {
  a <- colMeans(fit$a); d <- fit$d_mean; g <- colMeans(fit$gamma_stage1)
  w <- 1 / fit$s2_stage1_mean; v <- mean(fit$v_stage1)
  th <- colMeans(fit$theta)
  lp <- outer(th, a) - matrix(d, length(th), length(a), byrow = TRUE)
  I <- 1 + drop((dnorm(lp)^2 / (pnorm(lp) * pnorm(-lp))) %*% a^2)
  Sgw <- sum(g * w); Sw <- sum(w); Sg2w <- sum(g^2 * w)
  P11 <- I + eta * Sg2w; P12 <- -eta * Sgw; P22 <- eta / v + eta * Sw
  det <- P11 * P22 - P12^2
  dth <- eta * (P22 * Sgw - P12 * (-Sw)) / det        # [P^{-1} db]_1 per unit delta
  dta <- eta * (-P12 * Sgw + P11 * (-Sw)) / det       # [P^{-1} db]_2
  list(lambda = mean(dth), c = mean(dth / dta), lambda_i = dth, c_i = dth / dta)
}
