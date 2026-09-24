# tau_screen.R — score-based screen for speed heterogeneity along a covariate, from one full fit.
# Uses the casewise tau-location scores s_i = E[(tau_i - m_i) / v | data] returned by smi_rtirt
# (eta > 0). Scores for the location and for the covariates already in the tau mean (Xtau) are
# centred; the location score is projected off the Xtau scores (efficient score) and cumulated in
# the order of z. Under a correct tau mean the process is approximately a Brownian bridge:
# DM = sup |B| (double-max, any shape of drift), LM2 = median-split LM (a two-group comparison).
if (!exists("score_instability")) source("aghq_score.R")
tau_screen <- function(fit, z, Xtau = NULL) {
  s <- fit$tau_score
  if (is.null(s)) stop("fit has no tau_score (eta must be > 0)")
  if (!is.null(Xtau)) return(score_instability(scale(cbind(s, s * as.matrix(Xtau)), scale = FALSE), 1, z))
  s <- s - mean(s); n <- length(s)
  B <- cumsum(s[order(z)]) / sqrt(n * mean(s^2)); mid <- floor(n / 2)
  lm2 <- B[mid]^2 / ((mid / n) * (1 - mid / n))
  c(DM = max(abs(B)), p_DM = p_supbb(max(abs(B))), LM2 = lm2, p_LM2 = pchisq(lm2, 1, lower.tail = FALSE))
}
