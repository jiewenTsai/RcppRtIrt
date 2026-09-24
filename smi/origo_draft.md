# Origo draft — semi-modular inference for RT-IRT

Status: **Argument Gate NOT CLEAR** (see "Blocking issue"). Draft for discussion, not a fixed origo.

## P — problem (causal-necessity sentence)

Without a way to control how much response-time (RT) information enters ability estimation,
researchers who score examinees or compare groups with joint RT-IRT models must choose between the
full joint model and a two-stage (cut) analysis. The full model shifts the estimated ability gap
between groups that differ in speed but not in ability (here 0.11–0.16 SD), and the usual accuracy
criteria (overall RMSE, DIC) still prefer it. The two-stage analysis protects the scores but
attenuates the estimated speed–ability relation by the reliability of theta (≈ 0.80 here). Either
way, a reported quantity is wrong in a direction the standard workflow does not reveal.

## S — solution

Semi-modular inference for RT-IRT with:
1. marginal tempering of the RT module (tau integrated out), which makes eta a graded dial and
   avoids the impropriety of tempering a latent-variable prior;
2. target-specific choice of eta, using the cut posterior as the unbiased reference for the target
   (a Hausman-type conflict check and a plug-in posterior-risk rule);
3. separate eta for separate purposes: low / data-chosen eta for scores and group comparisons,
   eta → 1 for the speed–ability relation;
4. an efficient probit DA Gibbs implementation (z-scale PX-DA for items, C++ stage 2).

## E — evidence (empirical so-what chain)

- Finding: when one group is slower but equally able, the full model distorts the group gap by
  0.11–0.16 SD while its overall theta RMSE beats the cut (0.38 vs 0.42). The cut attenuates
  gamma to 0.80 of its value. The risk / Hausman rules pick eta ≈ 0.07–0.15 under a speed shift
  (gap RMSE 0.02 vs 0.11–0.16) and eta ≈ 0.7–0.9 without one (theta RMSE 0.37–0.38 vs 0.42).
- Consequence: the model-selection criterion, not the data, decides whether group comparisons are
  biased; a single joint model cannot serve both score reporting and speed–ability research.
- Decision: an analyst reporting group comparisons of RT-informed scores selects eta with the
  target-specific rule (and reports the share of RT information borrowed); an analyst estimating
  the speed–ability relation uses eta near 1 (or conditions on the needed variables).

## V — value (derived; provisional)

O-level: O2 (new boundary condition and procedure for an existing method), with an O1 component
(marginal tempering construction; impropriety of latent-prior tempering). Candidate outlets:
Psychometrika / BJMSP (if the tempering result and the risk rule are formalised), JEBS / JEM
(if the fairness framing leads).

Core papers (≤ 5; stances to confirm by reading):
1. van der Linden (2007), hierarchical framework — partial: assumes one speed distribution; group
   speed differences then leak into ability through the covariance.
2. Bolsinova & Tijmstra (2018), RT to improve ability precision — partial/oppose: the precision
   gain is real, but precision-based criteria hide group bias.
3. Frazier, Nott et al. (2023, JASA), posterior risk of modular / semi-modular inference —
   support/extend: we give a plug-in, target-specific risk rule for a latent-variable module.
4. Carmona & Nicholls (2020), semi-modular inference — partial: tempering with latent variables in
   the suspect module must be defined on the marginal likelihood.
5. Mislevy (1991), plausible values / conditioning models — support: the cut reproduces the
   attenuation of unconditioned plausible values. (Not checked via Consensus in this session.)
Related: Levy (2026, BJMSP) modular Bayesian IRT/SEM (no RT, no SMI).

Field-level value: without this, joint RT-IRT scoring either biases group comparisons or, when
cut, underestimates the speed–ability link, and neither error is visible to RMSE-based checks.

## Assumption audit (three weakest)

1. Group labels are known for the eta rule — (b) acknowledged; extension to continuous covariates
   via the score-based instability tests from `aghq_score.R`.
2. Hausman approximation Var(D) ≈ V(0) − V(eta) with posterior variances under priors —
   **(c) unaddressed**: no size study of the Hausman rule and no coverage study of the risk rule.
3. Only one misspecification type (group speed shift) plus mild item-level gamma heterogeneity —
   (b) acknowledged; would need at least cross-loading heterogeneity and a nonlinear
   speed–ability relation.
Also assumed: the trusted IRT module is correct (b).

## Update after `exp_hier_vs_cut.R` and `exp_choose_eta_cond.R`

- The cut is an unbiased reference for a group gap only when theta's prior conditions on that
  grouping (conditioning model). The eta rules are now run relative to cut+G.
- A hierarchical model with G in the theta and tau means is unbiased for G; the cut's added value
  is protection for comparisons on variables that are not in the RT model (secondary analyses).
- Risk rule: works for both targets (keeps eta ≈ 0.75 without misspecification, drops to ≈ 0.03
  under a Z speed shift). Hausman rule: oversized under the null (stops at eta = 0 in ~half of the
  datasets), because of Monte Carlo error at small eta. Drop it or correct it for MC error.
- Revised P candidate: RT-informed scores are reused for comparisons on variables the scoring
  model did not condition on; the joint model then leaks speed differences on those variables into
  ability, and no conditioning model can include every future analysis variable.

## Blocking issue

Assumption 2 is (c). Required revision: a null-scenario size study of the Hausman rule and a
coverage study of the eta chosen by the risk rule (≥ 100 datasets), before the procedure can be
claimed. P also needs an actor check: which analysts actually use RT-informed scores for group
comparisons (research vs operational), since most operational programs do not report them.

Gilbert positioning check not run: `references/gilbert-papers.md` of the mentor skill is not
available in this environment.
