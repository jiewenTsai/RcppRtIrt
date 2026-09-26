# Origo — semi-modular inference for RT-IRT

Status: **Argument Gate CLEAR** (gate block below). Earlier drafts and the reasons they
did not clear follow it for the record.

## Gate block

```
═══════════════════════════════════════════════
Argument Gate: CLEAR
───────────────────────────────────────────────
P  : Without a target-specific check on how much response-time information enters ability
     estimation, analysts who compare groups with joint RT-IRT scores (group gaps in large-scale
     assessments, RT-informed EAP scoring in CAT) will report gaps shifted by the groups' speed
     differences: 0.07-0.13 SD for person-level speed differences on variables outside the
     scoring model, and up to 0.2-0.44 SD when a group is slower on a few items (differential
     response time), even for variables the hierarchical model conditions on. Precision criteria
     (theta RMSE, reliability) favour the joint model in all of these cases, so the error is not
     visible in the standard workflow.

S  : Semi-modular RT-IRT (marginal tempering of the RT module, conditioning model for reporting
     variables) with a target-specific choice of the influence parameter eta: the eta that
     minimises the plug-in posterior risk of the target relative to the conditioned cut
     posterior; plus a one-fit score-based screen for person-level speed drift along a covariate.

E  : Empirical so-what chain (simulation). Finding: the risk rule keeps eta ≈ 0.6-0.8 when there
     is no speed difference (theta RMSE 0.38-0.39 vs 0.44 for the cut) and drops to ≈ 0-0.13 under
     person-level or item-level speed differences, with gap bias ≤ 0.013 (full model: 0.07-0.44)
     and coverage .88-.96 (full: 0-.68). RT adds almost no precision to group gaps (posterior SD
     0.042 vs 0.044). Consequence: the joint model's accuracy advantage does not certify group
     comparisons, and a hierarchical conditioning model does not protect against item-level speed
     differences. Decision: an analyst reporting a group comparison from joint RT-IRT scores
     selects eta for that comparison with the risk rule and reports it (default: the cut when no
     RT data or covariate is at hand); a program releasing scores for secondary analysis releases
     response-only (cut) plausible values for group comparisons, because RT buys no gap precision
     and post-hoc correction fails under differential response time.

V  : O2 contribution (boundary condition on RT as collateral information + a procedure).
     Core papers:
       1. van der Linden, Klein Entink & Fox (2010, APM) — partial — reason: collateral RT
          information improves accuracy when the speed model holds for everyone; group speed
          differences and DRT turn that information into group bias.
       2. Bolsinova & Tijmstra (2018, BJMSP) — partial/oppose — reason: precision gains are
          real, but RMSE-type evaluation hides the group bias shown here.
       3. Kern et al. (2021, APM) — partial — reason: J-EAP is evaluated on individual accuracy;
          its use for group comparisons needs the target-specific check.
       4. Frazier, Nott et al. (2023/2025, JASA) — support/extend — reason: posterior-risk
          selection of semi-modular posteriors; we give a plug-in, target-specific rule for a
          latent-variable RT module and show when it is needed.
       5. Levy (2024, JEBS; 2026, BJMSP) — support/extend — reason: measurement-preserving
          (cut) multistage Bayesian IRT; we add graded feedback via eta and a data-driven choice,
          and show the RT module needs it.
     Context (not core): plausible-value conditioning bias for omitted variables is attenuation
     (Monseur & Adams 2009; Bailey et al. 2023), whereas RT leakage is additive and survives
     conditioning under DRT; real speed differences by gender and multilingual status
     (Kapoor et al. 2024, JEM; Park et al. 2024, EMIP).
     Field-level value: joint RT-IRT scoring is currently evaluated by precision alone; without
     this paper, group comparisons made with RT-informed scores inherit speed differences
     invisibly, and the hierarchical model's usual remedy (conditioning) does not remove them.
───────────────────────────────────────────────
Assumption audit (three weakest):
  1. The conditioned cut is an unbiased reference for the target (IRT module and conditioning
     model correct for that target) — (b) limitation; the same assumption as plausible values.
  2. Post-selection coverage of the risk rule is .88-.94, slightly below nominal — (b) reported.
  3. Simulation conditions (lognormal RT, one speed factor, n = 500, p = 10) — (b); an empirical
     illustration (e.g. PISA 2018 multilingual-learner or gender gaps with log data) is the next
     required piece for APM/JEM.
  Resolved: Hausman rule dropped (oversized); post-hoc correction scoped out (fails under DRT).
Weakest defended assumption: assumption 2 (post-selection coverage).
Gilbert positioning: no duplication (Gilbert et al. 2026, EPM studies RT–discrimination, not RT
  as collateral information for scores); closest contrast is Gilbert, Soland & Domingue (2026,
  EMIP), which shows scoring decisions shift downstream inference — this paper supplies a
  decision rule for one such scoring decision. Same O2 level.
Journal: JEM (group comparisons / fairness framing, needs the empirical illustration) or APM.
Next: empirical illustration, then mvp-writer for the theory scaffold.
═══════════════════════════════════════════════
```

---

# Earlier drafts (not cleared)


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
