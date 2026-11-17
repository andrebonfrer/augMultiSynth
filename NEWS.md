## augMultiSynth 0.3.8
- `fit_weights_fw()` now converges. It previously took plain Frank–Wolfe steps and stopped when the
  objective changed by less than `tol` between iterations; `multiout_synth()` called it with
  `max_iter = 1000, tol = 1e-5`. In benchmarks on the simulated panel (100 and 500 donors) this left the
  objective up to 3.5% and 5.9% above its minimum and the weights up to 0.32 and 0.64 (L1) from the
  exact solution. The solver now uses pairwise Frank–Wolfe steps with exact line search and stops on the
  Frank–Wolfe duality gap, which bounds the distance to the optimal objective (default relative gap
  `1e-8`, `max_iter = 50000`). It matches an exact QP solution to about 1e-14 (relative objective) and
  is 3-12x faster per unit. **Weights, and therefore tau, change for existing fits.**
- `fit_weights_fw()` gains `w0` (warm start) and returns attributes `fw_gap`, `fw_rel_gap`, `fw_iter`,
  `fw_converged`.
- `multiout_synth()` gains `fw_max_iter` and `fw_tol`, returns `fw_diagnostics` (per-unit relative gap,
  iterations, convergence), and warns if any unit did not converge.
- `tune_lambda()` defaults `max_iter = 50000, tol = 1e-6`; `tune_nu()` uses `tol = 1e-6` (gap-based).
- Matrix is no longer imported (it was only used by the old solver).

## augMultiSynth 0.3.7
- New `tau_inference()`: unit-level estimate, standard error, t, p-value and interval for the average
  post-treatment effect `mean(fit$tau[j, m, post + 1])` (default `post = 1:K`, as in `eval_avg_1K()`).
  Inference conditions on donor paths and weights: the pseudo-panel likelihood factorises into
  p(Y_j | Y_D; tau, mu, sigma^2) p(Y_D), and the first factor is the regression of the gap series
  g(t) = Y_j(t) - w'Y_D(t) on an intercept and a post indicator over event times -L..-1 and `post`.
  Estimates equal the average of `fit$tau` to machine precision; standard errors equal OLS on the gap series.
  A weighted regression on the pseudo-panel (treated unit plus donors) with unit and event-time fixed
  effects reproduces the point estimate, but not the standard errors, because donor rows enter as
  observations; without the fixed effects the pooled regression is biased by the treated-vs-synthetic
  level gap.
- New `scm_gap_panel()`: the gap series in long form, for plotting or for any sampler with an identity
  weight matrix. For Bayesian inference, scmBayesPost >= 0.4.8 builds the same likelihood directly with
  `prepare_data_general(likelihood = "conditional", L, post)`.
- `multiout_synth()` now also returns `L` and `K`.
- `fit_weights_fw()`: the line-search denominator no longer calls `crossprod()` on a Matrix object, which
  failed on R < 4.4 ("requires numeric/complex matrix/vector arguments"). Results are unchanged.
- `run_demo(tune_nu = TRUE)` now fits with the selected `nu`. Previously `tune_nu()` ran but its `nu_opt`
  was discarded and the fit used the `nu` argument. `run_demo()` also returns `nu_used` and the
  `tune_nu()` result (`tune_nu`).
- `tune_nu()` is exported (it was documented with `@export` but missing from NAMESPACE); stale Rd files
  regenerated.
- testthat suite added (`tests/testthat/test-tau-inference.R`).

## augMultiSynth 0.3.6
- `multiout_synth(donors_list = )`: optional per-unit candidate donor sets (unit indices or character ids, keyed by
  treated unit id or aligned with `treated_units`). Screening runs inside the set; the eligibility rule still applies;
  units without an entry behave as before. Lets a caller impose design-specific pools, e.g. donors adopting within a
  band of weeks around the treated unit's relative adoption week (budgetsc's two-launch design). Default `NULL`
  reproduces 0.3.5 exactly.
- `donor_screen(candidates = )`: the same restriction at the screening level.
- `tune_lambda()`: with `intercept = "outcome"` the per-outcome intercept is now estimated on the training pre-periods
  and applied to the held-out periods. Previously the held-out block was re-centred on its own mean (via
  `demean_within_outcome_blocks()` on the validation rows), so a level divergence between fitted and held-out weeks
  did not count as validation error and the criterion scored shape only. Validation MSEs are therefore larger than
  before and lambda choices may differ; `intercept = "global"` was already correct. `tune_lambda()` also accepts
  `donors_list`.
- Internal `resolve_donors_list()`.

## augMultiSynth 0.3.1
- Fix scoping bug in parallel backend where M was not defined in multiout_synth().
