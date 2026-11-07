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
