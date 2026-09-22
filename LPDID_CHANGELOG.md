# Changelog for -lpdid- Stata command

This document lists all significant changes to `lpdid`, the Stata command implementing the Local Projections
Difference-in-Differences (LP-DiD) estimator of Dube, Girardi, Jordà and Taylor (*Journal of Applied
Econometrics*, 2025).

`lpdid` is written by Alexander Busch (Massachusetts Institute of Technology) and Daniele
Girardi (King's College London), in collaboration with Arin Dube, Òscar Jordà and Alan M.
Taylor. Versions 1.0.0 through 1.0.2 were written by both authors; versions 1.0.3 and 1.1.0
were written by Daniele Girardi.

Source and examples: <https://github.com/danielegirardi/lpdid>

To update to the latest release:

```stata
ssc install lpdid, replace
```

---

## Version 1.1.0 (September 2026)

The largest release since 1.0.0. It adds a new estimator of the overall average effect 
(in addition to the pooled LP-DiD estimates),  a pre-trend test, and a number of 
corrections to the definition of the estimation sample. Some
of those corrections can change results relative to previous versions.

### Main changes to existing behaviour

- **The pooled estimates are no longer reported by default.** They now require the new
  `pooled` option. By default, `lpdid` now prints only the event study. As a result,
  `e(pooled_results)` is no longer produced by default, 
  but only if the `pooled` option is specified.

  To restore the previous output (event study + pooled estimates), add `pooled`:

  ```stata
  lpdid Y, unit(id) time(t) treat(d) pre_window(10) post_window(10) pooled
  ```

  `pooled` is implied by `only_pooled`, `pre_pooled()` and `post_pooled()`, each of which is
  interpreted as a request for the pooled estimate. The `only_event` option is retained for backward compatibility
  and still suppresses the pooled estimates, so do-files that use it are unaffected; 
  by default it now does nothing.

- **`nocomp` is now stricter in holding the estimation sample fixed across horizons.** 
  The `nocomp` option rules out composition effects across the event window, ensuring that the set of contributing 
  observations is the same across every post- and pre-treatment horizon in the window. In previous
  releases, it ruled out composition changes caused by a control entering treatment, but let the sample vary
  whenever the outcome at a given horizon was unavailable — from a missing value, or
  from the window running off the end of the panel. It now additionally requires the outcome to
  be observed at every horizon, so a row enters every horizon or none.

  This makes the sample fully fixed across horizons also 
  in settings with missing values or truncated windows. 
  It comes with a cost in terms of sample size.
  
- **The clean-control sample is now stricter at the panel edge, under `nonabsorbing()`.** 
  To count as a clean control at horizon *h*, a row must have had no treatment switch in the preceding
  *L* periods (where *L* is the numerical argument of the `nonabsorbing()` option). 
  At the start of a unit's panel those lags do not exist, so they cannot be
  checked. 1.0.3 admitted such rows, effectively assuming that no treatment events occurred
  before the start of the sample period. This assumption is now not made by default, and such units 
  are now excluded because their clean control condition cannot be verified. 
  The new option `untreated_before` restores the previous behaviour and its underlying assumption.
  
- **The clean-control condition is now stricter in presence of missing values, under `nonabsorbing()`.** 
  Under persistent non-absorbing treatment, a unit whose treatment status was *missing* at some
  future period inside the window was admitted as a clean control. It is now excluded. 

- **An inverted pooled window is now refused.** A two-number argument given the wrong way round
  — `post_pooled(4 2)`, or `pre_pooled(5 2)` — previously ran and returned a (probably) incorrect
  estimate. Both now exit with r(198). This is the only change in the release that takes a
  specification which previously ran and makes it an error.

### New options

- **`aggregate_average`** — reports the aggregate average effect across all treated observations in 
  the post-treatment window, and the corresponding average over the pre-treatment window, each with a standard
  error, a p-value and a confidence interval. It averages the event-study coefficients,
  weighting each horizon by the number of treated observations behind it. With `rw` it targets the ATT;
  without `rw`, the variance-weighted counterpart (VWATT). Saved in `e(aggregate)`.

- **`pretrend_test`** — a joint *F* test of the null that all pre-treatment coefficients are
  zero, reported beneath the table and saved in four scalars. 

- **`pooled`** — requests the pooled LP-DiD estimates, which 1.0.3 and earlier reported by
  default. See *Main changes* above.

- **`untreated_before`** — assume every unit is untreated before it enters the panel,
  admitting observations whose earlier treatment history cannot be verified. Relevant only with
  `nonabsorbing()`. 

### New saved results

- `e(aggregate)` — the aggregate average estimates, in the same seven columns as
  `e(results)`, one row per window. Produced only if the new option `aggregate_average` is selected. 
  Note that its last (observation) column counts treated observations only, 
  where `e(results)` and `e(pooled_results)` count all observations in the regression.
- `e(pretrend_F)`, `e(pretrend_p)`, `e(pretrend_df)`, `e(pretrend_df_r)` — produced only if
  the new option `pretrend_test` is selected.

### Other changes that could change results

- **Horizons the estimator could not fit are now reported as missing.** 
  Where a horizon has no
  clean control observation, so that the treatment indicator is collinear with the time
  effects, 1.0.3 reported a coefficient of 0 with a standard error of 0. Those cells now come back missing and
  are named in a note, in the event-study and pooled tables alike. An aggregate average whose
  window contains a horizon that could not be estimated is suppressed.

  The case of a pre-treatment horizon that is zero by construction (for example because the specification 
  controls for lagged outcomes) is treated differently: it keeps its 0 value and produces an explanatory note.

- **Strengthened the checks to make sure that Regression Adjustment does not extrapolate into cells with no clean control.** 
  The regression-adjustment route (`rw` with covariates) imputes a counterfactual for treated
  observations, so it needs a clean control in every cell it conditions on. Where a cell has
  none, Regression Adjustment might extrapolate from the other cells and report a wrong number. 
  Previous versions already included an algorithm to avoid this. It has now been strengthened
  to catch cases that previous versions might have missed. In the vast majority of cases, this
  will produce no change at all.

- **Numerical precision.** 
  The LP-DiD outcome variable and, where `rw` is not used, the weight
  variables, were stored as `float`. They are now `double` throughout. 
  This can cause very small changes in point estimates.
  Nothing about inference changes.
  
- **A single unfeasible regression no longer aborts the whole command.**
  Horizons that cannot be estimated are now reported as missing, but the rest of the run completes
  and results for the remaining horizons are reported.

### Bug fixes

- When `nonabsorbing` was selected with both `oneoff` and `firsttreat`, and neither `notyet` nor `nevertreated`, 
  treatment episodes were miscounted, resulting in usable observations being discarded. 
  Fixed in version 1.1.0
- Under the `bootstrap()` option, controls written in factor notation other than the plainest 
  (eg, `ib2.cat`, `ibn.cat`, `io3.cat`, `i(2/5).cat`, `2.cat`) returned an empty table.
  They now estimate. Controls in simpler factor notation (eg, `i.cat`) were never affected.
- `rw` with `bootstrap()` and an `i.`-prefixed control entered through the `controls()` option 
  returned an empty table; it now returns estimates.
- `rw` with `bootstrap()` on a deterministic panel with no stochastic variation aborted the command with r(303); 
  it now returns results. This only affects deterministic panels, so it is irrelevant for any
  real-world application. It can be relevant for deterministic simulated test datasets without random variation.
- With `bootstrap()` option, a failure inside `boottest` aborted the entire run, discarding every other estimate.
  `boottest` can now decline to compute while still returning the other estimates.
- `pmd()` and pre-treatment pooled estimates could discard usable rows in presence of missing values. 
  In previous releases, the moving average was built with `egen … , filter()`, which drops a row whenever the outcome is missing
  at the *current* period — even when the current period is not in the window being averaged.
  `pmd(3)` averages the outcome at t−3, t−2 and t−1, so whether it is observed at t is
  irrelevant. Yet a missing value there removed the row. The window is now built arithmetically
  and the defect is gone. Post-treatment pooled estimates were never affected.
- The `e(dylags)` scalar, containing the number of first-differenced lags of the outcome used as covariates,
  mirroring `e(ylags)`, was promised by earlier versions but not actually produced. This was fixed in version 1.1.0.
- When `rw` uses a weighted regression (ie, when no covariates or additional fixed effects are added), 
  it could have some imperfection in presence of missing values. 
  The weights were computed in a sample that might not 
  be fully identical to the estimation sample in presence of missing values. 
  All this is now fixed, so the weights are always computed on the correct estimation sample also in presence of missing values.
- When `rw` uses a weighted regression (ie, when no covariates or additional fixed effects are added) and the user supplies their own weights, 
  the user-supplied weights were not incorporated into the re-weighting factor. Now they are (as they should).
- The level() option did not reach the pooled estimate's confidence interval on the Regression Adjustment specifications.
  On specifications that trigger regression adjustment, the pooled row's interval ignored
  `level()` and was always the 95% interval, while the event-study rows responded correctly. 
  Anyone who used `level()` with a value other than 95 on such a specification was given the 
  wrong pooled interval and should re-run. At `level(95)` nothing changes, and the point estimate, 
  standard error and p-value were never affected.

### Other changes

- **The test-statistic column is now labelled `z` on the regression-adjustment route.** With
  `rw` together with covariates and without `bootstrap()`, the p-value and confidence interval
  refer to the normal distribution, not to *t* — as they always have. The column in
  `e(results)` and `e(pooled_results)` was nonetheless named `t` on every route. It now reads
  `z` where the normal is used. Code that selects a column of those matrices *by name* is
  therefore route-dependent; code that selects by position is unaffected. No value changes.

- **Regression adjustment now exits when it can estimate nothing.** With `rw` and covariates
  and without `bootstrap()`, a run in which no horizon could be fitted previously returned
  r(0) together with an all-missing `e(results)`. It now exits with r(2001) and a message. 
  The default routes already behaved this way. 
  `rw` with `bootstrap()` is unchanged.

- **`egenmore` is no longer required.** The dependencies are `reghdfe`, `listreg` and
  `boottest`; `reghdfe` in turn needs `ftools` and `require`. If a previously working
  specification begins to fail after a Stata or package update, check that `reghdfe` runs on
  your data on its own.
  
- **Wild bootstrap now prints a message when repetitions are necessarily fewer than requested.** 
  When wild bootstrap is used (through the `bootstrap()` option), if clusters are too few to 
  perform the requested number of repetitions, the command now prints a note explaining that 
  `boottest` (the underlying command used for wild bootstrap inference) has enumerated the full 
  Rademacher universe, so the number of replications performed is necessarily smaller than requested.

---

## Version 1.0.3 (August 2026)

- **Bug fix:** the `pmd(max)` option could incorrectly compute the transformed outcome variable
  in panels with explicit missing values.
- **Bug fix:** using `pmd()` together with the `oneoff` suboption of `nonabsorbing()` imposed a
  slightly too restrictive definition of the clean control sample, possibly resulting in a loss
  of statistical power. This bug had been introduced in version 1.0.2.
- Improved numerical precision of the PMD and pooled outcome transformations, which are now
  computed in `double` precision rather than `float`.

---

## Version 1.0.2 (December 2025)

- The `rw` option with covariates or non-absorbing treatment runs much faster than in previous
  versions, and is now compatible with `bootstrap()` for wild bootstrap standard errors. This
  was achieved by (a) extending the set of cases in which a weighted regression is used,
  (b) using `listreg` to perform regression adjustment when wild cluster bootstrapping is not
  requested, and (c) using `margins` to perform regression adjustment with wild cluster
  bootstrap.
- Introduced the `oneoff` suboption within `nonabsorbing()`, for repeated one-off (or shock-type) 
  treatments (eg, hurricanes or other natural disasters).
- **Bug fix:** specifying both `pmd()` and `nonabsorbing()` gave a wrong sample definition,
  which resulted in most observations being dropped.
- **Bug fix:** specifying a control variable, `rw`, and a user-supplied weight together
  produced a syntax error.

---

## Version 1.0.1 (July 2024)

- Added the `absorb()` option, for additional absorbed fixed effects.
- Stata prefixes, including time-series operators, are now allowed in the arguments of
  `controls()`.
- Weights can now be supplied in the standard Stata format, `[pweight=varname]`.

---

## Version 1.0.0 (November 2023)

First release.
