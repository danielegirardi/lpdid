#
# This R script applies the Local Projections Difference-in-Differences (LP-DiD) estimator in a simulated example dataset.
# This example focuses on the reweighted version of LP-DiD, which yields the equally-weighted average treatment effect on the treated
# This script demonstrates two numerically equivalent ways to obtain the equally-weighted ATT using LP-DiD: 
# the 1st computes the variance weighting inherent in OLS, and then "undoes" them through a weighted regression, forcing equal weights on all treated units
# the 2nd uses a Regression-Adjustment LP-DiD specification, which requires fewer lines of code but can be slower to run
# (other scripts in this repository show how to implement variance-weighted LP-DiD, which yields a variance-weighted ATT with strictly positive weights)
# Related paper: Dube, Girardi, Jorda and Taylor (2025) "A Local Projections Approach to Difference-in-Differences", Journal of Applied Econometrics 40 (7), https://doi.org/10.1002/jae.70000
# This and other example files can be downloaded at https://github.com/danielegirardi/lpdid/ 
# Author: Daniele Girardi (King's College London), daniele.girardi@kcl.ac.uk
# 23 Feb 2026
# 

# Required packages
library(haven)            # loading .dta files
library(dplyr)            # data manipulation
library(fixest)           # regressions with FEs and clustered SEs
library(marginaleffects)  # regression adjustment estimator
library(tsibble)          # time-aware panel operations 
library(slider)           # moving averages (for "pooled" estimates)
library(ggplot2)          # plotting

# Match reghdfe small-sample correction for clustered SEs
setFixest_ssc(ssc(adj = TRUE, cluster.adj = TRUE))


# ============================================================
#
#   (1) Upload dataset and preparation
#
# ============================================================

# Upload simulated dataset with staggered binary treatment
# (treatment is absorbing: once a unit is treated, it stays treated forever)
# (other example files in this repository demonstrate the use of LP-DiD with non-absorbing treatment, where units can enter and exit treatment multiple times.)
df <- read_dta("http://fmwww.bc.edu/repec/bocode/l/lpdidtestdata1.dta")

# Set estimation window
post_window <- 5
pre_window  <- 5

# Convert to a tsibble (time-aware panel structure) and fill any gaps.
# (filling gaps makes the code compatible with unbalanced panels - although this specific simulated panel dataset is balanced).
df_ts <- df |>
  as_tsibble(key = unit, index = time) |>
  fill_gaps()

# Compute differenced treatment indicator (D_treat) and lag of outcome (L_Y)
df_ts <- df_ts |>
  group_by(unit) |>
  arrange(time, .by_group = TRUE) |>
  mutate(
    # First difference of treatment 
    D_treat = treat - dplyr::lag(treat, 1),
    # Lag of outcome 
    L_Y     = dplyr::lag(Y, 1)
  ) |>
  ungroup()

# Dynamically generate leads of treatment (to be used in the clean control condition) 
# and long differences of the outcome (to be used as the left-hand variables in event study estimates). 
# Using loops so the number of leads/lags automatically follows post_window / pre_window.
df_ts <- df_ts |>
  group_by(unit) |>
  arrange(time, .by_group = TRUE)

for (h in 0:post_window) {
  df_ts <- df_ts |>
    mutate(
      # Lead of treatment for clean control condition at horizon h
      !!paste0("F", h, "_treat") := dplyr::lead(treat, h),
      # long forward difference: D{h}y = F{h}.Y - L.Y
      !!paste0("D", h, "y")      := dplyr::lead(Y, h) - L_Y
    )
}

for (h in 2:pre_window) {
  df_ts <- df_ts |>
    mutate(
      # long pre-period difference: Dm{h}y = L{h}.Y - L.Y
      !!paste0("Dm", h, "y") := dplyr::lag(Y, h) - L_Y
    )
}

# Create pooled outcome variables (to be used in pooled estimates)
# slide_dbl computes a moving average within each unit group.
# .before / .after define how many rows before / after the current row to include.
# setting na.rm=FALSE and .complete=TRUE ensures that averages are not computed if there are missing values in the window.
df_ts <- df_ts |>
  mutate(
    # Post pooled: average of Y from current period to post_window ahead, minus Y at t-1
    pooled_post_y = slider::slide_dbl(Y, \(x) mean(x, na.rm = FALSE),
                                      .before = 0, .after = post_window, .complete = TRUE) - L_Y,
    # Pre pooled: average of Y from pre_window to 2 periods back, minus Y at t-1
    pooled_pre_y  = slider::slide_dbl(Y, \(x) mean(x, na.rm = FALSE),
                                      .before = pre_window, .after = -2, .complete = TRUE) - L_Y
  ) |>
  ungroup()

# Drop the artificially inserted gap rows (if there are any) and convert back to a plain data frame
df <- df_ts |>
  dplyr::filter(!is.na(Y)) |>
  as_tibble()

# Compute and store the true (equally-weighted) dynamic average treatment effects from the simulation
# (this will not be possible with real-world data!)
true_fx <- df |>
  dplyr::filter(!is.na(event_time), event_time >= 0, event_time <= post_window) |>
  group_by(event_time) |>
  summarise(
    true_att     = mean(effect, na.rm = TRUE),
    min_true_att = min(effect,  na.rm = TRUE),
    max_true_att = max(effect,  na.rm = TRUE),
    .groups = "drop"
  ) |>
  rename(horizon = event_time)

# Create results data frame indexed by event-time horizon
horizons <- (-pre_window):post_window   # -5, -4, ..., 0, ..., 5
results  <- data.frame(horizon = horizons) |>
  left_join(true_fx, by = "horizon") |>
  mutate(
    true_att     = replace(true_att,     horizon < 0, 0),
    min_true_att = replace(min_true_att, horizon < 0, 0),
    max_true_att = replace(max_true_att, horizon < 0, 0)
  )


# =============================================================================================
#
#   (2) Reweighted LP-DiD event study estimates - weighted regression implementation
#
# =============================================================================================

# Define a function that, at each given time horizon, computes the variance-based weights inherent in OLS and takes their inverse 
# (in subsequent regressions, we will weight observations by the inverse of the OLS variance-based weights, 
# thus "undoing" the variance-weighting inherent in OLS and yielding equal weights)
get_reweights <- function(df, h) {
  lead_col <- paste0("F", h, "_treat")
  clean_h  <- df |>
    dplyr::filter(!is.na(D_treat), D_treat == 1 | .data[[lead_col]] == 0)

  # Regress D_treat on time FEs within the clean sample; extract residuals
  m_w <- feols(D_treat ~ 1 | time, data = clean_h, vcov = "iid")
  clean_h$num_weights <- residuals(m_w)

  # Only keep the residuals for the treated observations
  clean_h <- clean_h |>
    mutate(num_weights = if_else(D_treat != 1, NA_real_, num_weights))

  # Compute the OLS variance-based weights = residual / sum-of-treated-residuals
  den <- sum(clean_h$num_weights, na.rm = TRUE)
  clean_h <- clean_h |>
    mutate(weight = num_weights / den) |>
    group_by(time) |>
    mutate(
      # Assign each control unit the weight of the corresponding treated units
      gweight = suppressWarnings(max(weight, na.rm = TRUE)),
      gweight = if_else(is.infinite(gweight), NA_real_, gweight),
      weight  = if_else(is.na(weight), gweight, weight)
    ) |>
    ungroup() |>
    mutate(reweight = 1 / weight)

  clean_h[, c("unit", "time", "reweight")]
}

# Compute and store weights for all post-treatment time horizons (h = 0, ..., post_window)
reweights <- lapply(0:post_window, function(h) get_reweights(df, h))

# Create empty vectors where to store the event study estimates
b_rw_lpdid  <- setNames(rep(NA_real_, length(horizons)), as.character(horizons))
se_rw_lpdid <- setNames(rep(NA_real_, length(horizons)), as.character(horizons))

# Post periods: merge horizon-specific reweights into clean sample, then run weighted feols
# (note: R starts counting at 1, not 0, so the weights for h=0 are stored in reweights[[1]], so we use reweights[[h+1]] below here)
for (h in 0:post_window) {
  lead_col <- paste0("F", h, "_treat")
  clean_h  <- df |>
    dplyr::filter(!is.na(D_treat), D_treat == 1 | .data[[lead_col]] == 0) |>
    left_join(reweights[[h + 1]], by = c("unit", "time"))
  m <- feols(as.formula(paste0("D", h, "y ~ D_treat | time")),
             data = clean_h, weights = ~reweight, vcov = ~unit)
  b_rw_lpdid[as.character(h)]  <- coef(m)["D_treat"]
  se_rw_lpdid[as.character(h)] <- se(m)["D_treat"]
}

# Pre periods: use reweight_0 (absorbing treatment: pre-period weights = h=0 weights)
clean_pre_rw <- df |>
  dplyr::filter(!is.na(D_treat), D_treat == 1 | treat == 0) |>
  left_join(reweights[[1]], by = c("unit", "time"))   # reweights[[1]] = h=0

for (h in 2:pre_window) {
  m <- feols(as.formula(paste0("Dm", h, "y ~ D_treat | time")),
             data = clean_pre_rw, weights = ~reweight, vcov = ~unit)
  b_rw_lpdid[as.character(-h)]  <- coef(m)["D_treat"]
  se_rw_lpdid[as.character(-h)] <- se(m)["D_treat"]
}

b_rw_lpdid["-1"]  <- 0
se_rw_lpdid["-1"] <- NA

results$b_rw_lpdid  <- b_rw_lpdid[as.character(horizons)]
results$se_rw_lpdid <- se_rw_lpdid[as.character(horizons)]

# =============================================================================================
#
#   (3) Reweighted LP-DiD event study estimates - Regression Adjustment (RA) implementation
#
# =============================================================================================

# This is an alternative and numerically equivalent way to implement reweighted LP-DiD, using RA instead of the weighted regression.
# Advantage: you don't need to explicitly compute the variance-based weights, so the code is much shorter and simpler.
# Disadvantage: because RA can be computationally heavier, it can be slower to run.

b_lpdid_ra  <- setNames(rep(NA_real_, length(horizons)), as.character(horizons))
se_lpdid_ra <- setNames(rep(NA_real_, length(horizons)), as.character(horizons))

# Post periods
for (h in 0:post_window) {
  message("Estimating regression adjustment at horizon ", h)
  dep_var  <- paste0("D", h, "y")
  lead_col <- paste0("F", h, "_treat")
  clean_h  <- df |>
    dplyr::filter(!is.na(D_treat), D_treat == 1 | .data[[lead_col]] == 0) |>
    mutate(D_treat_f = factor(D_treat))
  # D_treat_f * factor(time) allows separate time FEs per treatment group,
  # equivalent to fitting separate regressions (as teffects ra or margins do in Stata)
  m_ra  <- lm(as.formula(paste0(dep_var, " ~ D_treat_f * factor(time)")), data = clean_h)
  atet  <- avg_comparisons(m_ra, variables = "D_treat_f",
                           newdata = dplyr::filter(clean_h, D_treat == 1),
                           vcov    = ~unit)
  b_lpdid_ra[as.character(h)]  <- atet$estimate[1]
  se_lpdid_ra[as.character(h)] <- atet$std.error[1]
}

# Pre periods
for (h in 2:pre_window) {
  message("Estimating regression adjustment at horizon minus ", h)
  dep_var <- paste0("Dm", h, "y")
  clean_h <- df |>
    dplyr::filter(!is.na(D_treat), D_treat == 1 | treat == 0) |>
    mutate(D_treat_f = factor(D_treat))
  m_ra  <- lm(as.formula(paste0(dep_var, " ~ D_treat_f * factor(time)")), data = clean_h)
  atet  <- avg_comparisons(m_ra, variables = "D_treat_f",
                           newdata = dplyr::filter(clean_h, D_treat == 1),
                           vcov    = ~unit)
  b_lpdid_ra[as.character(-h)]  <- atet$estimate[1]
  se_lpdid_ra[as.character(-h)] <- atet$std.error[1]
}

b_lpdid_ra["-1"]  <- 0
se_lpdid_ra["-1"] <- NA

results$b_lpdid_ra  <- b_lpdid_ra[as.character(horizons)]
results$se_lpdid_ra <- se_lpdid_ra[as.character(horizons)]


# ===========================================================================
#
#   (4) Pooled LP-DiD estimates - Weighted regression implementation
#
# ===========================================================================

# ---- Reweighted pooled LP-DiD estimates ----

# Post: clean control condition is the most restrictive (clean through end of post window)
post_pool_clean <- df |> dplyr::filter(!is.na(D_treat), D_treat == 1 | .data[[paste0("F", post_window, "_treat")]] == 0)

# ---- Weighted regression method ----

# Post: use weights for h=post_window 
clean_rw_post_pool <- post_pool_clean |>
  left_join(reweights[[post_window + 1]], by = c("unit", "time"))
m <- feols(pooled_post_y ~ D_treat | time,
           data = clean_rw_post_pool, weights = ~reweight, vcov = ~unit)
b_rw_post_pool  <- coef(m)["D_treat"]
se_rw_post_pool <- se(m)["D_treat"]

# Pre: same clean control condition as pre event-study
m <- feols(pooled_pre_y ~ D_treat | time,
           data = clean_pre_rw, weights = ~reweight, vcov = ~unit)
b_rw_pre_pool  <- coef(m)["D_treat"]
se_rw_pre_pool <- se(m)["D_treat"]

# ===========================================================================
#
#   (5) Pooled LP-DiD estimates - Regression Adjustment (RA) implementation
#
# ===========================================================================

# ---- Regression adjustment (RA) method ----

# Post
clean_ra_post <- post_pool_clean |> mutate(D_treat_f = factor(D_treat))
m_ra <- lm(pooled_post_y ~ D_treat_f * factor(time), data = clean_ra_post)
atet <- avg_comparisons(m_ra, variables = "D_treat_f",
                        newdata = dplyr::filter(clean_ra_post, D_treat == 1),
                        vcov    = ~unit)
b_ra_post_pool  <- atet$estimate[1]
se_ra_post_pool <- atet$std.error[1]

# Pre
pre_pool_clean <- df |> dplyr::filter(!is.na(D_treat), D_treat == 1 | treat == 0)
clean_ra_pre <- pre_pool_clean |> mutate(D_treat_f = factor(D_treat))
m_ra <- lm(pooled_pre_y ~ D_treat_f * factor(time), data = clean_ra_pre)
atet <- avg_comparisons(m_ra, variables = "D_treat_f",
                        newdata = dplyr::filter(clean_ra_pre, D_treat == 1),
                        vcov    = ~unit)
b_ra_pre_pool  <- atet$estimate[1]
se_ra_pre_pool <- atet$std.error[1]

# Compute and store the true pooled average effect (true effect averaged over the post-treatment window)
# (this will not be possible with real-world data!)
true_pooled_post <- mean(df$effect[!is.na(df$event_time) &
                                     df$event_time >= 0 &
                                     df$event_time <= post_window], na.rm = TRUE)

# Collect pooled results into a table
pooled_results <- data.frame(
  `Post-Coef` = c(true_pooled_post, b_rw_post_pool, b_ra_post_pool),
  `Post-SE`   = c(NA,               se_rw_post_pool, se_ra_post_pool),
  `Pre-Coef`  = c(0,                b_rw_pre_pool,  b_ra_pre_pool),
  `Pre-SE`    = c(NA,               se_rw_pre_pool, se_ra_pre_pool),
  row.names   = c("True_ATT", "RW_LP-DiD", "RW_LP-DiD_RA"),
  check.names = FALSE
)


# ============================================================
#
#   (6) Display results
#
# ============================================================

# ---- Event study estimates table ----

cat("Reweighted LP-DiD Event Study Estimates\n")
es_table <- results[, c("horizon", "true_att",
                         "b_rw_lpdid", "se_rw_lpdid",
                         "b_lpdid_ra", "se_lpdid_ra")]
print(data.frame(lapply(es_table, function(x) if (is.numeric(x)) round(x, 2) else x)),
      row.names = FALSE)

# ---- Event study graph ----

ggplot(results, aes(x = horizon)) +
  geom_ribbon(aes(ymin = min_true_att, ymax = max_true_att, fill = "Range of True Treatment Effects"),
              alpha = 0.2) +
  geom_line(aes(y = true_att,   color = "True Equally Weighted ATT"),
            linewidth = 1.0) +
  geom_errorbar(aes(ymin  = b_rw_lpdid - 1.96 * se_rw_lpdid,
                    ymax  = b_rw_lpdid + 1.96 * se_rw_lpdid,
                    color = "Reweighted LP-DiD Estimate"),
                width = 0.3) +
  geom_point(aes(y = b_rw_lpdid, color = "Reweighted LP-DiD Estimate"),
             size = 1.5) +
  scale_fill_manual(name = NULL, values = c("Range of True Treatment Effects" = "blue")) +
  scale_color_manual(name = NULL, values = c(
    "True Equally Weighted ATT" = "gray50",
    "Reweighted LP-DiD Estimate"         = "orange"
  )) +
  labs(
    x     = "Event Time",
    y     = "Treatment Effect",
    title = "Actual and Estimated Treatment Effects"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom", legend.box = "vertical")

# ---- Pooled estimates table ----

cat("\nReweighted LP-DiD Pooled Estimates\n")
print(round(pooled_results, 2))
