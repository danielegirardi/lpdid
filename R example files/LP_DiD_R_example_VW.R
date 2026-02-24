#
# This R script applies the Local Projections Difference-in-Differences (LP-DiD) estimator in a simulated example dataset.
# Related paper: Dube, Girardi, Jorda and Taylor (2025) "A Local Projections Approach to Difference-in-Differences", Journal of Applied Econometrics 40 (7), https://doi.org/10.1002/jae.70000
# This and other example files can be downloaded at https://github.com/danielegirardi/lpdid/ 
# Author: Daniele Girardi (King's College London), daniele.girardi@kcl.ac.uk
# 23 Feb 2026
# 

# Required packages
library(haven)    # loading .dta files
library(dplyr)    # data manipulation
library(fixest)   # regressions with FEs and clustered SEs
library(tsibble)  # time-aware panel operations
library(slider)   # moving averages (for "pooled" estimates)
library(ggplot2)  # plotting

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

# Convert to a tsibble (time-aware panel structure) and fill any gaps
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
# and long differences of the outcome (to be used as the left-hand variables in event study estimates)
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
# Setting na.rm to FALSE and .complete to TRUE ensure that averages are not computed if there are missing values in the window.
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


# ============================================================
#
#   (2) LP-DiD event-study estimates
#
# ============================================================

# ---- Baseline (variance-weighted) LP-DiD event study estimates ----
# (other scripts in this repository show how to implement re-weighted LP-DiD for the equally-weighted ATT)

b_lpdid_vw  <- setNames(rep(NA_real_, length(horizons)), as.character(horizons))
se_lpdid_vw <- setNames(rep(NA_real_, length(horizons)), as.character(horizons))

# Estimate dynamic treatment effects in the post-treatment window
for (h in 0:post_window) {
  lead_col <- paste0("F", h, "_treat")
  clean_h  <- df |> dplyr::filter(!is.na(D_treat), D_treat == 1 | .data[[lead_col]] == 0)
  m <- feols(as.formula(paste0("D", h, "y ~ D_treat | time")), data = clean_h, vcov = ~unit)
  b_lpdid_vw[as.character(h)]  <- coef(m)["D_treat"]
  se_lpdid_vw[as.character(h)] <- se(m)["D_treat"]
}

# Estimate pre-trends in the pre-treatment periods
for (h in 2:pre_window) {
  clean_h <- df |> dplyr::filter(!is.na(D_treat), D_treat == 1 | treat == 0)
  m <- feols(as.formula(paste0("Dm", h, "y ~ D_treat | time")), data = clean_h, vcov = ~unit)
  b_lpdid_vw[as.character(-h)]  <- coef(m)["D_treat"]
  se_lpdid_vw[as.character(-h)] <- se(m)["D_treat"]
}

b_lpdid_vw["-1"]  <- 0    # reference period
se_lpdid_vw["-1"] <- NA

results$b_lpdid_vw  <- b_lpdid_vw[as.character(horizons)]
results$se_lpdid_vw <- se_lpdid_vw[as.character(horizons)]


# ============================================================
#
#   (3) Pooled LP-DiD estimates
#
# ============================================================

# ---- Baseline (variance-weighted) pooled LP-DiD estimates----
# (other scripts in this repository show how to implement re-weighted LP-DiD for the equally-weighted ATT)

# Post: clean control condition is the most restrictive (clean through end of post window)
post_pool_clean <- df |> dplyr::filter(!is.na(D_treat), D_treat == 1 | .data[[paste0("F", post_window, "_treat")]] == 0)
m <- feols(pooled_post_y ~ D_treat | time, data = post_pool_clean, vcov = ~unit)
b_vw_post_pool  <- coef(m)["D_treat"]
se_vw_post_pool <- se(m)["D_treat"]

# Pre: same clean control condition as pre event-study
pre_pool_clean <- df |> dplyr::filter(!is.na(D_treat), D_treat == 1 | treat == 0)
m <- feols(pooled_pre_y ~ D_treat | time, data = pre_pool_clean, vcov = ~unit)
b_vw_pre_pool  <- coef(m)["D_treat"]
se_vw_pre_pool <- se(m)["D_treat"]

# Compute and store the true pooled average effect (true effect averaged over the post-treatment window)
# (this will not be possible with real-world data!)
true_pooled_post <- mean(df$effect[!is.na(df$event_time) &
                                     df$event_time >= 0 &
                                     df$event_time <= post_window], na.rm = TRUE)

# Collect pooled results into a table
pooled_results <- data.frame(
  `Post-Coef` = c(true_pooled_post, b_vw_post_pool),
  `Post-SE`   = c(NA,               se_vw_post_pool),
  `Pre-Coef`  = c(0,                b_vw_pre_pool),
  `Pre-SE`    = c(NA,               se_vw_pre_pool),
  row.names   = c("True_ATT", "VW_LP-DiD"),
  check.names = FALSE
)


# ============================================================
#
#   (4) Display results
#
# ============================================================

# ---- Event study estimates table ----

cat("LP-DiD Event Study Estimates (Variance-Weighted)\n")
es_table <- results[, c("horizon", "true_att", "b_lpdid_vw", "se_lpdid_vw")]
print(data.frame(lapply(es_table, function(x) if (is.numeric(x)) round(x, 2) else x)),
      row.names = FALSE)

# ---- Event study graph ----

ggplot(results, aes(x = horizon)) +
  geom_ribbon(aes(ymin = min_true_att, ymax = max_true_att, fill = "Range of True Treatment Effects"),
              alpha = 0.2) +
  geom_line(aes(y = true_att,   color = "True Equally Weighted ATE"),
            linewidth = 1.0) +
  geom_errorbar(aes(ymin  = b_lpdid_vw - 1.96 * se_lpdid_vw,
                    ymax  = b_lpdid_vw + 1.96 * se_lpdid_vw,
                    color = "LP-DiD Estimate"),
                width = 0.3) +
  geom_point(aes(y = b_lpdid_vw, color = "LP-DiD Estimate"),
             size = 1.5) +
  scale_fill_manual(name = NULL, values = c("Range of True Treatment Effects" = "blue")) +
  scale_color_manual(name = NULL, values = c(
    "True Equally Weighted ATE" = "gray50",
    "LP-DiD Estimate"  = "darkgreen"
  )) +
  labs(
    x     = "Event Time",
    y     = "Treatment Effect",
    title = "Actual and Estimated Treatment Effects"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom", legend.box = "vertical")

# ---- Pooled estimates table ----

cat("\nLP-DiD Pooled Estimates (Variance-Weighted)\n")
print(round(pooled_results, 2))
