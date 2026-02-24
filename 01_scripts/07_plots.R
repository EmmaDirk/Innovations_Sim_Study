# 07_plots.R
# -----------------------------------------------------------------
# Create one combined 2x4 figure from the *_int simulation outputs
#
# Expected input:
#   out <- simulation_wrapper(...)
#   out$summary has columns: aU, method, beta_hat, se_hat, bias
#   out$draws   has columns: aU, rep, method, beta_hat, se_hat, bias
#
# Figure:
#   Row 1) Bias vs aU (lines per method), faceted into 4 panels
#   Row 2) RMSE vs aU (lines per method), faceted into 4 panels
#
# RMSE is computed from replication draws:
#   RMSE(aU, method) = sqrt( mean( (beta_hat - b1)^2 ) )
# which equals:
#   RMSE(aU, method) = sqrt( mean( bias^2 ) )
#
# Notes:
#   - Uses ggplot2 (returns ggplot objects so ggsave() works).
#   - Uses patchwork to combine plots into a single 2x4 figure
#     with one shared legend at the bottom.
#   - Uses a shared viridis palette so method colours are identical
#     across the Bias and RMSE rows.
# -----------------------------------------------------------------

# load here package for file paths
library("here")

# step 0: load packages
source(here("01_scripts", "00_packages.R"))

# step 1: load the *_int simulation outputs
# These were saved in your master script as:
#   save(out_U_005_int, file = here("03_output", "out_U_005_int.RData"))
#   save(out_U_010_int, file = here("03_output", "out_U_010_int.RData"))
#   save(out_U_020_int, file = here("03_output", "out_U_020_int.RData"))
#   save(out_U_040_int, file = here("03_output", "out_U_040_int.RData"))
load(here("03_output", "out_U_005_int.RData"))  # object: out_U_005_int
load(here("03_output", "out_U_010_int.RData"))  # object: out_U_010_int
load(here("03_output", "out_U_020_int.RData"))  # object: out_U_020_int
load(here("03_output", "out_U_040_int.RData"))  # object: out_U_040_int

# -----------------------------------------------------------------
# helper: compute RMSE curves from replication-level draws
# RMSE = sqrt(mean(bias^2)) for each (aU, method)
# -----------------------------------------------------------------
compute_rmse_summary <- function(sim_out) {

  # checks
  if (is.null(sim_out$draws)) stop("sim_out must contain $draws.")
  d <- sim_out$draws

  # check columns
  need <- c("aU", "method", "bias")

  # if any missing, stop
  miss <- setdiff(need, names(d))
  if (length(miss) > 0) stop(paste("draws missing:", paste(miss, collapse = ", ")))

  # RMSE = sqrt(mean(bias^2)) for each (aU, method)
  rmse_df <- aggregate(
    bias ~ aU + method,
    data = d,
    FUN  = function(x) sqrt(mean(x^2, na.rm = TRUE))
  )

  # rename
  names(rmse_df)[names(rmse_df) == "bias"] <- "rmse"

  # return
  rmse_df
}

# -----------------------------------------------------------------
# step 2: build a combined data frame for Bias (from out$summary)
# This binds together the 4 *_int simulation outputs and adds a
# scenario label so we can facet into 4 columns.
# -----------------------------------------------------------------
bias_df <- bind_rows(
  out_U_005_int$summary %>%
    transmute(aU, method, bias, scenario = "U-confounding = 0.05"),
  out_U_010_int$summary %>%
    transmute(aU, method, bias, scenario = "U-confounding = 0.10"),
  out_U_020_int$summary %>%
    transmute(aU, method, bias, scenario = "U-confounding = 0.20"),
  out_U_040_int$summary %>%
    transmute(aU, method, bias, scenario = "U-confounding = 0.40")
) %>%
  mutate(
    scenario = factor(
      scenario,
      levels = c(
        "U-confounding = 0.05",
        "U-confounding = 0.10",
        "U-confounding = 0.20",
        "U-confounding = 0.40"
      )
    ),
    method = factor(method)
  )

# -----------------------------------------------------------------
# step 3: build a combined data frame for RMSE (from out$draws)
# RMSE is computed per (aU, method) within each *_int output, and
# then we bind the 4 scenarios together for faceting.
# -----------------------------------------------------------------
rmse_df <- bind_rows(
  compute_rmse_summary(out_U_005_int) %>%
    mutate(scenario = "U-confounding = 0.05"),
  compute_rmse_summary(out_U_010_int) %>%
    mutate(scenario = "U-confounding = 0.10"),
  compute_rmse_summary(out_U_020_int) %>%
    mutate(scenario = "U-confounding = 0.20"),
  compute_rmse_summary(out_U_040_int) %>%
    mutate(scenario = "U-confounding = 0.40")
) %>%
  mutate(
    scenario = factor(
      scenario,
      levels = c(
        "U-confounding = 0.05",
        "U-confounding = 0.10",
        "U-confounding = 0.20",
        "U-confounding = 0.40"
      )
    ),
    method = factor(method)
  )

# -----------------------------------------------------------------
# step 4: define a shared method colour palette
# Using one palette for both Bias and RMSE ensures methods have the
# same colours in both rows of the final figure.
# -----------------------------------------------------------------

# keep a fixed legend order so the lines match what we want to compare
want_methods <- c(
  "SRS",
  "NPS (200)",
  "NPS (10000)",
  "Combined (200)",
  "Combined (10000)"
)

# keep any extra methods (if present) after the main five
all_methods <- union(want_methods, union(levels(bias_df$method), levels(rmse_df$method)))

# enforce method ordering across both data frames
bias_df$method <- factor(bias_df$method, levels = all_methods)
rmse_df$method <- factor(rmse_df$method, levels = all_methods)

pal_method <- viridis(
  n      = length(all_methods),
  option = "viridis",
  begin  = 0.1,
  end    = 0.9
)
names(pal_method) <- all_methods

# -----------------------------------------------------------------
# step 5: plot 1 (row 1): Bias vs aU
# One line per method; dashed zero line; faceted into 4 panels.
# -----------------------------------------------------------------
p_bias_row <- ggplot(bias_df, aes(x = aU, y = bias, color = method, group = method)) +

  # add horizontal zero line
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +

  # draw one line per method
  geom_line(linewidth = 0.9) +

  # show points at each aU
  geom_point(size = 2) +

  # facet into the 4 U-confounding scenarios (4 columns)
  facet_wrap(~ scenario, nrow = 1) +

  # colours (keep as-is) + keep legend order stable across scenarios
  scale_color_manual(
    values = pal_method,
    breaks = all_methods
  ) +

  # labels
  labs(
    x = "aU (effect of U on selection)",
    y = "Bias",
    color = NULL
  ) +

  # scientific theme
  theme_classic(base_size = 13) +
  theme(
    panel.spacing.x = unit(1.2, "cm"),
    strip.text      = element_text(size = 14),
    axis.title      = element_text(size = 13),
    axis.text       = element_text(size = 12),
    legend.position = "bottom",
    legend.text     = element_text(size = 12)
  )

# -----------------------------------------------------------------
# step 6: plot 2 (row 2): RMSE vs aU
# RMSE is computed from replication draws; one line per method;
# faceted into 4 panels with the same layout as the Bias row.
# -----------------------------------------------------------------
p_rmse_row <- ggplot(rmse_df, aes(x = aU, y = rmse, color = method, group = method)) +

  # draw one line per method
  geom_line(linewidth = 0.9) +

  # show points at each aU
  geom_point(size = 2) +

  # facet into the 4 U-confounding scenarios (4 columns)
  facet_wrap(~ scenario, nrow = 1) +

  # colours (keep as-is) + keep legend order stable across scenarios
  scale_color_manual(
    values = pal_method,
    breaks = all_methods
  ) +

  # labels
  labs(
    x = "aU (effect of U on selection probability)",
    y = "RMSE",
    color = NULL
  ) +

  # scientific theme (match the Bias plot styling)
  theme_classic(base_size = 13) +
  theme(
    panel.spacing.x = unit(1.2, "cm"),
    strip.text      = element_text(size = 14),
    axis.title      = element_text(size = 13),
    axis.text       = element_text(size = 12),
    legend.position = "bottom",
    legend.text     = element_text(size = 12)
  )

# -----------------------------------------------------------------
# step 7: combine into one 2x4 figure
# - Bias row on top, RMSE row below
# - One shared legend at the bottom
# -----------------------------------------------------------------
p_big_2x4 <- (p_bias_row / p_rmse_row) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

p_big_2x4

# -----------------------------------------------------------------
# step 8: save the combined plot
# -----------------------------------------------------------------
ggsave(
  filename = here("03_output", "big_2x4_bias_rmse_int.png"),
  plot     = p_big_2x4,
  width    = 16,
  height   = 8,
  dpi      = 300
)