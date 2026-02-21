# -----------------------------------------------------------------
# 03_analysis.R
#
# this scripts has three analysis functions:
# 1. Analyse the SRS using a linear regression of Y on X and Z.
# 2. Analyse the NPS using a linear regression of Y on X, Z, and U.
# 3. Analyse the combined sample using a two step approach:
#     3a. Fit a logistic regression of S on X, Z, and Y to estimate the selection model.
#     3b. Fit a weighted linear regression of Y on X, Z, and U,
#         using inverse odds weights from the estimated selection model.
#
# Since we use two models to estimate the effect in the two step approach, we need to
# correct for the uncertainty in the first step when we compute standard errors for the second step.
# As such we will do the first step multiple times on different random draws from the
# combined sample, and then compute the standard error as the standard deviation of the estimates.
#
# IMPORTANT:
# - Function names are kept EXACTLY as in your current engine:
#   analyze_srs(), analyze_nps(), analyze_combined_ipw_bootstrap()
# - The combined weights are inverse-odds: (1 - p_hat)/p_hat
# - Small numeric guards + optional trimming are included for stability.
# -----------------------------------------------------------------

# 1. analyse the SRS using a linear regression of Y on X and Z.
analyze_srs <- function(srs_df) {

  # check that required variables exist
  need <- c("Y", "X", "Z")

  # missing variables if any
  miss <- setdiff(need, names(srs_df))

  # if any missing, stop
  if (length(miss) > 0) stop(paste("SRS missing:", paste(miss, collapse = ", ")))

  # fit OLS: Y on X and Z
  fit <- lm(Y ~ X + Z, data = srs_df)

  # pull coefficient table from summary
  sm <- summary(fit)$coefficients

  # Return beta for X and its standard error (plus the full model object)
  list(
    beta_X = unname(sm["X", "Estimate"]),      # estimated effect of X on Y
    se_X   = unname(sm["X", "Std. Error"]),    # standard error of beta_X
    fit    = fit                               # keep the fitted model in case you want more details
  )
}

# 2. analyse the NPS using a linear regression of Y on X, Z, and U.
analyze_nps <- function(nps_df) {

  # check that required variables exist
  need <- c("Y", "X", "Z", "U")

  # missing variables if any
  miss <- setdiff(need, names(nps_df))

  # if any missing, stop
  if (length(miss) > 0) stop(paste("NPS missing:", paste(miss, collapse = ", ")))

  # fit OLS: Y on X, Z, and U
  fit <- lm(Y ~ X + Z + U, data = nps_df)

  # pull coefficient table from summary
  sm <- summary(fit)$coefficients

  # Return beta for X and its standard error (plus the full model object)
  list(
    beta_X = unname(sm["X", "Estimate"]),      # estimated effect of X on Y
    se_X   = unname(sm["X", "Std. Error"]),    # standard error of beta_X
    fit    = fit
  )
}

# 3. Analyze the combined sample using a two step approach:
#     3a. Fit a model to estimate the selection model, S on X, Z, and Y.
#     3b. Fit a weighted linear regression of Y on X, Z, and U,
#         to estimate the effect of X on Y. 
# combine SRS and NPS using inverse-odds weights
# and compute bootstrap standard errors for beta_X
analyze_two_step <- function(srs_df,                     # dataframe for SRS
                             nps_df,                     # dataframe for NPS
                             B = 500,                    # number of bootstrap replications
                             ps_formula = ~ X + Z + Y,   # propensity score formula
                             trim_ps = c(0.01, 0.99),    # trim extreme p-hats for stability
                             trim_w  = c(0.01, 0.99),    # optional weight trimming by quantiles
                             stabilize = TRUE,
                             seed = NULL) {

  # if no seed is provided, set the seed
  if (!is.null(seed)) set.seed(seed)

  # checks
  # outcome model vars must exist in NPS
  need_nps <- c("Y","X","Z","U")

  # missing variables if any
  miss_nps <- setdiff(need_nps, names(nps_df))

  # if any missing, stop
  if (length(miss_nps) > 0) stop(paste("NPS missing:", paste(miss_nps, collapse=", ")))

  # propensity vars must exist in BOTH
  ps_vars <- all.vars(ps_formula)

  # missing variables if any
  miss_srs <- setdiff(ps_vars, names(srs_df))
  miss_nps2 <- setdiff(ps_vars, names(nps_df))

  # if any missing, stop
  if (length(miss_srs) > 0) stop(paste("SRS missing (for ps_formula):", paste(miss_srs, collapse=", ")))
  if (length(miss_nps2) > 0) stop(paste("NPS missing (for ps_formula):", paste(miss_nps2, collapse=", ")))

  # store sample sizes
  n_srs <- nrow(srs_df)
  n_nps <- nrow(nps_df)

  # container for bootstrap estimates
  beta_vec <- numeric(B)

  # bootstrap loop
  for (b in seq_len(B)) {

    # resample SRS with replacement
    idx_srs <- sample.int(n_srs, size = n_srs, replace = TRUE)
    srs_b <- srs_df[idx_srs, , drop = FALSE]

    # resample NPS with replacement
    idx_nps <- sample.int(n_nps, size = n_nps, replace = TRUE)
    nps_b <- nps_df[idx_nps, , drop = FALSE]

    # stack data for sample-membership model
    srs_ps <- srs_b[, ps_vars, drop = FALSE]
    nps_ps <- nps_b[, ps_vars, drop = FALSE]

    comb <- rbind(
      data.frame(srs_ps, source_nps = 0L),
      data.frame(nps_ps, source_nps = 1L)
    )

    # fit propensity model: P(source = NPS | shared vars)
    ps_fit <- glm(update(ps_formula, source_nps ~ .),
                  data = comb, family = binomial())

    # predict P(source = NPS | shared vars) for NPS bootstrap sample
    p_hat <- predict(ps_fit, newdata = nps_ps, type = "response")

    # trim p-hat away from 0/1
    p_lo <- trim_ps[1]; p_hi <- trim_ps[2]
    p_hat <- pmin(pmax(p_hat, p_lo), p_hi)

    # inverse-odds weights
    w <- (1 - p_hat) / p_hat

    # stabilize (keeps mean weight ~ 1-ish)
    if (stabilize) {
      pi_nps <- mean(comb$source_nps == 1L)
      pi_srs <- 1 - pi_nps
      w <- w * (pi_nps / pi_srs)
    }

    # optional trimming of extreme weights
    if (!is.null(trim_w)) {
      q <- quantile(w, probs = trim_w, na.rm = TRUE)
      w <- pmin(pmax(w, q[1]), q[2])
    }

    # normalize weights to sum to NPS sample size
    w <- w * (n_nps / sum(w))

    # weighted outcome model on bootstrap NPS
    fit_w <- lm(Y ~ X + Z + U, data = nps_b, weights = w)

    # store bootstrap estimate of beta_X
    beta_vec[b] <- coef(fit_w)["X"]
  }

  # return bootstrap mean and bootstrap SE
  list(
    beta_X = mean(beta_vec),         # bootstrap mean estimate
    se_X   = sd(beta_vec),           # bootstrap standard error
    boot_dist = beta_vec             # full bootstrap distribution
  )
}


