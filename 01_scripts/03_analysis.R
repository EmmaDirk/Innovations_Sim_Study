# this scripts has three analysis functions:
# 1. Analyse the SRS using a linear regression of Y on X and Z. 
# 2. Analyse the NPS using a linear regression of Y on X, Z, and U.
# 3. Analyse the combined sample using a two step approach:
#     3a. Fit a logistic regression of S on X, Z to estimate the selection model.
#     3b. Fit a weighted linear regression of Y on X and Z, 
#         using inverse probability weights from the estimated selection model.
#
# Since we use two models to estimate the effect in the two step approach, we need to 
# correct for the uncertainty in the first step when we compute standard errors for the second step.
# As such we will do the first step multiple times on different random draws from the 
# combined sample, and then compute the standard error as the standard deviation of the estimates.
# ------------------------------------------------------------------------------

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

# 3. analyse the combined sample using a two step approach
analyze_combined_ipw_bootstrap <- function(srs_df,                      # data frame of SRS
                                          nps_df,                       # data frame of NPS
                                          B = 500,                      # number of bootstrap samples
                                          seed = NULL) {

  # if no seed is provided, set the seed
  if (!is.null(seed)) set.seed(seed)

  # checks for selection model (needs X and Z in both samples)
  need_sel <- c("X", "Z")
  for (v in need_sel) {
    if (!(v %in% names(srs_df))) stop(paste("SRS missing:", v))
    if (!(v %in% names(nps_df))) stop(paste("NPS missing:", v))
  }

  # checks for outcome model (needs Y, X, Z, U in NPS)
  need_out <- c("Y", "X", "Z", "U")
  for (v in need_out) {
    if (!(v %in% names(nps_df))) stop(paste("NPS missing:", v))
  }

  # stack the data, create the selection variable
  # (note: we only need X and Z for the selection model; Y is not needed in SRS)
  srs_df$S <- 0
  nps_df$S <- 1
  comb <- rbind(
    srs_df[, c("X", "Z", "S")],
    nps_df[, c("X", "Z", "S")]
  )

  # total sample sizes
  n_srs <- nrow(srs_df)
  n_nps <- nrow(nps_df)

  # fit the two step model
  two_step_fit <- function(srs_d, nps_d) {

    # add selection indicator (in case bootstrap removed it / for safety)
    srs_d$S <- 0
    nps_d$S <- 1

    # stack the data for the selection model (needs only X, Z, S)
    d_sel <- rbind(
      srs_d[, c("X", "Z", "S")],
      nps_d[, c("X", "Z", "S")]
    )

    # 1) selection model to predict S from X and Z (U is unobserved in SRS)
    sel_fit <- glm(S ~ X + Z, data = d_sel, family = binomial())

    # 2) predicted probabilities of selection for NPS units
    p_hat_nps <- as.numeric(predict(sel_fit, newdata = nps_d, type = "response"))

    # avoid exactly 0 or 1 probabilities (prevents infinite weights)
    eps <- 1e-6
    p_hat_nps <- pmin(pmax(p_hat_nps, eps), 1 - eps)

    # IPW weights for NPS: weight = 1 / p_hat
    w_nps <- 1 / p_hat_nps

    # 3) weighted outcome model on NPS 
    out_fit <- lm(Y ~ X + Z + U, data = nps_d, weights = w_nps)
    sm <- summary(out_fit)$coefficients

    # save beta_X and se_X
    # note that these are naive SEs that ignore the uncertainty in the selection model, 
    # which is why we will do the bootstrap to get corrected SEs
    list(
      beta_X = unname(sm["X", "Estimate"]),
      se_naive = unname(sm["X", "Std. Error"]), 
      selection_model = sel_fit,
      outcome_model = out_fit
    )
  }

  # get the point estimate once 
  point <- two_step_fit(srs_df, nps_df)

  # save the point estimate
  beta_hat <- point$beta_X

  # get the bootstrap SEs
  beta_boot <- numeric(B)

  # for b in the number of bootstrap samples
  for (b in seq_len(B)) {

    # bootstrap within each sample (keeps the two-sample structure intact)
    srs_b <- srs_df[sample.int(n_srs, size = n_srs, replace = TRUE), , drop = FALSE]
    nps_b <- nps_df[sample.int(n_nps, size = n_nps, replace = TRUE), , drop = FALSE]

    # If bootstrap draw has only one class of S in the stacked data, glm can't run -> mark NA
    # (rare, but possible if one sample size is tiny)
    d_sel_b <- rbind(
      transform(srs_b[, c("X", "Z")], S = 0),
      transform(nps_b[, c("X", "Z")], S = 1)
    )

    if (length(unique(d_sel_b$S)) < 2) {
      beta_boot[b] <- NA_real_
    } else {
      beta_boot[b] <- tryCatch(two_step_fit(srs_b, nps_b)$beta_X, error = function(e) NA_real_)
    }
  }

  # if any bootstrap estimates are NA or infinite, we remove them before computing SE and CI
  beta_boot_ok <- beta_boot[is.finite(beta_boot)]

  # compute bootstrap SE as the standard deviation of the bootstrap estimates
  se_boot <- sd(beta_boot_ok)

  # compute percentile 95% CI
  ci_95 <- as.numeric(quantile(beta_boot_ok, probs = c(0.025, 0.975), na.rm = TRUE))

  # keep the results
  list(
    beta_X = beta_hat,              # 2-step point estimate
    se_boot = se_boot,              # bootstrap SE (accounts for 2-step uncertainty)
    ci_95 = ci_95,                  # percentile 95% CI
    se_naive = point$se_naive,      # naive SE from weighted lm only (ignores step 1)
    B = B,
    n_boot_ok = length(beta_boot_ok),
    beta_boot = beta_boot_ok,
    selection_model = point$selection_model,
    outcome_model = point$outcome_model
  )
}
