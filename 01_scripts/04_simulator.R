# -----------------------------------------------------------------
# Simulation function: bias curves vs aU for 3 approaches
#
# Saves BOTH raw estimates (beta_hat, se_hat) AND bias (beta_hat - b1)
#
# Inputs:
#   pop_df  : population data frame (must contain Y, X, Z, U)
#   b1      : true causal effect of X on Y (used to compute bias)
#   n_srs   : SRS sample size
#   n_nps   : NPS sample size
#   aX,aY,aZ: fixed selection effects for NPS
#   aU_vals : vector of effect of U on selection
#   R       : number of Monte Carlo replications per aU value
#   boot_B  : number of bootstrap draws used inside the combined 2-step estimator
#   base_seed      : base seed for reproducibility
#   show_progress  : show a progress bar (TRUE/FALSE)
#
# Returns:
#   list(summary=..., draws=..., settings=...)
#     - draws: replication-level rows with beta_hat, se_hat, bias
#     - summary: averages by aU and method (mean beta_hat, mean se_hat, mean bias)
# -----------------------------------------------------------------

simulate_bias_curves <- function(pop_df,                     # population data frame
                                 b1,                         # TRUE effect for bias calc
                                 n_srs,                      # sample size for SRS
                                 n_nps,                      # sample sizes for NPS
                                 aX,                         # effect of X on selection
                                 aY,                         # effect of Y on selection
                                 aZ,                         # effect of Z on selection
                                 aU_vals,                    # vector of effect of U on selection
                                 R = 200,                    # number of replications
                                 boot_B = 200,               # number of bootstrap draws for combined estimator
                                 base_seed = 1L,             # base seed for reproducibility
                                 show_progress = TRUE) {     # show progress bar

  # checks: population must have all variables needed
  need_pop <- c("Y", "X", "Z", "U")
  miss <- setdiff(need_pop, names(pop_df))
  if (length(miss) > 0) stop(paste("Population df missing:", paste(miss, collapse = ", ")))

  # check: b1 must be numeric scalar
  if (!is.numeric(b1) || length(b1) != 1 || !is.finite(b1)) {
    stop("b1 must be a single finite numeric value.")
  }

  # check: base_seed must be numeric scalar
  if (!is.numeric(base_seed) || length(base_seed) != 1 || !is.finite(base_seed)) {
    stop("base_seed must be a single finite numeric value.")
  }

  # convert base_seed to integer (if not already)
  base_seed <- as.integer(base_seed)

  # storage for replication-level results
  # we have length(aU_vals) * R replications, and 3 methods per replication
  out_list <- vector("list", length = length(aU_vals) * R * 3L)

  # counter for storing results in out_list
  k <- 0L

  # progress bar setup
  # total number of iterations
  total_iters <- length(aU_vals) * R

  # set up progress bar if requested
  if (show_progress) {

    # use txtProgressBar for a simple text-based progress bar
    pb <- txtProgressBar(min = 0, max = total_iters, style = 3)

    # counter for updating progress bar
    iter <- 0L

    # ensure progress bar is closed on exit (even if error occurs)
    on.exit(close(pb), add = TRUE)
  }

  # loop over aU values and replications
  for (i in seq_along(aU_vals)) {

    # current aU value for this iteration
    aU <- aU_vals[i]

    # for each replication within this aU value:
    for (r in seq_len(R)) {

      # 1) draw samples
      # base seeds + reproducible seeds for each replication
      seed_srs <- base_seed + 100000L * i + 10L * r + 1L
      seed_nps <- base_seed + 100000L * i + 10L * r + 2L

      # draw samples
      srs_df <- srs(pop_df, n = n_srs, seed = seed_srs)
      nps_df <- nps(pop_df, n = n_nps, aX = aX, aY = aY, aZ = aZ, aU = aU, seed = seed_nps)

      # 2) analyze SRS (Y ~ X + Z)
      res_srs <- analyze_srs(srs_df)

      # update counter and store results for SRS
      k <- k + 1L
      out_list[[k]] <- data.frame(
        aU = aU,
        rep = r,
        method = "SRS: OLS(Y~X+Z)",
        beta_hat = res_srs$beta_X,           # raw estimate
        se_hat = res_srs$se_X,               # raw SE
        bias = res_srs$beta_X - b1           # bias vs provided b1
      )

      # 3) analyze NPS (Y ~ X + Z + U)
      res_nps <- analyze_nps(nps_df)

      # update counter and store results for NPS
      k <- k + 1L
      out_list[[k]] <- data.frame(
        aU = aU,
        rep = r,
        method = "NPS: OLS(Y~X+Z+U)",
        beta_hat = res_nps$beta_X,
        se_hat = res_nps$se_X,
        bias = res_nps$beta_X - b1
      )

      # 4) combined 2-step with bootstrap SE (selection on X,Z; outcome on NPS with U)
      # reproducible seeds for bootstrap (different from sample draws)
      seed_boot <- base_seed + 100000L * i + 10L * r + 3L

      # analyze combined sample with bootstrap SE
      res_ipw <- analyze_combined_ipw_bootstrap(srs_df, nps_df, B = boot_B, seed = seed_boot)

      # update counter and store results for combined estimator
      k <- k + 1L
      out_list[[k]] <- data.frame(
        aU = aU,
        rep = r,
        method = "Combined: IPW + U (boot SE)",
        beta_hat = res_ipw$beta_X,
        se_hat = res_ipw$se_boot,           
        bias = res_ipw$beta_X - b1
      )

      # progress update
      if (show_progress) {
        iter <- iter + 1L
        setTxtProgressBar(pb, iter)
      }
    }
  }

  # convert list to data frame
  draws <- do.call(rbind, out_list)

  # summarize for plotting (mean lines)
  summary <- aggregate(
    cbind(beta_hat, se_hat, bias) ~ aU + method,
    data = draws,
    FUN = mean,
    na.rm = TRUE
  )

  # save results
  list(
    summary = summary,
    draws = draws,
    settings = list(
      b1 = b1,
      n_srs = n_srs,
      n_nps = n_nps,
      aX = aX, aY = aY, aZ = aZ,
      aU_vals = aU_vals,
      R = R,
      boot_B = boot_B,
      base_seed = base_seed
    )
  )
}
