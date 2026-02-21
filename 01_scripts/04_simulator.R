# -----------------------------------------------------------------
# Simulation function: bias curves vs aU for 3 approaches
#
# Saves BOTH raw estimates (beta_hat, se_hat) AND bias (beta_hat - b1)
#
# UPDATED to also save the nonlinear term b6 (coefficient on X2c = X^2 - 1):
#   - beta_hat_b1, se_hat_b1, bias_b1  (for X)
#   - beta_hat_b6, se_hat_b6, bias_b6  (for X2c = X^2 - 1)
#
# Inputs:
#   pop_df  : population data frame (must contain Y, X, Z, U)
#   b1      : true linear coefficient on X (used to compute bias for b1)
#   b6      : true coefficient on (X^2 - 1) (used to compute bias for b6)
#   n_srs   : SRS sample size
#   n_nps   : NPS sample size
#   aX,aY,aZ: fixed selection effects for NPS
#   aU_vals : vector of effect of U on selection
#   R       : number of Monte Carlo replications per aU value
#   boot_B  : number of bootstrap draws used inside the combined 2-step estimator
#   base_seed      : base seed for reproducibility
#   show_progress  : show a progress bar (TRUE/FALSE)
#   n_cores        : number of parallel workers (NULL uses detectCores()-1)
#
# Returns:
#   list(summary=..., draws=..., settings=...)
#     - draws: replication-level rows with beta_hat, se_hat, bias
#     - summary: averages by aU and method (mean beta_hat, mean se_hat, mean bias)
# -----------------------------------------------------------------

simulation_wrapper <- function(pop_df,                       # population data frame
                               b1,                           # TRUE linear effect for bias calc
                               b6,                           # TRUE nonlinear effect for bias calc (X^2 - 1)
                               n_srs,                        # sample size for SRS
                               n_nps,                        # sample sizes for NPS
                               aX,                           # effect of X on selection
                               aY,                           # effect of Y on selection
                               aZ,                           # effect of Z on selection
                               aU_vals,                      # vector of effect of U on selection
                               R = 200,                      # number of replications
                               boot_B = 200,                 # number of bootstrap draws for combined estimator
                               base_seed = 1L,               # base seed for reproducibility
                               show_progress = TRUE,         # show progress bar
                               n_cores = NULL) {             # number of parallel workers

  # -----------------------------
  # 0) Basic input checks
  # -----------------------------

  # population must contain these variables
  need_pop <- c("Y", "X", "Z", "U")

  # stop if any are missing
  miss <- setdiff(need_pop, names(pop_df))
  if (length(miss) > 0) stop(paste("Population df missing:", paste(miss, collapse = ", ")))

  # b1 must be a single finite number
  if (!is.numeric(b1) || length(b1) != 1 || !is.finite(b1)) {
    stop("b1 must be a single finite numeric value.")
  }

  # b6 must be a single finite number
  if (!is.numeric(b6) || length(b6) != 1 || !is.finite(b6)) {
    stop("b6 must be a single finite numeric value.")
  }

  # base_seed must be a single finite number
  if (!is.numeric(base_seed) || length(base_seed) != 1 || !is.finite(base_seed)) {
    stop("base_seed must be a single finite numeric value.")
  }

  # convert seed to integer (just to be safe)
  base_seed <- as.integer(base_seed)

  # pbapply must be installed
  if (!requireNamespace("pbapply", quietly = TRUE)) {
    stop("Package 'pbapply' is required. Install it with install.packages('pbapply').")
  }

  # -----------------------------
  # 1) Decide how many workers
  # -----------------------------

  # if n_cores is NULL, use all logical cores minus 1
  if (is.null(n_cores)) {
    n_cores <- max(1L, parallel::detectCores(logical = TRUE) - 1L)
  }

  # enforce integer >= 1
  n_cores <- as.integer(n_cores)
  if (n_cores < 1L) stop("n_cores must be >= 1")

  # -----------------------------
  # 2) Create the task grid
  # -----------------------------

  # each task is one pair (aU index i, replication r)
  grid <- expand.grid(
    i = seq_along(aU_vals),
    rep = seq_len(R),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  # total number of tasks
  n_tasks <- nrow(grid)

  # -----------------------------
  # 3) Start a PSOCK cluster
  # -----------------------------

  # PSOCK works on Windows/Mac/Linux
  cl <- parallel::makeCluster(n_cores)

  # always stop cluster at the end (also on error)
  on.exit(parallel::stopCluster(cl), add = TRUE)

  # set a reproducible RNG stream across workers
  parallel::clusterSetRNGStream(cl, iseed = base_seed)

  # -----------------------------
  # 4) Export objects to workers
  # -----------------------------

  # workers need:
  #   - the population data
  #   - constants and tuning parameters
  #   - your sampling + analysis functions
  parallel::clusterExport(
    cl,
    varlist = c(
      "pop_df", "b1", "b6", "n_srs", "n_nps", "aX", "aY", "aZ", "aU_vals", "boot_B", "base_seed",
      "srs", "nps",
      "analyze_srs", "analyze_nps", "analyze_two_step"
    ),
    envir = environment()
  )

  # -----------------------------
  # 5) One-task worker function
  # -----------------------------

  # this runs exactly one (aU, rep) combination and returns 3 rows:
  #   - SRS estimate
  #   - NPS estimate
  #   - Combined/IPW estimate
  worker_fun <- function(task_row) {

    # unpack the task identifiers
    i <- task_row[["i"]]
    r <- task_row[["rep"]]

    # current aU value
    aU <- aU_vals[i]

    # deterministic seeds per task:
    # ensures identical results regardless of number of cores / scheduling
    seed_srs  <- base_seed + 100000L * i + 10L * r + 1L
    seed_nps  <- base_seed + 100000L * i + 10L * r + 2L
    seed_boot <- base_seed + 100000L * i + 10L * r + 3L

    # draw the SRS
    srs_df <- srs(pop_df, n = n_srs, seed = seed_srs)

    # draw the NPS with selection parameters (including current aU)
    nps_df <- nps(pop_df, n = n_nps, aX = aX, aY = aY, aZ = aZ, aU = aU, seed = seed_nps)

    # fit SRS outcome model
    res_srs <- analyze_srs(srs_df)

    # extract b1/b6 pieces if your analyze_* functions return them
    # (recommended: analyze_* fits Y ~ X + I(X^2 - 1) + ... and returns both coefficients + SEs)
    b1hat_srs <- if (!is.null(res_srs$beta_X)) res_srs$beta_X else NA_real_
    se1_srs   <- if (!is.null(res_srs$se_X))   res_srs$se_X   else NA_real_
    b6hat_srs <- if (!is.null(res_srs$beta_X2c)) res_srs$beta_X2c else NA_real_
    se6_srs   <- if (!is.null(res_srs$se_X2c))   res_srs$se_X2c   else NA_real_

    # store SRS result
    row_srs <- data.frame(
      aU = aU,
      rep = r,
      method = "SRS: OLS(Y~X+I(X^2-1)+Z)",
      beta_hat_b1 = b1hat_srs,
      se_hat_b1 = se1_srs,
      bias_b1 = b1hat_srs - b1,
      beta_hat_b6 = b6hat_srs,
      se_hat_b6 = se6_srs,
      bias_b6 = b6hat_srs - b6
    )

    # fit NPS outcome model (includes U)
    res_nps <- analyze_nps(nps_df)

    # extract b1/b6 pieces if your analyze_* functions return them
    b1hat_nps <- if (!is.null(res_nps$beta_X)) res_nps$beta_X else NA_real_
    se1_nps   <- if (!is.null(res_nps$se_X))   res_nps$se_X   else NA_real_
    b6hat_nps <- if (!is.null(res_nps$beta_X2c)) res_nps$beta_X2c else NA_real_
    se6_nps   <- if (!is.null(res_nps$se_X2c))   res_nps$se_X2c   else NA_real_

    # store NPS result
    row_nps <- data.frame(
      aU = aU,
      rep = r,
      method = "NPS: OLS(Y~X+I(X^2-1)+Z+U)",
      beta_hat_b1 = b1hat_nps,
      se_hat_b1 = se1_nps,
      bias_b1 = b1hat_nps - b1,
      beta_hat_b6 = b6hat_nps,
      se_hat_b6 = se6_nps,
      bias_b6 = b6hat_nps - b6
    )

    # fit combined two-step estimator (bootstrap SE)
    res_ipw <- analyze_two_step(srs_df, nps_df, B = boot_B, seed = seed_boot)

    # extract b1/b6 pieces if your combined function returns them
    b1hat_ipw <- if (!is.null(res_ipw$beta_X)) res_ipw$beta_X else NA_real_
    se1_ipw   <- if (!is.null(res_ipw$se_boot_X)) res_ipw$se_boot_X else if (!is.null(res_ipw$se_boot)) res_ipw$se_boot else NA_real_
    b6hat_ipw <- if (!is.null(res_ipw$beta_X2c)) res_ipw$beta_X2c else NA_real_
    se6_ipw   <- if (!is.null(res_ipw$se_boot_X2c)) res_ipw$se_boot_X2c else NA_real_

    # store combined result
    row_ipw <- data.frame(
      aU = aU,
      rep = r,
      method = "Combined: IPW + U (boot SE)",
      beta_hat_b1 = b1hat_ipw,
      se_hat_b1 = se1_ipw,
      bias_b1 = b1hat_ipw - b1,
      beta_hat_b6 = b6hat_ipw,
      se_hat_b6 = se6_ipw,
      bias_b6 = b6hat_ipw - b6
    )

    # return 3-row data.frame for this task
    rbind(row_srs, row_nps, row_ipw)
  }

  # -----------------------------
  # 6) Run tasks in parallel (with pbapply progress bar)
  # -----------------------------

  # pbapply controls whether the progress bar is shown
  pbapply::pboptions(type = if (show_progress) "timer" else "none")

  # pblapply will use the cluster and show a progress bar in the master
  res_list <- pbapply::pblapply(
    X = split(grid, seq_len(n_tasks)),   # list of 1-row data.frames
    FUN = worker_fun,
    cl = cl
  )

  # -----------------------------
  # 7) Combine all replication draws
  # -----------------------------

  draws <- do.call(rbind, res_list)

  # -----------------------------
  # 8) Summarize for plotting
  # -----------------------------

  summary <- aggregate(
    cbind(beta_hat_b1, se_hat_b1, bias_b1, beta_hat_b6, se_hat_b6, bias_b6) ~ aU + method,
    data = draws,
    FUN = mean,
    na.rm = TRUE
  )

  # -----------------------------
  # 9) Return everything
  # -----------------------------

  list(
    summary = summary,
    draws = draws,
    settings = list(
      b1 = b1,
      b6 = b6,
      n_srs = n_srs,
      n_nps = n_nps,
      aX = aX, aY = aY, aZ = aZ,
      aU_vals = aU_vals,
      R = R,
      boot_B = boot_B,
      base_seed = base_seed,
      n_cores = n_cores
    )
  )
}