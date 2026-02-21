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
#   n_cores        : number of parallel workers (NULL uses detectCores()-1)
#
# Returns:
#   list(summary=..., draws=..., settings=...)
#     - draws: replication-level rows with beta_hat, se_hat, bias
#     - summary: averages by aU and method (mean beta_hat, mean se_hat, mean bias)
# -----------------------------------------------------------------

simulation_wrapper <- function(pop_df,                       # population data frame
                               b1,                           # TRUE effect for bias calc
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
      "pop_df", "b1", "n_srs", "n_nps", "aX", "aY", "aZ", "aU_vals", "boot_B", "base_seed",
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

    # store SRS result
    row_srs <- data.frame(
      aU = aU,
      rep = r,
      method = "SRS: OLS(Y~X+Z)",
      beta_hat = res_srs$beta_X,
      se_hat = res_srs$se_X,
      bias = res_srs$beta_X - b1
    )

    # fit NPS outcome model (includes U)
    res_nps <- analyze_nps(nps_df)

    # store NPS result
    row_nps <- data.frame(
      aU = aU,
      rep = r,
      method = "NPS: OLS(Y~X+Z+U)",
      beta_hat = res_nps$beta_X,
      se_hat = res_nps$se_X,
      bias = res_nps$beta_X - b1
    )

    # fit combined two-step estimator (bootstrap SE)
    res_ipw <- analyze_two_step(srs_df, nps_df, B = boot_B, seed = seed_boot)

    # store combined result
    row_ipw <- data.frame(
      aU = aU,
      rep = r,
      method = "Combined: IPW + U (boot SE)",
      beta_hat = res_ipw$beta_X,
      se_hat = res_ipw$se_X,
      bias = res_ipw$beta_X - b1
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
    cbind(beta_hat, se_hat, bias) ~ aU + method,
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
