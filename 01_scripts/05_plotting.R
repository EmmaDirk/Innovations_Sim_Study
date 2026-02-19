# -----------------------------------------------------------------
# Plotting functions for simulation_wrapper() output
#
# Expected input:
#   out <- simulation_wrapper(...)
#   out$summary has columns: aU, method, beta_hat, se_hat, bias
#   out$draws   has columns: aU, rep, method, beta_hat, se_hat, bias
#
# Plots:
#   1) Bias vs aU (lines per method)
#   2) RMSE vs aU (lines per method)
#      - RMSE is computed from replication draws:
#          RMSE(aU, method) = sqrt( mean( (beta_hat - b1)^2 ) )
#        which equals sqrt( mean( bias^2 ) )
#
# Notes:
#   - Uses base R only (no ggplot).
#   - One plot per call (you’ll call it 4 times for your 4 populations).
# -----------------------------------------------------------------

# helper: compute RMSE curves from replication-level draws
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

    # RMSE = sqrt(mean(bias^2))
    bias ~ aU + method,

    # data
    data = d,

    # function
    FUN = function(x) sqrt(mean(x^2, na.rm = TRUE))
  )

  # rename
  names(rmse_df)[names(rmse_df) == "bias"] <- "rmse"

  # return
  rmse_df
}

# plot 1: bias vs aU
plot_bias_curves <- function(sim_out,
                             main = "Bias vs aU",
                             xlab = "aU (effect of U on selection)",
                             ylab = "Bias (beta_hat - b1)",
                             lwd = 2) {

  # checks
  if (is.null(sim_out$summary)) stop("sim_out must contain $summary.")

  # get summary
  s <- sim_out$summary

  # check columns
  need <- c("aU", "method", "bias")

  # if any missing, stop
  miss <- setdiff(need, names(s))
  if (length(miss) > 0) stop(paste("summary missing:", paste(miss, collapse = ", ")))

  # order x values
  aU_vals <- sort(unique(s$aU))
  methods <- unique(s$method)

  # choose colours (one per method)
  cols <- grDevices::hcl.colors(length(methods), palette = "Dark 3")
  names(cols) <- methods

  # y-range over all methods
  ylim <- range(s$bias, finite = TRUE)

  # set up empty plot
  plot(aU_vals, rep(NA_real_, length(aU_vals)),
       type = "n", main = main, xlab = xlab, ylab = ylab, ylim = ylim)

  # add horizontal zero line
  abline(h = 0, lty = 2)

  # draw one line per method
  for (m in methods) {

    # subset
    sm <- s[s$method == m, ]

    # order
    sm <- sm[order(sm$aU), ]

    # draw
    lines(sm$aU, sm$bias, lwd = lwd, col = cols[m])
  }

  # legend with matching colours
  legend("topright", legend = methods, lwd = lwd, col = cols, bty = "n")
}


# plot 2: RMSE vs aU (computed from draws)
plot_rmse_curves <- function(sim_out,
                             main = "RMSE vs aU",
                             xlab = "aU (effect of U on selection)",
                             ylab = "RMSE",
                             lwd = 2) {

  rmse_df <- compute_rmse_summary(sim_out)

  # order x values
  aU_vals <- sort(unique(rmse_df$aU))
  methods <- unique(rmse_df$method)

  # choose colours (one per method)
  cols <- grDevices::hcl.colors(length(methods), palette = "Dark 3")
  names(cols) <- methods

  # y-range
  ylim <- range(rmse_df$rmse, finite = TRUE)

  # set up empty plot
  plot(aU_vals, rep(NA_real_, length(aU_vals)),
       type = "n", main = main, xlab = xlab, ylab = ylab, ylim = ylim)

  # draw one line per method
  for (m in methods) {

    # subset
    dm <- rmse_df[rmse_df$method == m, ]

    # order
    dm <- dm[order(dm$aU), ]

    # draw
    lines(dm$aU, dm$rmse, lwd = lwd, col = cols[m])
  }

  # legend with matching colours
  legend("topright", legend = methods, lwd = lwd, col = cols, bty = "n")
}
