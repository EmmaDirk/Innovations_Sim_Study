# -----------------------------------------------------------------
# Plotting functions for simulation_wrapper() output (ggplot2 version)
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
#   - Uses ggplot2 (returns ggplot objects so ggsave() works).
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
                             lwd = 1) {

  # checks
  if (is.null(sim_out$summary)) stop("sim_out must contain $summary.")

  # get summary
  s <- sim_out$summary

  # check columns
  need <- c("aU", "method", "bias")

  # if any missing, stop
  miss <- setdiff(need, names(s))
  if (length(miss) > 0) stop(paste("summary missing:", paste(miss, collapse = ", ")))

  # ensure method is treated as a discrete series
  s$method <- factor(s$method)

  # build plot
  ggplot2::ggplot(s, ggplot2::aes(x = aU, y = bias, color = method)) +

    # add horizontal zero line
    ggplot2::geom_hline(yintercept = 0, linetype = 2) +

    # draw one line per method
    ggplot2::geom_line(linewidth = lwd) +

    # (optional) show points at each aU
    ggplot2::geom_point(size = 2) +

    # labels
    ggplot2::labs(title = main, x = xlab, y = ylab, color = "Method") +

    # clean theme
    ggplot2::theme_minimal()
}


# plot 2: RMSE vs aU (computed from draws)
plot_rmse_curves <- function(sim_out,
                             main = "RMSE vs aU",
                             xlab = "aU (effect of U on selection)",
                             ylab = "RMSE",
                             lwd = 1) {

  # compute RMSE table
  rmse_df <- compute_rmse_summary(sim_out)

  # ensure method is treated as a discrete series
  rmse_df$method <- factor(rmse_df$method)

  # build plot
  ggplot2::ggplot(rmse_df, ggplot2::aes(x = aU, y = rmse, color = method)) +

    # draw one line per method
    ggplot2::geom_line(linewidth = lwd) +

    # (optional) show points at each aU
    ggplot2::geom_point(size = 2) +

    # labels
    ggplot2::labs(title = main, x = xlab, y = ylab, color = "Method") +

    # clean theme
    ggplot2::theme_minimal()
}