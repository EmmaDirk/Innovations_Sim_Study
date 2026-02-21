# this script is used to create the samples:
# - Simple random sampling (SRS)
# - Non-probability sampling (NPS) with selection bias
# -----------------------------------------------------------------

# function to create SRS:
srs <- function(df,                      # data frame (population)
                n,                       # sample size
                seed = NULL) {
  
  # if no seed is provided, set the seed
  if (!is.null(seed)) set.seed(seed)

  # if n is larger than nrow(df), stop
  N <- nrow(df)
  if (n > N) stop("n cannot be larger than nrow(df).")
  
  # sample with replacement
  idx <- sample.int(N, size = n, replace = FALSE)
  out <- df[idx, , drop = FALSE]
  
  # drop U (if present)
  if ("U" %in% names(out)) out$U <- NULL
  
  out
}

# a function to create NPS. 
nps <- function(df,                       # data frame (population)
                n,                        # sample size
                aX = 0,                   # effect of X on selection
                aY = 0,                   # effect of Y on selection
                aZ = 0,                   # effect of Z on selection
                aU = 0,                   # effect of U on selection
                seed = NULL) {
  
  # if no seed is provided, set the seed
  if (!is.null(seed)) set.seed(seed)
  
  # if n is larger than nrow(df), stop
  N <- nrow(df)
  if (n > N) stop("n cannot be larger than nrow(df).")
  
  # linear selection index
  eta <- aX*df$X + aY*df$Y + aZ*df$Z + aU*df$U
  
  # exponential weights (monotone in eta)
  w <- exp(eta)
  
  # draw exactly n without replacement
  idx <- sample.int(N, size = n, replace = FALSE, prob = w)
  
  # sample with replacement
  out <- df[idx, , drop = FALSE]
  
  # store weights for diagnostics
  attr(out, "selection") <- list(
    a = c(aX=aX, aY=aY, aZ=aZ, aU=aU),
    weights_all = w
  )
  
  return(out)
}

# checks 
# first simulate data
# df <- simulate_pop_data(N = 100000, 
#                  b1 = 0.15, 
#                  b2 = 0.15, 
#                  b3 = 0.15,
#                  b4 = 0.15, 
#                  b5 = 0.15)

# sample with the functions
# df_srs <- srs(df, n = 1000)
# df_nps <- nps(df, 
#              n = 1000,
#              aX = 0.15, 
#              aY = 0.15, 
#              aZ = 0.15, 
#              aU = 0.15)

# check that the linear model on the SRS is biased because U is omitted
# lm(Y ~ X + Z, data = df_srs)

# check that the linear model on the NPS is biased because selection
# lm(Y ~ X + Z + U, data = df_nps)
