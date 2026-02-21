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


