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

# a function to create NPS. All variables predict the selection probability.
nps <- function(df,                     # data frame (population)
                n,                      # sample size
                aX,                     # effect of X on selection (log-odds scale)
                aY,                 
                aZ, 
                aU,
                seed = NULL) {

  # if no seed is provided, set the seed
  if (!is.null(seed)) set.seed(seed)

  # if n is larger than nrow(df), stop
  N <- nrow(df)
  if (n > N) stop("n cannot be larger than nrow(df).")

  # step 0: helper: 
  # compute expected number of selected units for a given intercept a0
  expected_n <- function(a0) {

    # compute inclusion probabilities
    p <- plogis(a0 + aX*df$X + aY*df$Y + aZ*df$Z + aU*df$U)

    # expected number is the sum of the inclusion probabilities
    sum(p)
  }

  # find intercept a0 so that sum(p_i) ~= n 
  # expected_n(a0) is monotone increasing in a0, so uniroot works.
  a0 <- uniroot(function(t) expected_n(t) - n, lower = -50, upper = 50)$root

  # compute inclusion probabilities p_i
  p <- plogis(a0 + aX*df$X + aY*df$Y + aZ*df$Z + aU*df$U)

  # realize inclusion indicators S_i (Poisson/Bernoulli sampling)
  S <- rbinom(N, size = 1, prob = p)
  idx <- which(S == 1)

  # we want exactly n units selected 
  if (length(idx) >= n) {

    # too many selected: randomly keep n of them
    idx <- sample(idx, size = n, replace = FALSE)

    # otherwise, too few selected
  } else {

    # top up from non-selected, more likely from high p_i units
    remaining <- setdiff(seq_len(N), idx)
    need <- n - length(idx)

    # compute sampling weights for remaining units (proportional to p_i)
    w <- p[remaining]

    # normalize
    w <- w / sum(w)

    # sample the remaining units with probabilities proportional to p_i
    idx <- c(idx, sample(remaining, size = need, replace = FALSE, prob = w))
  }

  # Return the sample
  out <- df[idx, , drop = FALSE]

  # Store selection info for debugging/analysis
  attr(out, "selection") <- list(
    a0 = a0,
    a = c(aX = aX, aY = aY, aZ = aZ, aU = aU),
    p_all = p,   
    S_all = S    
  )

  out
}
