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

  # sampling package is required
  if (!requireNamespace("sampling", quietly = TRUE)) {
    stop("Package 'sampling' is required. Install it with install.packages('sampling').")
  }

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

  # compute "target" inclusion probabilities p_i
  p <- plogis(a0 + aX*df$X + aY*df$Y + aZ*df$Z + aU*df$U)


  # We want a fixed sample size exactly n.
  # Many fixed-size algorithms in {sampling} assume sum(pik) == n.
  # So we (i) scale to sum n and (ii) lightly clip away from 0/1.

  # scale to sum exactly n (up to floating error)
  pik <- p * (n / sum(p))

  # clip to (eps, 1-eps) to keep algorithms stable
  eps <- 1e-6
  pik <- pmin(pmax(pik, eps), 1 - eps)

  # re-scale again after clipping
  pik <- pik * (n / sum(pik))

  # if clipping caused tiny drift outside [0,1], clip once more
  pik <- pmin(pmax(pik, eps), 1 - eps)

  # ensure sum(pik) is as close as possible to n (and integer-sized designs work)
  # final tiny correction on the largest-probability unit
  drift <- n - sum(pik)
  if (abs(drift) > 1e-8) {
    j <- which.max(pik)
    pik[j] <- pik[j] + drift
    pik[j] <- min(max(pik[j], eps), 1 - eps)
  }

  # draw a fixed-size sample:
  # UPsystematic returns a 0/1 vector of length N with exactly sum(pik) ones (when sum is integer).
  S <- sampling::UPsystematic(pik)
  idx <- which(S == 1)

  # (very rare) safety: if due to numeric quirks we didn't get exactly n, repair deterministically
  if (length(idx) != n) {
    if (length(idx) > n) {
      idx <- sample(idx, size = n, replace = FALSE)
    } else {
      remaining <- setdiff(seq_len(N), idx)
      need <- n - length(idx)
      w <- pik[remaining]
      w <- w / sum(w)
      idx <- c(idx, sample(remaining, size = need, replace = FALSE, prob = w))
    }
    S <- integer(N); S[idx] <- 1L
  }

  # Return the sample
  out <- df[idx, , drop = FALSE]

  # Store selection info for debugging/analysis
  attr(out, "selection") <- list(
    a0 = a0,
    a  = c(aX = aX, aY = aY, aZ = aZ, aU = aU),
    p_all   = p,     # original logistic probs (sum ~= n)
    pik_all = pik,   # fixed-size design inclusion probs (sum == n)
    S_all   = S      # realized sample membership
  )

  out
}

