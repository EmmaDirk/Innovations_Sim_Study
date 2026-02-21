# this script is used to generate population data
# -----------------------------------------------------------------

simulate_pop_data <- function(N,              # sample size
                              b1,             # X -> Y  (linear part)
                              b2,             # Z -> X
                              b3,             # Z -> Y
                              b4,             # U -> Y
                              b5,             # U -> X
                              b6 = 0,         # NEW: curvature term for X^2 -> Y (heterogeneous effect over X)
                              seed = NULL) {

  # mean-0 Normal variables with Var(Z)=Var(U)=Var(X)=Var(Y)=1 (population).
  #
  # DGM:
  # Z ~ N(0,1), U ~ N(0,1), independent
  # X = b2*Z + b5*U + e_x,  Var(e_x)=1-(b2^2+b5^2)
  #
  # UPDATED Y to allow effect of X to vary over the range of X (nonlinear marginal effect):
  # Y = b1*X + b6*(X^2 - 1) + b3*Z + b4*U + e_y
  #
  # Notes:
  # - (X^2 - 1) is a centered quadratic term (since E[X^2]=1 when Var(X)=1), so it has mean 0 in population.
  # - The marginal effect is dY/dX = b1 + 2*b6*X, so it changes with X.
  #
  # Var(e_y)= 1 - ( b1^2 + b3^2 + b4^2
  #                + 2*b1*b3*b2
  #                + 2*b1*b4*b5
  #                + 2*b6^2 )
  #
  # The extra "+ 2*b6^2" comes from Var(X^2 - 1) = 2 (for mean-0 Var-1 normal X),
  # and cross-covariances with X, Z, U are 0 under this normal DGM.

  if (!is.null(seed)) set.seed(seed)

  # step 1: simulate exogenous variables Z and U
  Z <- rnorm(N, mean = 0, sd = 1)
  U <- rnorm(N, mean = 0, sd = 1)

  # step 2: simulate endogenous variable X
  # first we need to compute how much variance is left for e_x given b2 and b5
  # var(e_x) = sigma2_x
  sx2 <- 1 - (b2^2 + b5^2)

  # sanity check: if sx2 < 0, then Var(e_x) would be negative, which is impossible
  if (sx2 < 0) {
    stop(sprintf("Invalid b's: b2^2 + b5^2 = %.4f > 1, so Var(e_x) would be negative.",
                 b2^2 + b5^2))
  }

  # simulate random errors with mean 0 and variance sx2
  ex <- rnorm(N, mean = 0, sd = sqrt(sx2))

  # simulate X from the DGM
  X <- b2 * Z + b5 * U + ex

  # step 3: simulate endogenous variable Y
  # first we need to compute how much variance is left for e_y given b1, b3, b4, b2, and b5
  #
  # NEW: add centered quadratic term (X^2 - 1) so the effect of X varies over X
  # (marginal effect dY/dX = b1 + 2*b6*X)
  X2c <- X^2 - 1

  sys_var_y <- b1^2 + b3^2 + b4^2 + 2 * b1 * b3 * b2 + 2 * b1 * b4 * b5 +
    2 * b6^2

  # var(e_y) = sigma2_y
  sy2 <- 1 - sys_var_y

  # sanity check: if sy2 < 0, then Var(e_y) would be negative, which is impossible
  if (sy2 < 0) {
    stop(sprintf("Invalid b's: systematic Var(Y) = %.4f > 1, so Var(e_y) would be negative.",
                 sys_var_y))
  }

  # simulate random errors with mean 0 and variance sy2
  ey <- rnorm(N, mean = 0, sd = sqrt(sy2))

  # simulate Y from the DGM
  Y <- b1 * X + b6 * X2c + b3 * Z + b4 * U + ey

  # step 4: combine into a data frame
  out <- data.frame(Z = Z, U = U, X = X, Y = Y)

  # step 5: save parameters as attributes for later reference
  attr(out, "params") <- list(
    b = c(b1 = b1, b2 = b2, b3 = b3, b4 = b4, b5 = b5, b6 = b6),
    sx2 = sx2,
    sy2 = sy2
  )

  # step 6: return the simulated data
  return(out)
}

# example use
# set.seed(123)
# df <- simulate_pop_data(N = 100000,
# b1 = 0.15,
# b2 = 0.15,
# b3 = 0.15,
# b4 = 0.15,
# b5 = 0.15,
# b6 = 0.10)

# checking data
# var=1, and covariances
# cov(df[, c("Z","U","X","Y")])

# check if parameters are recoverable under the correct model
# b1 = 0.15, b3 = 0.15, b4 = 0.15, b6 = 0.10, so we expect to recover these values when regressing
# Y on X, I(X^2 - 1), Z, and U
# my <- lm(Y ~ X + I(X^2 - 1) + Z + U, data = df)
# coef(my)

# check if parameters are recoverable under the correct model
# b2 = 0.15, b5 = 0.15, so we expect to recover these values when regressing X on Z and U
# mX <- lm(X ~ Z + U, data = df)
# coef(mX)


