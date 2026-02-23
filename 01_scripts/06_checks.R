# this script is meant to check the functions
# -----------------------------------------------------------------

# libs
library(here)

# -----------------------------------------------------------------
# check script 0: packages
# should load tidyverse and sampling
source(here("01_scripts", "00_packages.R"))

# should be TRUE
"isNamespaceLoaded"("sampling")
"isNamespaceLoaded"("tidyverse")

# -----------------------------------------------------------------
# check script 1: population data generation
source(here("01_scripts", "01_population_data_generation.R"))

# example use
set.seed(123)
df <- simulate_pop_data(
  N  = 100000,
  b1 = 0.15,
  b2 = 0.15,
  b3 = 0.15,
  b4 = 0.15,
  b5 = 0.15,
  b6 = 0.10
)

# checking data
# var=1, and covariances
cov(df[, c("Z","U","X","Y")])

# check if parameters are recoverable under the correct model
# b1 = 0.15, b3 = 0.15, b4 = 0.15, b6 = 0.10, so we expect to recover these values when regressing
# Y on X, I(X^2 - 1), Z, and U
my <- lm(Y ~ X + I(X^2 - 1) + Z + U, data = df)
coef(my)

# check if parameters are recoverable under the correct model
# b2 = 0.15, b5 = 0.15, so we expect to recover these values when regressing X on Z and U
mX <- lm(X ~ Z + U, data = df)
coef(mX)

# -----------------------------------------------------------------
# check script 2: sampling
source(here("01_scripts", "02_sampling.R"))

# first simulate data
df <- simulate_pop_data(
  N  = 100000,
  b1 = 0.15,
  b2 = 0.15,
  b3 = 0.15,
  b4 = 0.15,
  b5 = 0.15
)

# sample with the functions
df_srs <- srs(df, n = 1000)
df_nps <- nps(
  df,
  n  = 1000,
  aX = 0.15,
  aY = 0.15,
  aZ = 0.15,
  aU = 0.15
)

# check that the linear model on the SRS is biased because U is omitted
lm(Y ~ X + Z, data = df_srs)

# check that the linear model on the NPS is biased because selection
lm(Y ~ X + Z + U, data = df_nps)

# -----------------------------------------------------------------
# check script 3: analysis
source(here("01_scripts", "03_analysis.R"))

# generate data
df_pop <- simulate_pop_data(
  N  = 100000,
  b1 = 0.3,  
  b2 = 0.3,
  b3 = 0,
  b4 = 0.4,
  b5 = 0.4,
  b6 = 0.1
)

# check if parameters are recoverable under the correct model
true_fit <- lm(Y ~ X + Z + U, data = df_pop)
true_b1  <- coef(true_fit)["X"]
cat("True b1 (population):", true_b1, "\n\n")

# draw samples
srs_df <- srs(df_pop, n = 1000, seed = 1)
nps_df <- nps(
  df_pop,
  n  = 1000,
  aX = 0.8,
  aY = 0.9,
  aZ = 0.7,
  aU = 0.1,
  seed = 2
)

# analyze SRS: we expect bias
res_srs <- analyze_srs(srs_df)
cat("SRS results\n")
cat("beta_X:", res_srs$beta_X, "\n")
cat("se_X  :", res_srs$se_X, "\n\n")

# analyze NPS: we expect bias
res_nps <- analyze_nps(nps_df)
cat("NPS results\n")
cat("beta_X:", res_nps$beta_X, "\n")
cat("se_X  :", res_nps$se_X, "\n\n")

# analyze two-setp: we expect less bias
res_two <- analyze_two_step(
  srs_df,
  nps_df,
  B = 300,   
  ps_formula = ~ X + Z + Y,
  seed = 3
)

cat("Two-step IPW results\n")
cat("beta_X:", res_two$beta_X, "\n")
cat("se_X  :", res_two$se_X, "\n\n")