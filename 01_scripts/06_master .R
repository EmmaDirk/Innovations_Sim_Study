# this script is used to call all other scripts
# -----------------------------------------------------------------

# load here package for file paths
library("here")

# step 0: load packages
# load packages
source(here("01_scripts", "00_packages.R"))

# step 1: create population data
# load population data generation function
source(here("01_scripts", "01_population_data_generation.R"))

# reproducibility
set.seed(123)

# we want 4 populations. 
# 1a. population where U has an effect of 0.05 on X and Y
df_U_005 <- simulate_pop_data(N = 1000000,                   # population of 1M
                              b1 = 0.1,                      # causal estimand X -> Y
                              b2 = 0.15,                     # effect of Z -> X
                              b3 = 0.15,                     # effect of Z -> Y
                              b4 = 0.05,                     # effect of U -> Y
                              b5 = 0.05)                     # effect of U -> X
  
# 1b. population where U has an effect of 0.1 on X and Y
df_U_010 <- simulate_pop_data(N = 1000000,                   # population of 1M
                              b1 = 0.1,                      # causal estimand X -> Y
                              b2 = 0.15,                     # effect of Z -> X
                              b3 = 0.15,                     # effect of Z -> Y
                              b4 = 0.1,                      # effect of U -> Y
                              b5 = 0.1)                      # effect of U -> X

# 1c. population where U has an effect of 0.2 on X and Y
df_U_020 <- simulate_pop_data(N = 1000000,                   # population of 1M
                              b1 = 0.1,                      # causal estimand X -> Y
                              b2 = 0.15,                     # effect of Z -> X
                              b3 = 0.15,                     # effect of Z -> Y
                              b4 = 0.2,                      # effect of U -> Y
                              b5 = 0.2)                      # effect of U -> X

# 1d. population where U has an effect of 0.4 on X and Y
df_U_040 <- simulate_pop_data(N = 1000000,                   # population of 1M
                              b1 = 0.1,                      # causal estimand X -> Y
                              b2 = 0.15,                     # effect of Z -> X
                              b3 = 0.15,                     # effect of Z -> Y
                              b4 = 0.4,                      # effect of U -> Y
                              b5 = 0.4)                      # effect of U -> X

# save the dataframes in the data folder
save(df_U_005, file = here("02_data", "df_U_005.RData"))
save(df_U_010, file = here("02_data", "df_U_010.RData"))
save(df_U_020, file = here("02_data", "df_U_020.RData"))
save(df_U_040, file = here("02_data", "df_U_040.RData"))

# step 2: run the analyses
# load the required functions
source(here("01_scripts", "02_sampling.R"))
source(here("01_scripts", "03_analysis.R"))
source(here("01_scripts", "04_simulator.R"))

# we will run the analyses on the 4 datasets
# 2a. 
out_U_005 <- simulation_wrapper(
  pop_df = df_U_005,                                               # population data frame
  b1 = 0.1,                                                        # TRUE effect for bias calc
  n_srs = 1000,                                                    # sample size for SRS
  n_nps = 200,                                                     # sample sizes for NPS
  aX = 0.1,                                                        # selection effect of X (small)
  aY = 0.15,                                                       # selection effect of Y (small)
  aZ = 0.5,                                                        # selection effect of Z (large)
  aU_vals = seq(0, 1, by = 0.2),                                   # vector of effect of U on selection
  R = 1000,                                                        # number of replications
  boot_B = 50,                                                     # number of bootstrap draws
  base_seed = 1L,                                                  # seed for reproducibility
  show_progress = TRUE,                                            # show progress bar
  n_cores = 7                                                      # parallel workers (NULL = detectCores()-1)
)

# save the results
save(out_U_005, file = here("03_output", "out_U_005.RData"))

# 2b.
out_U_010 <- simulation_wrapper(
  pop_df = df_U_010,                                               # population data frame
  b1 = 0.1,                                                        # TRUE effect for bias calc
  n_srs = 1000,                                                    # sample size for SRS
  n_nps = 200,                                                     # sample sizes for NPS
  aX = 0.1,                                                        # selection effect of X (small)
  aY = 0.15,                                                       # selection effect of Y (small)
  aZ = 0.5,                                                        # selection effect of Z (large)
  aU_vals = seq(0, 1, by = 0.2),                                   # vector of effect of U on selection
  R = 1000,                                                        # number of replications
  boot_B = 50,                                                     # number of bootstrap draws
  base_seed = 1L,                                                  # seed for reproducibility
  show_progress = TRUE,                                            # show progress bar
  n_cores = 7                                                      # parallel workers (NULL = detectCores()-1)
)

# save the results
save(out_U_010, file = here("03_output", "out_U_010.RData"))

# 2c.
out_U_020 <- simulation_wrapper(
  pop_df = df_U_020,                                               # population data frame
  b1 = 0.1,                                                        # TRUE effect for bias calc
  n_srs = 1000,                                                    # sample size for SRS
  n_nps = 200,                                                     # sample sizes for NPS
  aX = 0.1,                                                        # selection effect of X (small)
  aY = 0.15,                                                       # selection effect of Y (small)
  aZ = 0.5,                                                        # selection effect of Z (large)
  aU_vals = seq(0, 1, by = 0.2),                                   # vector of effect of U on selection
  R = 1000,                                                        # number of replications
  boot_B = 50,                                                     # number of bootstrap draws
  base_seed = 1L,                                                  # seed for reproducibility
  show_progress = TRUE,                                            # show progress bar
  n_cores = 7                                                      # parallel workers (NULL = detectCores()-1)
)

# save the results
save(out_U_020, file = here("03_output", "out_U_020.RData"))

# 2d.
out_U_040 <- simulation_wrapper(
  pop_df = df_U_040,                                               # population data frame
  b1 = 0.1,                                                        # TRUE effect for bias calc
  n_srs = 1000,                                                    # sample size for SRS
  n_nps = 200,                                                     # sample sizes for NPS
  aX = 0.1,                                                        # selection effect of X (small)
  aY = 0.15,                                                       # selection effect of Y (small)
  aZ = 0.5,                                                        # selection effect of Z (large)
  aU_vals = seq(0, 1, by = 0.2),                                   # vector of effect of U on selection
  R = 1000,                                                        # number of replications
  boot_B = 50,                                                     # number of bootstrap draws
  base_seed = 1L,                                                  # seed for reproducibility
  show_progress = TRUE,                                            # show progress bar
  n_cores = 7                                                      # parallel workers (NULL = detectCores()-1)
)

# save the results
save(out_U_040, file = here("03_output", "out_U_040.RData"))

# step 3: create plots
source(here("01_scripts", "05_plotting.R"))

# 3a. U-confounding = 0.05
plot_bias_curves(out_U_005, main = "Bias vs aU (U-confounding = 0.05)")
plot_rmse_curves(out_U_005, main = "RMSE vs aU (U-confounding = 0.05)")

# 3b. U-confounding = 0.1
plot_bias_curves(out_U_010, main = "Bias vs aU (U-confounding = 0.1)")
plot_rmse_curves(out_U_010, main = "RMSE vs aU (U-confounding = 0.1)")

# 3c. U-confounding = 0.2
plot_bias_curves(out_U_020, main = "Bias vs aU (U-confounding = 0.2)")
plot_rmse_curves(out_U_020, main = "RMSE vs aU (U-confounding = 0.2)")

# 3d. U-confounding = 0.4
plot_bias_curves(out_U_040, main = "Bias vs aU (U-confounding = 0.4)")
plot_rmse_curves(out_U_040, main = "RMSE vs aU (U-confounding = 0.4)")