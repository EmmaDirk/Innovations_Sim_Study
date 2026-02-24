# 05_master .R
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
df_U_005_int <- simulate_pop_data(N = 1000000,                   # population of 1M
                                  b1 = 0.1,                      # causal estimand X -> Y
                                  b2 = 0.15,                     # effect of Z -> X
                                  b3 = 0.15,                     # effect of Z -> Y
                                  b4 = 0.05,                     # effect of U -> Y
                                  b5 = 0.05,                     # effect of U -> X
                                  b6 = 0.1)                        # curvature term for X^2 -> Y
  
# 1b. population where U has an effect of 0.1 on X and Y
df_U_010_int <- simulate_pop_data(N = 1000000,                   # population of 1M
                                  b1 = 0.1,                      # causal estimand X -> Y
                                  b2 = 0.15,                     # effect of Z -> X
                                  b3 = 0.15,                     # effect of Z -> Y
                                  b4 = 0.1,                      # effect of U -> Y
                                  b5 = 0.1,                      # effect of U -> X
                                  b6 = 0.1)                      # curvature term for X^2 -> Y

# 1c. population where U has an effect of 0.2 on X and Y
df_U_020_int <- simulate_pop_data(N = 1000000,                   # population of 1M
                                  b1 = 0.1,                      # causal estimand X -> Y
                                  b2 = 0.15,                     # effect of Z -> X
                                  b3 = 0.15,                     # effect of Z -> Y
                                  b4 = 0.2,                      # effect of U -> Y
                                  b5 = 0.2,                      # effect of U -> X
                                  b6 = 0.1)                      # curvature term for X^2 -> Y

# 1d. population where U has an effect of 0.4 on X and Y
df_U_040_int <- simulate_pop_data(N = 1000000,                   # population of 1M
                                  b1 = 0.1,                      # causal estimand X -> Y
                                  b2 = 0.15,                     # effect of Z -> X
                                  b3 = 0.15,                     # effect of Z -> Y
                                  b4 = 0.4,                      # effect of U -> Y
                                  b5 = 0.4,                      # effect of U -> X
                                  b6 = 0.1)                      # curvature term for X^2 -> Y

# save the dataframes in the data folder
save(df_U_005_int, file = here("02_data", "df_U_005_int.RData"))
save(df_U_010_int, file = here("02_data", "df_U_010_int.RData"))
save(df_U_020_int, file = here("02_data", "df_U_020_int.RData"))
save(df_U_040_int, file = here("02_data", "df_U_040_int.RData"))

# step 2: run the analyses
# load the required functions
source(here("01_scripts", "02_sampling.R"))
source(here("01_scripts", "03_analysis.R"))
source(here("01_scripts", "04_simulator.R"))

# we will run the analyses on the 4 datasets
# 2a. 
out_U_005_int <- simulation_wrapper(
  pop_df = df_U_005_int,                                           # population data frame
  b1 = 0.1,                                                        # TRUE effect for bias calc
  n_srs = 1000,                                                    # sample size for SRS
  n_nps = c(200, 5000),                                            # sample sizes for NPS (small NPS + big data)
  aX = 0.3,                                                        # selection effect of X (mild selection   - OR=1.35)
  aY = 0.5,                                                        # selection effect of Y (medium selection - OR=1.65)
  aZ = 0.8,                                                        # selection effect of Z (strong selection - OR=2.23)
  aU_vals = seq(0, 2, by = 0.1),                                   # vector of effect of U on selection
  R = 200,                                                         # number of replications
  boot_B = 50,                                                     # number of bootstrap draws
  base_seed = 1L,                                                  # seed for reproducibility
  show_progress = TRUE,                                            # show progress bar
  n_cores = 6                                                      # parallel workers (NULL = detectCores()-1)
)

# save the results
save(out_U_005_int, file = here("03_output", "out_U_005_int.RData"))

# 2b.
out_U_010_int <- simulation_wrapper(
  pop_df = df_U_010_int,                                           # population data frame
  b1 = 0.1,                                                        # TRUE effect for bias calc
  n_srs = 1000,                                                    # sample size for SRS
  n_nps = c(200, 5000),                                            # sample sizes for NPS (small NPS + big data)
  aX = 0.1,                                                        # selection effect of X (small)
  aY = 0.15,                                                       # selection effect of Y (small)
  aZ = 0.5,                                                        # selection effect of Z (large)
  aU_vals = seq(0, 2, by = 0.1),                                   # vector of effect of U on selection
  R = 200,                                                         # number of replications
  boot_B = 50,                                                     # number of bootstrap draws
  base_seed = 1L,                                                  # seed for reproducibility
  show_progress = TRUE,                                            # show progress bar
  n_cores = 7                                                      # parallel workers (NULL = detectCores()-1)
)

# save the results
save(out_U_010_int, file = here("03_output", "out_U_010_int.RData"))

# 2c.
out_U_020_int <- simulation_wrapper(
  pop_df = df_U_020_int,                                           # population data frame
  b1 = 0.1,                                                        # TRUE effect for bias calc
  n_srs = 1000,                                                    # sample size for SRS
  n_nps = c(200, 5000),                                            # sample sizes for NPS (small NPS + big data)
  aX = 0.1,                                                        # selection effect of X (small)
  aY = 0.15,                                                       # selection effect of Y (small)
  aZ = 0.5,                                                        # selection effect of Z (large)
  aU_vals = seq(0, 2, by = 0.1),                                   # vector of effect of U on selection
  R = 200,                                                         # number of replications
  boot_B = 50,                                                     # number of bootstrap draws
  base_seed = 1L,                                                  # seed for reproducibility
  show_progress = TRUE,                                            # show progress bar
  n_cores = 7                                                      # parallel workers (NULL = detectCores()-1)
)

# save the results
save(out_U_020_int, file = here("03_output", "out_U_020_int.RData"))

# 2d.
out_U_040_int <- simulation_wrapper(
  pop_df = df_U_040_int,                                           # population data frame
  b1 = 0.1,                                                        # TRUE effect for bias calc
  n_srs = 1000,                                                    # sample size for SRS
  n_nps = c(200, 5000),                                            # sample sizes for NPS (small NPS + big data)
  aX = 0.1,                                                        # selection effect of X (small)
  aY = 0.15,                                                       # selection effect of Y (small)
  aZ = 0.5,                                                        # selection effect of Z (large)
  aU_vals = seq(0, 2, by = 0.1),                                   # vector of effect of U on selection
  R = 200,                                                         # number of replications
  boot_B = 50,                                                     # number of bootstrap draws
  base_seed = 1L,                                                  # seed for reproducibility
  show_progress = TRUE,                                            # show progress bar
  n_cores = 7                                                      # parallel workers (NULL = detectCores()-1)
)

# save the results
save(out_U_040_int, file = here("03_output", "out_U_040_int.RData"))