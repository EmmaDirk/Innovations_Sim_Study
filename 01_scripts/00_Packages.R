# 00_packages.R
# use this script to load any package needed 
# -----------------------------------------------------------------

# required packages
pkgs <- c("tidyverse", "sampling", "viridis",
          "patchwork", "pbapply")

# load packages
lapply(pkgs, library, character.only = TRUE)