# ============================== 00_config.R ==============================

# --- package setup ---
# check if packages are installed; if not, install
pkgs <- c("readr", "dplyr")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
}

# load packages
library(readr)
library(dplyr)
library(tidyr)

# --- load helper code ---
#source("code/utils.R")

# --- import experiment data ---
# read the csv of experiment jumbles and take a small random sample
# note that with no seed set, samples will not be identical
exp_data <- read_csv("data/exp_df.csv") |> 
  #filter(week %in% c(3, 4)) |> 
  filter(jumbler == "mp")
  

# --- import real data ---
real_data <- read_csv("data/combos_release_2_20241024.csv")

# --- import pin data --- # not actually helpful - abandoned this signal
pin_data <- read_csv("data/4_pin_leaks.csv")
