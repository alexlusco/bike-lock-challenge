# ============================== 00_config.R ==============================

# --- package setup ---
# check if packages are installed; if not, install
pkgs <- c("readr", "dplyr", "stringr", "class")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
}

# load packages
library(readr)
library(dplyr)
library(stringr)
library(class)

# --- load helper code ---
source("code/utils.R")

# --- import experiment data ---
# read the csv of experiment jumbles and take a small random sample
# note that with no seed set, samples will not be identical
exp_data <- read_csv("data/exp_df.csv")

source("code/ground_truth.R") # excluded from git repo for security reasons :) if you wish to adapt this code you will need your own jumble sets for "training"

if(gt){
exp_data <- exp_data |> mutate(ground_truth = gt) # this is only available with access to the above sourced ground_truth.R script
}

# --- import real data ---
real_data <- read_csv("data/combos_release_4_20251114.csv") |> rename("index" = `...1`)

# --- import pin data --- # not actually helpful - abandoned this signal
pin_data <- read_csv("data/4_pin_leaks.csv")

# import data for NN template prediction (optional)
templates <- read_csv("data/jumble_templates_with_stats_updated.csv")
game <- read_csv("data/combos_release_3_with_stats.csv")
