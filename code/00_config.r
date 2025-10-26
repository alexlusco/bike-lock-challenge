# ------------------------------------------------------------
# author: alex luscombe & m.p.
# date: october 2025
#
# description:
# loads required packages, pulls in helper code, 
# reads a small experimental data set used to build the 
# intial model
# ------------------------------------------------------------

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
source("code/utils.R")

# --- import experiment data ---
# read the csv of experiment jumbles and take a small random sample
# note that with no seed set, samples will not be identical
exp_data <- read_csv("data/exp_df.csv") |> 
  slice_sample(n = 20)

# --- import real data ---
real_data <- read_csv("data/combos_release_1_20241024.csv") |> mutate(week = 1)

# --- import pin data --- # not actually helpful - abandoned this signal
#pin_data <- read_csv("data/4_pin_leaks.csv")
