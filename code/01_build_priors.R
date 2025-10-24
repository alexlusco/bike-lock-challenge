# ------------------------------------------------------------
# author: alex luscombe & m.p.
# date: october 2025
#
# this script builds probabilistic priors for each wheel of a 
# 4-digit bike lock based on observed "jumble" data. it uses 
# entropy to estimate how much each wheel tends to spin and 
# applies a kernel to weight digits near frequently stopped 
# positions. from those priors, it scores all 10,000 possible 
# combinations, optionally banning certain starting digits.
#
# after each real-world guess, the script takes feedback on 
# how many digits were in the correct position (k) and filters 
# the candidate pool accordingly. multiple guesses can be 
# chained to iteratively shrink the search space. the final 
# output is a ranked list of remaining likely combinations.
# ------------------------------------------------------------

# build priors for each wheel
# if spin_amount = NULL, it will estimate spin automatically from entropy
# lambda_sharp and lambda_flat define the sharpness range for that mapping
# row_weights can downweight certain rows (e.g., rushed or gloved observations)
priors <- build_priors(
  jumbles      = exp_data,
  spin_amount  = NULL,        # override with c(a,b,c,d) to set spins manually
  lambda_sharp = 1.1,
  lambda_flat  = 0.5,
  row_weights  = NULL
)

# score all 10,000 possible codes using the priors
# you can optionally ban starting digits (here 7, 8, 9)
w <- code_weights_from_priors(
  priors,
  combos_universe,
  first_digit_prohibit = c(7,8,9)
)

# list of guesses made so far and feedback k-values (how many digits correct)
guesses <- list(c(2,3,6,7))
ks      <- c(1)

# rank and filter the universe based on the guesses and feedback
# this function keeps only codes consistent with all guesses/k-values
# weights are subset and renormalized to match remaining rows
ranked <- rank_universe(
  combos_universe,
  w,
  guesses = guesses,
  ks = ks
)

# number of possible codes left after filtering
nrow(ranked)

# show the top 10 highest weighted surviving combinations
head(ranked, 10)