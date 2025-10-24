# ------------------------------------------------------------
# author: alex luscombe & m.p.
# date: october 2025
#
# what this file does:
# - takes the released “jumble” digits (per wheel) and builds a soft prior for
#   each wheel. idea is: if you often stop near a digit, true digit is more likely
#   to be near there too. we estimate how “spun” each wheel is using entropy,
#   then turn that into a smoothing sharpness (lambda).
# - scores every possible 4-digit code (total of 10,000) by multiplying the 4 wheel priors.
# - optionally bans certain starting digits (e.g., 7/8/9 if folks avoid them).
# - supports a sequence of guesses with k-feedback (“k digits exactly right”)
#   to shrink the candidate set, w/ weight re-normalization.
# ------------------------------------------------------------

# ============================================================
# helpers
# ============================================================

digits <- 0:9

# circular distance on a ring of 10 digits (0..9). max distance is 5.
circ_dist <- function(a, b) {
  d <- abs(a - b)
  pmin(d, 10 - d)
}

# shannon entropy for a length-10 histogram (lower = peaky, higher = flat)
entropy <- function(tab10) {
  p <- tab10 / sum(tab10)
  -sum(ifelse(p > 0, p * log(p), 0))
}

# generate the full 10,000-code universe as a 4-col matrix of digits
# (digits not strings here; we’ll format later if needed)
combos_universe <- as.matrix(
  expand.grid(digit_1 = digits, digit_2 = digits, digit_3 = digits, digit_4 = digits)
)

# map spin in [0,1] → lambda (how sharp the kernel is)
# spin = 0 → lambda_sharp (tight), spin = 1 → lambda_flat (loose)
lambda_from_spin <- function(spin, lambda_sharp = 1.0, lambda_flat = 0.5) {
  stopifnot(length(spin) == 4, all(spin >= 0 & spin <= 1))
  lambda_sharp * (1 - spin) + lambda_flat * spin
}

# estimate spin per wheel from entropy of observed stops.
# higher entropy → more spin (flatter kernel)
auto_spin_from_entropy <- function(jumbles, ent_min = 0.6, ent_max = 2.3) {
  stopifnot(all(c("digit_1","digit_2","digit_3","digit_4") %in% names(jumbles)))
  
  # simple rescaler to [0,1] with clipping
  rescale01 <- function(x, xmin, xmax) {
    z <- (x - xmin) / (xmax - xmin)
    pmin(pmax(z, 0), 1)
  }
  
  # compute entropy for each wheel’s 0..9 histogram
  ents <- vapply(1:4, function(i) {
    tab <- tabulate(jumbles[[paste0("digit_", i)]] + 1, nbins = 10)
    entropy(tab)
  }, numeric(1))
  
  # turn entropies into spins in [0,1]
  rescale01(ents, ent_min, ent_max)
}


# ============================================================
# priors
# ============================================================

# single-wheel prior using a symmetric circular kernel.
# obs_digits: vector of observed end digits (0..9) for one wheel
# lambda: how sharp to weight “nearby” digits (bigger lambda = sharper)
# row_weights: optional per-row weights same length as obs_digits. 
position_prior_kernel <- function(obs_digits, lambda, row_weights = NULL) {
  # histogram of observed stops for digits 0..9 (unweighted or weighted)
  hist <- tabulate(obs_digits + 1, nbins = 10)
  if (!is.null(row_weights)) {
    # re-build histogram as a weighted sum by digit (keeps 10 bins even if some zero)
    hist <- tapply(row_weights, factor(obs_digits, levels = 0:9), sum, default = 0)
    hist <- as.numeric(hist)
  }
  
  # kernel weights for circular distances 0..5. exponential decay.
  K <- exp(-lambda * (0:5))
  
  # score each candidate true digit t by smoothing the histogram around the ring
  scores <- sapply(digits, function(t) {
    s <- 0
    for (x in digits) {
      d <- circ_dist(t, x)
      s <- s + hist[x + 1] * K[d + 1]
    }
    s
  })
  
  # turn scores into a proper probability vector (sums to 1)
  scores / sum(scores)
}

# build priors for all four wheels.
# if spin_amount=NULL, estimate it from entropy; otherwise use the given 4-vector.
# lambda_sharp / lambda_flat define the map from spin → lambda.
build_priors <- function(jumbles,
                         spin_amount = NULL,
                         lambda_sharp = 1.1,
                         lambda_flat  = 0.5,
                         row_weights  = NULL) {
  stopifnot(all(c("digit_1","digit_2","digit_3","digit_4") %in% names(jumbles)))
  
  # decide spin per wheel
  if (is.null(spin_amount)) {
    spin_amount <- auto_spin_from_entropy(jumbles)  # 4 numbers in [0,1]
  } else {
    stopifnot(length(spin_amount) == 4, all(spin_amount >= 0 & spin_amount <= 1))
  }
  
  # map spin → lambda sharpness per wheel
  lambda_vec <- lambda_from_spin(spin_amount, lambda_sharp, lambda_flat)
  
  # build a length-10 prob vector for each wheel
  p1 <- position_prior_kernel(jumbles$digit_1, lambda_vec[1], row_weights)
  p2 <- position_prior_kernel(jumbles$digit_2, lambda_vec[2], row_weights)
  p3 <- position_prior_kernel(jumbles$digit_3, lambda_vec[3], row_weights)
  p4 <- position_prior_kernel(jumbles$digit_4, lambda_vec[4], row_weights)
  
  # return priors and the spins/lambdas we used (handy for logging)
  list(p1 = p1, p2 = p2, p3 = p3, p4 = p4,
       spin_amount = spin_amount, lambda_vec = lambda_vec)
}


# ============================================================
# scoring + optional rule for first digit
# ============================================================

# multiply the four wheel priors to score every row of combos_universe.
# you can optionally prohibit certain starting digits (hard 0 weight).
code_weights_from_priors <- function(priors, combos_universe, first_digit_prohibit = integer(0)) {
  M <- if (is.data.frame(combos_universe)) as.matrix(combos_universe[,1:4]) else combos_universe
  storage.mode(M) <- "integer"
  
  # base weight = product of per-wheel probabilities at each candidate’s digits
  w <- priors$p1[M[,1] + 1] *
    priors$p2[M[,2] + 1] *
    priors$p3[M[,3] + 1] *
    priors$p4[M[,4] + 1]
  
  # ban starts we don’t like (e.g., 7/8/9) by zeroing them out
  if (length(first_digit_prohibit) > 0) {
    w[M[,1] %in% first_digit_prohibit] <- 0
  }
  
  # normalize to sum 1 (avoid divide-by-zero if everything got banned)
  s <- sum(w)
  if (s == 0) stop("all weights zero after first-digit rule; relax the prohibition.")
  w / s
}


# ============================================================
# ranking + sequential filtering by guess/k feedback
# ============================================================

# rank_universe:
# - takes combos (matrix or df), weights, and an optional sequence of up to 5 guesses with k feedback.
# - for each guess/k, keeps only rows with exactly k positional matches vs that guess.
# - weights are subset + renormalized after every filter.
# - returns a ranked data.frame with digits, combo_str, numeric combo, and weight.
rank_universe <- function(combos, weights, guesses = NULL, ks = NULL) {
  # make sure we’re working with an integer matrix
  if (is.data.frame(combos)) combos <- as.matrix(combos[, 1:4])
  storage.mode(combos) <- "integer"
  n <- nrow(combos)
  
  # basic checks on weights
  stopifnot(length(weights) == n)
  w <- as.numeric(weights)
  
  # optional: apply each guess/k in order (for the game, we have max of 5)
  if (!is.null(guesses) || !is.null(ks)) {
    stopifnot(!is.null(guesses), !is.null(ks))
    stopifnot(is.list(guesses), length(guesses) == length(ks), length(guesses) <= 5)
    
    for (i in seq_along(guesses)) {
      g <- guesses[[i]]
      k <- ks[[i]]
      stopifnot(length(g) == 4, all(g %in% 0:9), k %in% 0:4)
      
      n_now <- nrow(combos)
      if (n_now == 0) break  # nothing left to do
      
      # build a matrix repeating the guess across rows, then count positional matches
      G <- matrix(rep(g, each = n_now), ncol = 4, byrow = FALSE)
      keep <- rowSums(combos == G) == k
      
      # filter combos and weights to consistent candidates
      combos <- combos[keep, , drop = FALSE]
      w <- w[keep]
      
      # if no rows left, bail with a friendly error (probably a bad k)
      if (length(w) == 0) {
        stop(sprintf("no candidates remain after guess #%d with k=%d. check inputs.", i, k))
      }
      
      # renormalize weights after each filter (so the next step is stable)
      w <- w / sum(w)
    }
  }
  
  # final renorm + ranking
  if (sum(w) <= 0) stop("all weights are zero. nothing to rank.")
  w <- w / sum(w)
  
  # build display fields
  combo_str   <- apply(combos, 1, paste0, collapse = "")
  combination <- combos[,1]*1000 + combos[,2]*100 + combos[,3]*10 + combos[,4]
  
  out <- data.frame(
    digit_1 = combos[,1],
    digit_2 = combos[,2],
    digit_3 = combos[,3],
    digit_4 = combos[,4],
    combo_str = combo_str,
    combination = combination,
    weight = w,
    stringsAsFactors = FALSE
  )
  
  # return sorted by weight (desc). ties are allowed
  out[order(out$weight, decreasing = TRUE), ]
}
