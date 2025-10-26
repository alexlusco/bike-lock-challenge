# ------------------------------------------------------------
# author: alex luscombe & m.p.
# date: october 2025
#
# what this file does:
# - builds wheel-wise priors from observed "jumble" stops using a circular kernel.
#   more concentrated stops → lower entropy → sharper kernel around popular digits.
# - scores all 10,000 4-digit codes by multiplying the four wheel priors.
# - optionally bans selected starting digits (first wheel) by zeroing their mass.
# - supports behavioral priors:
#     (a) a full-code prior from a ranked pin csv (from haveibeenpwned)
#     (b) a parametric pattern prior for bike locks (repeats, sequences, dates/years)
# - ranks candidates and applies sequential k-feedback filters (pos or anywhere).
# - includes fast shrink scoring utilities (expected remaining / expected entropy).
# - includes cv helpers to calibrate kernel sharpness and entropy mapping.
# ------------------------------------------------------------

# ============================================================
# misc + constants
# ============================================================

# digits used throughout (0..9)
digits <- 0:9 

# summarize prevalence of first digits in a ranked pin csv
# method "rank" downweights later rows as r^(-alpha); "uniform" treats all equal
leading_digit_summary_from_pin_csv <- function(csv_path,
                                               method = c("rank","uniform"),
                                               alpha = 0.8) {
  method <- match.arg(method)
  stopifnot(file.exists(csv_path))
  df <- utils::read.csv(csv_path, stringsAsFactors = FALSE)
  stopifnot("pins" %in% names(df))
  
  # clean to 4-digit strings (strip non-digits, left-pad with zeros)
  pins_chr <- gsub("\\D", "", df$pins)
  pins_chr <- sprintf("%04d", suppressWarnings(as.integer(pins_chr)))
  pins_chr <- pins_chr[!is.na(pins_chr)]
  if (!length(pins_chr)) stop("no 4-digit values found in 'pins' column.")
  
  # row weights derived from rank (1 = most common) or uniform
  w <- if (method == "rank") {
    r <- seq_along(pins_chr)
    r^(-alpha)
  } else {
    rep(1, length(pins_chr))
  }
  
  # aggregate by first digit → shares sum to 1
  lead <- as.integer(substr(pins_chr, 1, 1))
  agg  <- tapply(w, factor(lead, levels = 0:9), sum, default = 0)
  agg  <- as.numeric(agg)
  share <- agg / sum(agg)
  
  out <- data.frame(
    digit = 0:9,
    weight = agg,
    share = share,
    stringsAsFactors = FALSE
  )
  # smaller share = rarer start → lower rank value
  out$rank_least_common <- rank(out$share, ties.method = "min")
  out[order(out$share, out$digit), ]
}

# convenience wrapper to return the rarest starting digits
# either take the bottom-k by share or everything below a given share threshold
least_common_start_digits <- function(csv_path,
                                      k = 3,
                                      method = c("rank","uniform"),
                                      alpha = 0.8,
                                      max_share = NULL) {
  tab <- leading_digit_summary_from_pin_csv(csv_path, method = method[1], alpha = alpha)
  if (!is.null(max_share)) {
    return(tab$digit[tab$share <= max_share])
  }
  head(tab$digit, k)
}


# ============================================================
# geometry + information helpers
# ============================================================

# circular distance on a ring of 10 digits (0..9); max distance is 5
# used by the kernel so that 9 is close to 0 (wrap-around)
circ_dist <- function(a, b) {
  d <- abs(a - b)
  pmin(d, 10 - d)
}

# shannon entropy of a 10-bin histogram (natural log)
# lower = more peaky; used to map a wheel's "spin" to kernel sharpness
entropy <- function(tab10) {
  p <- tab10 / sum(tab10)
  -sum(ifelse(p > 0, p * log(p), 0))
}

# enumerate the full 10,000-code universe as a 4-column integer matrix
combos_universe <- as.matrix(
  expand.grid(digit_1 = digits, digit_2 = digits, digit_3 = digits, digit_4 = digits)
)

# small utilities
.clamp01 <- function(x) pmin(pmax(x, 0), 1)

# normalize a probability-like vector (fails if all mass is zero)
.normalize <- function(x) { s <- sum(x); if (s <= 0) stop("all weights zero"); x / s }

# partial pooling across positions (optional): shrink subsets toward their mean
.partial_pool <- function(x, groups, strength = 0) {
  if (is.null(groups) || strength <= 0) return(x)
  out <- x
  for (g in groups) {
    m <- mean(x[g])
    out[g] <- (1 - strength) * x[g] + strength * m
  }
  out
}

# map an observed distance metric to spin in [0,1]; not used by default
spin_from_dist <- function(d, d_lo = 1.5, d_hi = 2.7) {
  z <- (pmin(pmax(d, d_lo), d_hi) - d_lo) / (d_hi - d_lo)
  pmin(pmax(z, 0), 1)
}

# turn spin in [0,1] into a kernel sharpness (lambda)
# spin=0 → sharp (lambda_sharp); spin=1 → flat (lambda_flat)
lambda_from_spin <- function(spin, lambda_sharp = 1.0, lambda_flat = 0.5) {
  stopifnot(length(spin) == 4, all(spin >= 0 & spin <= 1))
  lambda_sharp * (1 - spin) + lambda_flat * spin
}

# estimate per-wheel spin via entropy of the observed stop histogram
# rescale entropy linearly into [0,1] between (ent_min, ent_max)
auto_spin_from_entropy <- function(jumbles, ent_min = 0.6, ent_max = 2.3) {
  stopifnot(all(c("digit_1","digit_2","digit_3","digit_4") %in% names(jumbles)))
  rescale01 <- function(x, xmin, xmax) {
    z <- (x - xmin) / (xmax - xmin)
    pmin(pmax(z, 0), 1)
  }
  ents <- vapply(1:4, function(i) {
    tab <- tabulate(jumbles[[paste0("digit_", i)]] + 1, nbins = 10)
    entropy(tab)
  }, numeric(1))
  rescale01(ents, ent_min, ent_max)
}


# ============================================================
# kernel priors (wheel-wise)
# ============================================================

# single-wheel prior using a symmetric circular kernel over digits 0..9
# hist counts are convolved with exp(-lambda * distance) out to distance 5
position_prior_kernel <- function(obs_digits, lambda, row_weights = NULL) {
  hist <- tabulate(obs_digits + 1, nbins = 10)
  if (!is.null(row_weights)) {
    # if per-row weights are given (e.g., time decay), aggregate with weights
    hist <- tapply(row_weights, factor(obs_digits, levels = 0:9), sum, default = 0)
    hist <- as.numeric(hist)
  }
  K <- exp(-lambda * (0:5))  # precompute kernel by distance
  scores <- sapply(digits, function(t) {
    s <- 0
    for (x in digits) {
      d <- circ_dist(t, x)
      s <- s + hist[x + 1] * K[d + 1]
    }
    s
  })
  scores / sum(scores)
}

# build 4 wheel priors by mapping entropy→spin→lambda and calling the kernel
build_priors <- function(jumbles,
                         spin_amount  = NULL,
                         lambda_sharp = 1.1,
                         lambda_flat  = 0.5,
                         row_weights  = NULL,
                         spin_target        = NULL,
                         spin_weight        = 0,
                         spin_groups        = NULL,
                         spin_group_strength= 0,
                         spin_multipliers   = NULL,
                         # pass-through for auto spin mapping
                         ent_min = 0.6,
                         ent_max = 2.3) {
  stopifnot(all(c("digit_1","digit_2","digit_3","digit_4") %in% names(jumbles)))
  
  # choose spin: either supplied or inferred from entropy
  if (is.null(spin_amount)) {
    spin <- auto_spin_from_entropy(jumbles, ent_min = ent_min, ent_max = ent_max)
  } else {
    stopifnot(length(spin_amount) == 4, all(spin_amount >= 0 & spin_amount <= 1))
    spin <- spin_amount
  }
  
  # optional: per-wheel multipliers and partial pooling toward group means
  if (!is.null(spin_multipliers)) {
    stopifnot(length(spin_multipliers) == 4, all(spin_multipliers >= 0))
    spin <- .clamp01(spin * spin_multipliers)
  }
  if (!is.null(spin_target) && spin_weight > 0) {
    stopifnot(length(spin_target) == 4, all(spin_target >= 0 & spin_target <= 1))
    stopifnot(spin_weight >= 0 & spin_weight <= 1)
    spin <- (1 - spin_weight) * spin + spin_weight * spin_target
    spin <- .clamp01(spin)
  }
  if (!is.null(spin_groups) && spin_group_strength > 0) {
    stopifnot(spin_group_strength >= 0 & spin_group_strength <= 1)
    spin <- .partial_pool(spin, spin_groups, strength = spin_group_strength)
    spin <- .clamp01(spin)
  }
  
  # convert spin to per-wheel sharpness
  lambda_vec <- lambda_from_spin(spin, lambda_sharp, lambda_flat)
  
  # build priors for each wheel
  p1 <- position_prior_kernel(jumbles$digit_1, lambda_vec[1], row_weights)
  p2 <- position_prior_kernel(jumbles$digit_2, lambda_vec[2], row_weights)
  p3 <- position_prior_kernel(jumbles$digit_3, lambda_vec[3], row_weights)
  p4 <- position_prior_kernel(jumbles$digit_4, lambda_vec[4], row_weights)
  
  list(p1 = p1, p2 = p2, p3 = p3, p4 = p4,
       spin_amount = spin, lambda_vec = lambda_vec,
       ent_min = ent_min, ent_max = ent_max)
}


# ============================================================
# code scoring (+ optional first-digit ban)
# ============================================================

# compute p(code) by multiplying wheel priors; optionally zero out disallowed starts
code_weights_from_priors <- function(priors, combos_universe, first_digit_prohibit = integer(0)) {
  M <- if (is.data.frame(combos_universe)) as.matrix(combos_universe[,1:4]) else combos_universe
  storage.mode(M) <- "integer"
  w <- priors$p1[M[,1] + 1] * priors$p2[M[,2] + 1] * priors$p3[M[,3] + 1] * priors$p4[M[,4] + 1]
  if (length(first_digit_prohibit) > 0) w[M[,1] %in% first_digit_prohibit] <- 0
  .normalize(w)
}

# optional temperature flattening for exploration in early game
temper_weights <- function(w, rho = 0.8) {
  w2 <- w^rho
  w2 / sum(w2)
}


# ============================================================
# behavioral priors (pin csv + bike-lock patterns)
# ============================================================

# full-code prior from a ranked pin csv (column 'pins' as strings)
# converts rank to score via r^{-alpha}, maps onto the 10k universe, fills tail uniformly
load_pin_prior_csv <- function(csv_path,
                               combos_universe,
                               alpha = 0.8,
                               tail_extra = 1.0) {
  stopifnot(file.exists(csv_path))
  df <- utils::read.csv(csv_path, stringsAsFactors = FALSE)
  stopifnot("pins" %in% names(df))
  
  # normalize to 4-digit strings
  pins_chr <- gsub("\\D", "", df$pins)
  pins_chr <- sprintf("%04d", as.integer(pins_chr))
  pins_chr <- pins_chr[!is.na(pins_chr)]
  
  # 1-based ranks → decreasing scores
  r <- seq_along(pins_chr)
  score <- r^(-alpha)
  
  # map scores into a full 10k vector keyed by "d1d2d3d4"
  key_all <- apply(combos_universe, 1, paste0, collapse = "")
  pin_p <- rep(NA_real_, length(key_all))
  names(pin_p) <- key_all
  
  pin_p[pins_chr] <- score[match(pins_chr, names(pin_p))]
  
  # for codes not present in csv, assign a small "tail" mass so nothing is zero
  if (anyNA(pin_p)) {
    tail_val <- ( (length(r) + 1)^(-alpha) ) * tail_extra
    pin_p[is.na(pin_p)] <- tail_val
  }
  .normalize(pin_p)
}

# feature engineering for bike-lock patterns (no keypad geometry)
# booleans for repeats, alternating patterns, wrap-around sequences, mmdd/ddmm, and years
pattern_features_matrix <- function(C) {
  if (is.data.frame(C)) C <- as.matrix(C[,1:4])
  storage.mode(C) <- "integer"
  d1 <- C[,1]; d2 <- C[,2]; d3 <- C[,3]; d4 <- C[,4]
  
  AAAA <- (d1==d2 & d2==d3 & d3==d4)
  AABB <- (d1==d2 & d3==d4 & d1!=d3)
  ABAB <- (d1==d3 & d2==d4 & d1!=d2)
  ABBA <- (d1==d4 & d2==d3 & d1!=d2)
  
  # wrap-around +1 and -1 sequences (e.g., 8-9-0-1)
  inc12 <- (d2 == (d1 + 1) %% 10)
  inc23 <- (d3 == (d2 + 1) %% 10)
  inc34 <- (d4 == (d3 + 1) %% 10)
  seq_inc <- (inc12 & inc23 & inc34)
  
  dec12 <- (d2 == (d1 + 9) %% 10)
  dec23 <- (d3 == (d2 + 9) %% 10)
  dec34 <- (d4 == (d3 + 9) %% 10)
  seq_dec <- (dec12 & dec23 & dec34)
  
  # simple calendar-like patterns
  mm   <- 10*d1 + d2
  dd   <- 10*d3 + d4
  dd2  <- 10*d1 + d2
  mm2  <- 10*d3 + d4
  
  is_mm  <- (mm  >= 1 & mm  <= 12)
  is_dd  <- (dd  >= 1 & dd  <= 31)
  is_dd2 <- (dd2 >= 1 & dd2 <= 31)
  is_mm2 <- (mm2 >= 1 & mm2 <= 12)
  
  mmdd <- (is_mm & is_dd)
  ddmm <- (is_dd2 & is_mm2)
  
  # two-digit year suffixes (approx ranges)
  last <- 10*d3 + d4
  year_19 <- (d1==1 & d2==9 & last >= 50 & last <= 99)
  year_20 <- (d1==2 & d2==0 & last >= 0  & last <= 24)
  
  cbind(AAAA=AAAA*1L, AABB=AABB*1L, ABAB=ABAB*1L, ABBA=ABBA*1L,
        seq_inc=seq_inc*1L, seq_dec=seq_dec*1L,
        mmdd=mmdd*1L, ddmm=ddmm*1L, year_19=year_19*1L, year_20=year_20*1L)
}

# default log-linear weights for the pattern prior (positive favors that pattern)
default_bikelock_pattern_weights <- function() {
  c(AAAA=1.5, AABB=0.8, ABAB=0.6, ABBA=0.6,
    seq_inc=1.2, seq_dec=1.0,
    mmdd=0.9, ddmm=0.6, year_19=0.4, year_20=0.3)
}

# convert features into a probability prior via p ∝ exp(w · x)
# safe for large w via max-subtraction
pattern_prior_from_weights <- function(C, weights_named) {
  X <- pattern_features_matrix(C)
  # align weight vector to feature columns; missing weights default to 0
  w <- rep(0, ncol(X)); names(w) <- colnames(X)
  for (nm in names(weights_named)) if (nm %in% colnames(X)) w[nm] <- weights_named[[nm]]
  s <- as.numeric(X %*% w)
  p <- exp(s - max(s))
  .normalize(p)
}

# blend multiple priors multiplicatively (log-space addition)
# exponents let you dial influence without editing the priors themselves
compose_code_weights <- function(base_w,
                                 pin_p   = NULL, eta_pin = 0,
                                 patt_p  = NULL, eta_patt = 0,
                                 eps = 1e-300) {
  w <- base_w
  if (!is.null(pin_p)  && eta_pin  != 0)  w <- w * (pmax(pin_p,  eps) ^ eta_pin)
  if (!is.null(patt_p) && eta_patt != 0)  w <- w * (pmax(patt_p, eps) ^ eta_patt)
  .normalize(w)
}

# ============================================================
# ranking + sequential filtering by feedback
# ============================================================

# rank all surviving codes given weights and optional guess/feedback constraints
# match_mode "pos" = exact slot matches; "any" = mastermind-style anywhere matches
rank_universe <- function(combos, weights, guesses = NULL, ks = NULL,
                          match_mode = c("any", "pos")) {
  match_mode <- match.arg(match_mode)
  
  if (is.data.frame(combos)) combos <- as.matrix(combos[, 1:4])
  storage.mode(combos) <- "integer"
  
  n <- nrow(combos)
  stopifnot(length(weights) == n)
  w <- as.numeric(weights)
  
  # exact-position matches for a given guess g
  .k_pos_for_guess <- function(C, g) {
    n <- nrow(C)
    G <- matrix(rep(as.integer(g), each = n), ncol = 4, byrow = FALSE)
    rowSums(C == G)
  }
  # anywhere (bag-of-digits) matches for a given guess g
  .k_any_for_guess <- function(C, g) {
    n <- nrow(C)
    counts_mat <- matrix(0L, nrow = n, ncol = 10)
    for (d in 0:9) counts_mat[, d + 1] <- rowSums(C == d)
    gcount <- tabulate(g + 1, 10)
    rowSums(pmin(counts_mat, matrix(rep(gcount, each = n), nrow = n)))
  }
  
  # apply sequential filters if guesses/ks provided; renormalize after each filter
  if (!is.null(guesses) || !is.null(ks)) {
    stopifnot(!is.null(guesses), !is.null(ks))
    stopifnot(is.list(guesses), length(guesses) == length(ks), length(guesses) <= 5)
    
    for (i in seq_along(guesses)) {
      g <- as.integer(guesses[[i]]); k <- as.integer(ks[[i]])
      stopifnot(length(g) == 4L, all(g %in% 0:9), k %in% 0:4)
      if (nrow(combos) == 0) break
      kvec <- if (match_mode == "any") .k_any_for_guess(combos, g) else .k_pos_for_guess(combos, g)
      keep <- (kvec == k)
      combos <- combos[keep, , drop = FALSE]
      w <- w[keep]
      if (!length(w)) stop(sprintf("no candidates remain after guess #%d with k=%d. check inputs.", i, k))
      w <- w / sum(w)
    }
  }
  
  if (sum(w) <= 0) stop("all weights are zero. nothing to rank.")
  w <- w / sum(w)
  
  combo_str   <- apply(combos, 1, paste0, collapse = "")
  combination <- combos[, 1] * 1000 + combos[, 2] * 100 + combos[, 3] * 10 + combos[, 4]
  
  out <- data.frame(
    digit_1 = combos[, 1],
    digit_2 = combos[, 2],
    digit_3 = combos[, 3],
    digit_4 = combos[, 4],
    combo_str = combo_str,
    combination = combination,
    weight = w,
    stringsAsFactors = FALSE
  )
  out[order(out$weight, decreasing = TRUE), ]
}

# ============================================================
# fast shrink scoring + entropy objective (any/pos)
# ============================================================

# precompute per-survivor digit counts for the anywhere rule
row_digit_counts_10 <- function(survivors) {
  n <- nrow(survivors)
  C <- matrix(0L, n, 10)
  for (d in 0:9) C[, d + 1] <- rowSums(survivors == d)
  C
}

# collect counts/mass of k under the anywhere rule for a candidate guess g
k_counts_mass_any <- function(C10, surv_w, g) {
  g10 <- tabulate(g + 1, 10)
  kvec <- rowSums(pmin(C10, matrix(g10, nrow(C10), 10, byrow = TRUE)))
  counts <- tabulate(kvec + 1, nbins = 5)
  mass   <- tapply(surv_w, factor(kvec, levels = 0:4), sum)
  mass[is.na(mass)] <- 0
  list(counts_k = counts, mass_k = as.numeric(mass), kvec = kvec)
}

# same as above but for exact-position matches
k_counts_mass_pos <- function(surv_combos, surv_w, g) {
  n <- nrow(surv_combos)
  G <- matrix(rep(as.integer(g), each = n), ncol = 4, byrow = FALSE)
  kvec <- rowSums(surv_combos == G)
  counts <- tabulate(kvec + 1, nbins = 5)
  mass   <- tapply(surv_w, factor(kvec, levels = 0:4), sum)
  mass[is.na(mass)] <- 0
  list(counts_k = counts, mass_k = as.numeric(mass), kvec = kvec)
}

# expected remaining survivors after playing g (any / pos)
expected_remaining_any <- function(C10, surv_w, g) {
  km <- k_counts_mass_any(C10, surv_w, g)
  sum(km$mass_k * km$counts_k)
}
expected_remaining_pos <- function(surv_combos, surv_w, g) {
  km <- k_counts_mass_pos(surv_combos, surv_w, g)
  sum(km$mass_k * km$counts_k)
}

# expected posterior entropy after playing g (lower is better)
expected_entropy_any <- function(C10, surv_w, g) {
  km <- k_counts_mass_any(C10, surv_w, g)
  Hk <- numeric(5)
  for (k in 0:4) {
    idx <- which(km$kvec == k)
    if (!length(idx)) { Hk[k+1] <- 0; next }
    wk <- surv_w[idx] / sum(surv_w[idx])
    Hk[k+1] <- -sum(ifelse(wk > 0, wk * log(wk), 0))
  }
  sum(km$mass_k * Hk)
}
expected_entropy_pos <- function(surv_combos, surv_w, g) {
  km <- k_counts_mass_pos(surv_combos, surv_w, g)
  Hk <- numeric(5)
  for (k in 0:4) {
    idx <- which(km$kvec == k)
    if (!length(idx)) { Hk[k+1] <- 0; next }
    wk <- surv_w[idx] / sum(surv_w[idx])
    Hk[k+1] <- -sum(ifelse(wk > 0, wk * log(wk), 0))
  }
  sum(km$mass_k * Hk)
}

# evaluate shrink quality for a slate of candidates; supports topk/topn/all pools
# repeat_penalty > 0 slightly prefers distinct-digit guesses under the anywhere rule
shrink_best_among <- function(ranked_survivors,
                              candidate_pool = c("topK","topN","all"),
                              K = 10L, N = 2000L,
                              objective = c("expected_remaining","expected_entropy"),
                              match_mode = c("any","pos"),
                              repeat_penalty = 0.0) {
  candidate_pool <- match.arg(candidate_pool)
  objective <- match.arg(objective)
  match_mode <- match.arg(match_mode)
  
  stopifnot(all(c("digit_1","digit_2","digit_3","digit_4","weight") %in% names(ranked_survivors)))
  
  S <- as.matrix(ranked_survivors[, c("digit_1","digit_2","digit_3","digit_4")])
  w <- ranked_survivors$weight
  Nsurv <- nrow(S)
  
  # choose candidate set
  cand <- switch(candidate_pool,
                 topK = S[seq_len(min(K, Nsurv)), , drop = FALSE],
                 topN = S[seq_len(min(N, Nsurv)), , drop = FALSE],
                 all  = combos_universe
  )
  
  # pick scoring function based on mode and objective
  if (match_mode == "any") {
    C10 <- row_digit_counts_10(S)
    scorer <- switch(objective,
                     expected_remaining = function(g) expected_remaining_any(C10, w, g),
                     expected_entropy   = function(g) expected_entropy_any(C10, w, g)
    )
  } else {
    scorer <- switch(objective,
                     expected_remaining = function(g) expected_remaining_pos(S, w, g),
                     expected_entropy   = function(g) expected_entropy_pos(S, w, g)
    )
  }
  
  # evaluate each candidate; also track its posterior p(correct) if it's a survivor
  M <- nrow(cand)
  exp_val <- numeric(M)
  p_corr  <- numeric(M)
  rep_pen <- numeric(M)
  
  key <- function(M) apply(M, 1, paste0, collapse = "")
  wmap <- stats::setNames(w, key(S))
  
  for (i in seq_len(M)) {
    g <- cand[i, ]
    exp_val[i] <- scorer(g)
    p_corr[i]  <- wmap[key(matrix(g, nrow=1))] %||% 0
    rep_pen[i] <- repeat_penalty * (4 - length(unique(g)))
  }
  
  df <- data.frame(
    digit_1 = cand[,1], digit_2 = cand[,2], digit_3 = cand[,3], digit_4 = cand[,4],
    expected_value = exp_val,
    prob_correct   = p_corr,
    repeats        = 4 - apply(cand, 1, function(x) length(unique(x))),
    stringsAsFactors = FALSE
  )
  df$combination <- df$digit_1*1000 + df$digit_2*100 + df$digit_3*10 + df$digit_4
  
  # primary sort by (objective + repeat penalty), tie-break by prob_correct then lexicographic
  ord <- order(df$expected_value + rep_pen, -df$prob_correct, df$combination)
  df[ord, , drop = FALSE]
}

# backwards-compatible helper that scores only the top-k survivors
shrink_best_among_topK <- function(ranked_survivors,
                                   topK = 10L,
                                   method = c("analytic","mc"),
                                   B = 2000L) {
  shrink_best_among(ranked_survivors,
                    candidate_pool = "topK",
                    K = topK,
                    objective = "expected_remaining",
                    match_mode = "any",
                    repeat_penalty = 0.0)
}

# ============================================================
# two-signal policy wrapper (map vs shrink)
# ============================================================

# returns both the most probable code (map) and the best shrink candidate, then picks one
recommend_next <- function(priors,
                           guesses = NULL, ks = NULL,
                           match_mode = "any",
                           tau_map = 0.20,
                           candidate_pool = "topN",
                           N = 2000L, K = 10L,
                           objective = "expected_entropy",
                           repeat_penalty = 0.05,
                           first_digit_prohibit = integer(0)) {
  
  # build overall posterior over all codes (respecting any banned starts)
  w_all <- code_weights_from_priors(priors, combos_universe, first_digit_prohibit)
  
  # apply any guess/feedback constraints and rank survivors
  ranked <- rank_universe(combos_universe, w_all, guesses, ks, match_mode = match_mode)
  
  # signal a: map (top posterior mass)
  MAP <- unname(as.integer(ranked[1, c("digit_1","digit_2","digit_3","digit_4")]))
  MAP_p <- ranked$weight[1]
  
  # signal b: best shrink among a candidate slate
  shrink_tab <- shrink_best_among(
    ranked_survivors = ranked,
    candidate_pool   = candidate_pool,
    K = K, N = N,
    objective = objective,
    match_mode = match_mode,
    repeat_penalty = repeat_penalty
  )
  SHR <- unname(as.integer(shrink_tab[1, c("digit_1","digit_2","digit_3","digit_4")]))
  
  # policy: if map prob is high enough, shoot; else take shrink pick
  pick <- if (MAP_p >= tau_map) MAP else SHR
  
  list(
    survivors = list(n = nrow(ranked), ranked = ranked),
    signal_A  = list(code = MAP, prob = MAP_p),
    signal_B  = list(code = SHR, table = shrink_tab, objective = objective),
    pick      = pick
  )
}


# ============================================================
# cv calibration (entropy mapping + kernel lambdas)
# ============================================================

# log-likelihood of held-out digits under wheel-wise priors (unsupervised)
.loglik_wheels <- function(test_df, priors) {
  ll <- 0
  for (j in 1:4) {
    p <- priors[[paste0("p", j)]]
    idx <- test_df[[paste0("digit_", j)]] + 1
    ll <- ll + sum(log(p[idx] + 1e-12))
  }
  ll
}

# helper to build priors from a parameter row
.build_with_params <- function(train_df, pars, row_weights = NULL) {
  build_priors(
    jumbles      = train_df,
    row_weights  = row_weights,
    lambda_sharp = pars$lambda_sharp,
    lambda_flat  = pars$lambda_flat,
    ent_min      = pars$ent_min,
    ent_max      = pars$ent_max
  )
}

# rolling cv across an ordered grouping column (e.g., "week")
# default fallback is random 5-fold when no group is provided
cv_calibrate_kernels <- function(all_data,
                                 group_col = NULL,
                                 grid = NULL,
                                 min_train_groups = 1,
                                 verbose = FALSE) {
  stopifnot(is.data.frame(all_data))
  if (is.null(grid)) {
    grid <- expand.grid(
      ent_min      = c(0.5, 0.6, 0.7),
      ent_max      = c(2.0, 2.3, 2.6),
      lambda_sharp = c(1.0, 1.1, 1.3),
      lambda_flat  = c(0.4, 0.5, 0.6),
      KEEP.OUT.ATTRS = FALSE
    )
  }
  
  # build rolling splits if group_col exists; else random folds
  if (!is.null(group_col) && group_col %in% names(all_data)) {
    g <- all_data[[group_col]]
    ug <- sort(unique(g))   # assumes group_col increases over time
    splits <- list()
    for (i in seq_along(ug)) {
      if (i == 1) next
      train_groups <- ug[seq_len(i-1)]
      test_group   <- ug[i]
      train_idx <- g %in% train_groups
      test_idx  <- g == test_group
      if (sum(train_idx) == 0 || sum(test_idx) == 0) next
      splits[[length(splits)+1]] <- list(train = train_idx, test = test_idx)
    }
  } else {
    set.seed(1)
    k <- 5L
    fold <- sample(rep(1:k, length.out = nrow(all_data)))
    splits <- lapply(1:k, function(i) list(train = fold != i, test = fold == i))
  }
  
  scores <- numeric(nrow(grid))
  counts <- integer(nrow(grid))
  
  # evaluate each grid row across splits and average
  for (s in splits) {
    train <- all_data[s$train, , drop = FALSE]
    test  <- all_data[s$test,  , drop = FALSE]
    if (nrow(train) < min_train_groups || nrow(test) == 0) next
    
    for (i in seq_len(nrow(grid))) {
      pars <- grid[i, , drop = FALSE]
      pri  <- .build_with_params(train, pars)
      ll   <- .loglik_wheels(test, pri)
      scores[i] <- scores[i] + ll
      counts[i] <- counts[i] + 1L
    }
  }
  
  avg <- ifelse(counts > 0, scores / counts, -Inf)
  res <- cbind(grid, mean_loglik = avg, n_splits = counts)
  res <- res[order(-res$mean_loglik, -res$n_splits), , drop = FALSE]
  rownames(res) <- NULL
  res
}

# convenience: run cv (optional) and then build priors using the best params
build_priors_calibrated <- function(jumbles,
                                    do_cv = TRUE,
                                    group_col = NULL,
                                    grid = NULL,
                                    row_weights = NULL,
                                    spin_amount = NULL,
                                    spin_target = NULL, spin_weight = 0,
                                    spin_groups = NULL, spin_group_strength = 0,
                                    spin_multipliers = NULL) {
  if (do_cv) {
    cvtab <- cv_calibrate_kernels(jumbles, group_col = group_col, grid = grid)
    best  <- cvtab[1, ]
    pri <- build_priors(
      jumbles      = jumbles,
      row_weights  = row_weights,
      lambda_sharp = best$lambda_sharp,
      lambda_flat  = best$lambda_flat,
      ent_min      = best$ent_min,
      ent_max      = best$ent_max,
      spin_amount  = spin_amount,
      spin_target  = spin_target,
      spin_weight  = spin_weight,
      spin_groups  = spin_groups,
      spin_group_strength = spin_group_strength,
      spin_multipliers    = spin_multipliers
    )
    attr(pri, "cv_table") <- cvtab
    pri
  } else {
    build_priors(
      jumbles      = jumbles,
      row_weights  = row_weights,
      lambda_sharp = 1.1,
      lambda_flat  = 0.5
    )
  }
}

# safe infix helper for null coalescing
`%||%` <- function(a, b) if (is.null(a)) b else a