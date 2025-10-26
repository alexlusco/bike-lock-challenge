# ------------------------------------------------------------
# author: alex luscombe & m.p.
# date: october 2025
#
# what this file does:
# - builds wheel-wise priors from observed "jumble" digits, using entropy to
#   infer per-wheel spin and a circular kernel to spread mass to nearby digits.
# - optionally calibrates kernel hyperparameters by rolling cv over a "week"
#   column to avoid leaking future weeks into past tuning.
# - constructs a base posterior over all 10,000 codes, with optional bans on
#   leading digits (e.g., to avoid rare starts from pin leak data).
# - optionally blends in behavioral priors: a full-code pin csv and a simple
#   pattern model (repeats/sequences/dates/years). exponents control influence.
# - supports two signals for picking the next guess:
#     (a) the map code (highest posterior mass)
#     (b) a shrink-maximizer that aims to reduce the survivor set the most
# - includes knobs for “early game” (explore/shrink hard) vs “later game”
#   (exploit higher-probability codes).
# ------------------------------------------------------------


# ============================================================
# 0) inputs and quick diagnostics
# ============================================================

#jumble_data <- exp_data
jumble_data <- real_data

# pin leak summary (optional): helps decide banned starting digits
# tweak:
# - method: "rank" downweights lower rows by r^(-alpha); "uniform" treats all rows equally
# - alpha (0.6–1.2 reasonable): higher → concentrates on the top of the leak list
lead_tab <- leading_digit_summary_from_pin_csv(
  "data/4_pin_leaks.csv",
  method = "rank",
  alpha  = 0.8
)
print(lead_tab)

# set any leading digits you want to prohibit (0–9). empty vector means "ban nothing".
# idea:
# - early game: consider banning very-rare starters to avoid low-yield shots
# - later game: lift bans if posterior mass starts to concentrate there
# example to ban the bottom-3 from the leak: least_common_start_digits("data/4_pin_leaks.csv", k=3)
prohibited_numbers <- c(9, 8, 7, 6)


# ============================================================
# 1) optional time-decay row weights (weekly arrivals)
# ============================================================

# row weights are used when building wheel histograms so newer rows matter more.
# tweak:
# - alpha in (0,1): geometric decay per step of "age"; closer to 1.0 = slower decay (remember more history)
# - early game: use stronger decay (e.g., 0.90–0.93) if weekly behavior drifts
# - later game: relax decay (e.g., 0.95–0.98) once behavior looks stable
# note: if you don't have a helper like row_weights_from_age, set rw <- NULL to disable time decay.
rw <- if (exists("row_weights_from_age")) {
  row_weights_from_age(jumble_data, alpha = 0.93)   # set to NULL if you don’t want decay
} else {
  NULL
}


# ============================================================
# 2) rolling cv grouping (temporal splits)
# ============================================================

# only use a column literally named "week" for rolling cv; otherwise no grouping.
# tweak:
# - keep group_col = "week" if you want pre/post temporal generalization checks
# - set do_cv = FALSE below if you want faster startup with default params
group_col <- if ("week" %in% names(jumble_data)) "week" else NULL


# ============================================================
# 3) build wheel priors (with optional cv calibration)
# ============================================================

# suppose you believe, on average, the wheels land ~this far from their modal stop:
# distances are per wheel (units: digits on the ring; 0..5; wrap-around aware)
dist_per_wheel <- c(2.6, 1.9, 2.3, 1.7)

# map distances -> spin in [0,1] (tune d_lo/d_hi if your notion of "little" vs "a lot" differs)
spin_override <- spin_from_dist(dist_per_wheel, d_lo = 1.5, d_hi = 2.7)

# build_priors_calibrated args (how to think about them):
# - jumbles: your data frame with columns digit_1..digit_4 (required).
# - do_cv (bool): true to grid-search lambda and entropy scaling via rolling cv.
#     early game: true is safer if you have enough weeks; false for speed.
# - group_col: set to "week" to do rolling splits; else random k-folds fallback.
# - row_weights: optional per-row weights (e.g., time decay). stronger decay if drift.
# - grid: custom tuning grid if you want finer search over (ent_min, ent_max, lambdas).
# - spin_* passthroughs: advanced knobs to override or bias spin estimation.
priors <- build_priors_calibrated(
  jumbles     = jumble_data,
  do_cv       = TRUE,          # early: true (if data allows); later: can set false for speed; if do_cv FALSE, and want to change lamba sharp and flat need to adjust in utils.R
  group_col   = group_col,     # only used if "week" exists
  row_weights = NULL,          # set to NULL to disable time decay
  spin_target = NULL,
  spin_weight = NULL,
  spin_groups = NULL,
)


# ============================================================
# 4) base posterior over all 10,000 codes (with optional bans)
# ============================================================

# code_weights_from_priors args:
# - priors: the 4 wheel priors list returned above.
# - combos_universe: the 10,000 codes matrix.
# - first_digit_prohibit: integer vector of digits to ban (e.g., c(9,8,6)).
#     note: banned codes are zeroed before normalization, so they cannot resurface later.
w_kernel <- code_weights_from_priors(
  priors,
  combos_universe,
  first_digit_prohibit = prohibited_numbers
)


# ============================================================
# 5) behavioral priors (toggleable)
# ============================================================

# 5a) full-code pin csv prior
# load_pin_prior_csv args:
# - csv_path: path to a csv with a 'pins' column (strings or ints).
# - alpha: converts rank r to score r^(-alpha). higher = more peaked toward top pins.
# - tail_extra: small boost for unseen pins vs a pure power tail.
# tweak:
# - early game: set eta_pin (below) higher (e.g., 0.6–1.0) if you strongly trust leaked pins.
# - later game: reduce eta_pin toward 0.2–0.4 to let your observed priors dominate.
use_pin_prior <- TRUE
pin_csv_path  <- "data/4_pin_leaks.csv"
eta_pin       <- 1.2 # this is the blend weight of PIN data into final code weights

pin_prior <- if (use_pin_prior) {
  load_pin_prior_csv(pin_csv_path, combos_universe, alpha = 0.8)
} else NULL

# 5b) simple pattern prior (bike-lock heuristics)
# pattern_prior_from_weights args:
# - combos_universe: all codes to score.
# - weights_named: named vector over features (aaaa, aabb, abab, abba, seq_inc/dec, mmdd, ddmm, years).
# tweak:
# - increase weights on patterns you believe users prefer (e.g., sequences/dates).
# - eta_pattern below scales the influence on the final posterior.
use_pattern_prior <- TRUE
eta_pattern       <- 0.3
pattern_weights   <- default_bikelock_pattern_weights()
pattern_prior     <- if (use_pattern_prior) {
  pattern_prior_from_weights(combos_universe, pattern_weights)
} else NULL


# ============================================================
# 6) blend priors multiplicatively (then normalize)
# ============================================================

# compose_code_weights args:
# - base_w: the kernel/posterior over codes (already normalized).
# - pin_p, patt_p: optional priors over full codes.
# - eta_pin, eta_patt: exponents controlling influence (0 disables).
#     intuition: log p_total = log base + eta_pin*log pin + eta_patt*log patt
# tweak:
# - early game: larger etas to inject stronger human-behavior signals.
# - later game: shrink etas toward 0 so the data-driven kernel dominates.
w_all <- compose_code_weights(
  base_w   = w_kernel,
  pin_p    = pin_prior,      eta_pin   = if (use_pin_prior)    eta_pin    else 0,
  patt_p   = pattern_prior,  eta_patt  = if (use_pattern_prior) eta_pattern else 0
)

# optional exploration temperature
# temper_weights arg:
# - rho in (0,1]: lower flattens the distribution (more exploration),
#                 1.0 leaves as-is, >1.0 (not typical here) would sharpen.
# tweak:
# - early game: try rho ~ 0.75–0.9 to diversify candidates.
# - later game: keep rho = 1 for faithful probabilities.
# w_all <- temper_weights(w_all, rho = 0.8)


# ============================================================
# 7) game state (guesses + feedback) and survivor ranking
# ============================================================

# match_mode choices:
# - "pos": k counts exact position matches (strict).
# - "any": k is mastermind-style anywhere matches (bag-of-digits).
# tweak:
# - use "any" when you only know how many digits overlap irrespective of position.
# - use "pos" when feedback is exact positional matches.
match_mode_now <- "pos"

# provide your played guesses as a list of integer vectors, and matching k’s
# example: guesses <- list(c(8,4,5,9), c(9,3,4,7)); ks <- c(0,2)
guesses <- list()
ks      <- c()

# rank_universe args:
# - combos_universe: candidate codes to evaluate.
# - w_all: weights over the full universe (after blending).
# - guesses, ks: optional sequential constraints; each step filters survivors and renormalizes.
# - match_mode: "any" or "pos" as above.
ranked <- rank_universe(
  combos_universe,
  w_all,
  guesses    = if (length(guesses)) guesses else NULL,
  ks         = if (length(ks))      ks      else NULL,
  match_mode = match_mode_now
)

cat("survivors:", nrow(ranked), "of 10,000\n")


# ============================================================
# 8) signal a: most-likely (map) code
# ============================================================

MAP <- unname(as.integer(ranked[1, c("digit_1","digit_2","digit_3","digit_4")]))
MAP_prob <- ranked$weight[1]
cat(sprintf("signal a — map guess: %s (pr=%.2f%%)\n",
            paste(MAP, collapse = ""), 100*MAP_prob))


# ============================================================
# 9) signal b: shrink-maximizer among a candidate slate
# ============================================================

# shrink_best_among args:
# - ranked_survivors: data frame from rank_universe (digits + weight).
# - candidate_pool: "topN", "topK", or "all" (controls slate size).
# - N, K: sizes for topN/topK.
# - objective: "expected_entropy" (information) or "expected_remaining" (mass-weighted survivors).
# - match_mode: must match your feedback semantics ("any" or "pos").
# - repeat_penalty: small penalty to prefer guesses with more distinct digits (useful under "any").
# tweak:
# - early game: objective="expected_entropy", candidate_pool="topN", N ~ 1000–3000, repeat_penalty ~ 0.03–0.08.
# - later game: reduce N (e.g., 500–1500), maybe set repeat_penalty → 0 to lean into p(correct).
shrink_table <- shrink_best_among(
  ranked_survivors = ranked,
  candidate_pool   = "topN",
  N                = 2000L,
  K                = 10L,
  objective        = "expected_entropy",
  match_mode       = match_mode_now,
  repeat_penalty   = if (match_mode_now == "any") 0.05 else 0.0
)

shrink_best <- unname(as.integer(shrink_table[1, c("digit_1","digit_2","digit_3","digit_4")]))
best_val    <- shrink_table$expected_value[1]
best_p_corr <- shrink_table$prob_correct[1]

cat(sprintf(
  "signal b — best shrink (topN): %s | expected=%.3f | pr(correct)=%.2f%%\n",
  paste(shrink_best, collapse=""), best_val, 100*best_p_corr
))


# ============================================================
# 10) policy: shoot map if confident, else shrink
# ============================================================

# recommend_next args:
# - priors: wheel priors (the kernel part).
# - guesses, ks, match_mode: same as rank_universe.
# - tau_map: probability threshold; if map ≥ tau_map, take it, else take shrink pick.
#     early game: higher (e.g., 0.20–0.30) to avoid premature shots.
#     later game: lower (e.g., 0.12–0.18) to capitalize on concentrated mass.
# - candidate_pool, N, K, objective, repeat_penalty: passed through to shrink scoring.
# - first_digit_prohibit: ensure bans are respected inside the policy as well.
policy <- recommend_next(
  priors              = priors,
  guesses             = if (length(guesses)) guesses else NULL,
  ks                  = if (length(ks))      ks      else NULL,
  match_mode          = match_mode_now,
  tau_map             = 0.20,
  candidate_pool      = "topN",
  N                   = 2000L,
  K                   = 10L,
  objective           = "expected_entropy",
  repeat_penalty      = if (match_mode_now == "any") 0.05 else 0.0,
  first_digit_prohibit= prohibited_numbers
)

cat("\n--- two-signal recommendation ---\n")
cat("map      :", paste(policy$signal_A$code, collapse=""),
    sprintf("(pr=%.2f%%)", 100*policy$signal_A$prob), "\n")
cat("shrink   :", paste(policy$signal_B$code, collapse=""),
    sprintf("| objective=%s", policy$signal_B$objective), "\n")
cat("pick     :", paste(policy$pick, collapse=""), "\n")


# ============================================================
# 11) quick diagnostics
# ============================================================

cat("\nTop 10 survivors (posterior):\n")
print(utils::head(policy$survivors$ranked, min(10, policy$survivors$n)))

cat("\nTop 10 shrink candidates (by objective):\n")
print(utils::head(policy$signal_B$table, 10))
