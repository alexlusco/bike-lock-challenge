# ============================ 01_build_guesses.R ============================

source("code/utils.R")

# ---------------- Load your jumbles ----------------
df <- exp_data # or real_data

# ---------------- Optional leak prior (build once) ----------------
# If you want leak available later, build the table now; keep beta=0 until you enable it.
# leak_prior_table <- build_leak_prior("data/4_pin_leaks.csv", temperature = 600)
leak_prior_table <- build_leak_prior(pin_data, temperature = 600)

# ---------------- Tunable params (edit anytime) -------------------
params <- list( # MP
  laplace        = 1,
  alpha_pmi      = 0.7, # 0.9 for
  w12            = 1.1,
  w34            = 1.25,
  w_cross        = 0.20,
  socs_bonus_log = log(1.7),
  leak_prior     = NULL,   # set to leak_prior_table when you want it ON
  beta_leak      = 0,      # e.g., 0.15 when enabling
  use_cross_pool = TRUE, K_prefix = 80, K_suffix = 80,
  use_mass_split = TRUE,
  sample_limit   = 800,
  respect_locks  = TRUE, lock_threshold = 0.88,
  always_offer_isolators = TRUE,
  top_n          = 20
)

params <- list( # AL
  laplace        = 1,
  alpha_pmi      = 0.6,          # strong PMI usage given coupling
  w12            = 1.1,          # max the prefix block (1,2)
  w34            = 1.25,          # max the suffix block (3,4)
  w_cross        = 0.2,          # keep a meaningful cross term
  socs_bonus_log = log(1.4),     # gentle suffix heuristic
  # leak: keep OFF initially; turn ON from round 3 if k is lagging
  #leak_prior     = leak_prior_table,  # build once; keep NULL if not loaded
  leak_prior     = NULL,
  #beta_leak      = 0.12,          # enable from Round 3+ only (set 0 for R1–R2)
  beta_leak      = 0,
  use_cross_pool = TRUE, K_prefix = 80, K_suffix = 80,
  use_mass_split = TRUE,
  sample_limit   = 900,           # you can bump to 1200 if you want a wider probe pool
  respect_locks  = TRUE, lock_threshold = 0.88,
  always_offer_isolators = TRUE,
  top_n          = 20
)

# ---------------- Helper printer -------------------
.code_str <- function(v) paste0(v, collapse = "")
print_options <- function(tag, res) {
  cat("\n", tag, "\n", sep = "")
  cat("Remaining:", res$remaining, "\n")
  if (any(!is.na(res$locks))) {
    lk <- ifelse(is.na(res$locks), ".", res$locks)
    cat("Locks (NA='.'):", paste(lk, collapse=" "), "\n")
  }
  cat("\nSuggested plays\n")
  cat("  SHRINK:   ", .code_str(res$next_guess), "\n", sep = "")
  if (!is.null(res$isolator_tests)) {
    cat("  ISOL A:   ", .code_str(res$isolator_tests$testA), "\n", sep = "")
    cat("  ISOL B:   ", .code_str(res$isolator_tests$testB), "\n", sep = "")
  } else cat("  Isolators: (none)\n")
  cat("  MAP:      ", .code_str(as.integer(res$top[1, 1:4])), "\n", sep = "")
  cat("\nTop MAP candidates:\n"); print(res$top)
}

# ---------------- State you append each round --------------
guesses <- list()
ks      <- integer(0)

# ============================ ROUND 1 — COMPUTE ============================
set.seed(1)
res <- do.call(propose_next_guess, c(list(df = df, guesses = guesses, ks = ks, seed = 1), params))
print_options("Round 1 — COMPUTE", res)

# ===== YOU PLAY Round 1: pick ONE of the following, then set k =====
#guesses[[1]] <- res$next_guess                         # SHRINK
guesses[[1]] <- as.integer(res$top[1, 1:4])            # MAP
# if (!is.null(res$isolator_tests)) guesses[[1]] <- as.integer(res$isolator_tests$testA)  # ISOL A
# if (!is.null(res$isolator_tests)) guesses[[1]] <- as.integer(res$isolator_tests$testB)  # ISOL B
# guesses[[1]] <- c(1,2,3,4)                             # CUSTOM
 ks[1]        <- 0                                      # returned k (0..4)

# ============================ ROUND 2 — COMPUTE ============================
# Example: enable leak starting here (optional)
# params$leak_prior <- leak_prior_table; params$beta_leak <- 0.15
set.seed(2)
res <- do.call(propose_next_guess, c(list(df = df, guesses = guesses, ks = ks, seed = 2), params))
print_options("Round 2 — COMPUTE", res)

# ===== YOU PLAY Round 2 =====
 guesses[[2]] <- res$next_guess
# guesses[[2]] <- as.integer(res$top[1, 1:4])
# if (!is.null(res$isolator_tests)) guesses[[2]] <- as.integer(res$isolator_tests$testA)
# if (!is.null(res$isolator_tests)) guesses[[2]] <- as.integer(res$isolator_tests$testB)
# guesses[[2]] <- c(1,2,3,4)
 ks[2]        <- 0

# ============================ ROUND 3 — COMPUTE ============================
set.seed(3)
res <- do.call(propose_next_guess, c(list(df = df, guesses = guesses, ks = ks, seed = 3), params))
print_options("Round 3 — COMPUTE", res)

# ===== YOU PLAY Round 3 =====
 guesses[[3]] <- res$next_guess
# guesses[[3]] <- as.integer(res$top[1, 1:4])
# if (!is.null(res$isolator_tests)) guesses[[3]] <- as.integer(res$isolator_tests$testA)
# if (!is.null(res$isolator_tests)) guesses[[3]] <- as.integer(res$isolator_tests$testB)
# guesses[[3]] <- c(1,2,3,4)
 ks[3]        <- 1

# ============================ ROUND 4 — COMPUTE ============================
set.seed(4)
res <- do.call(propose_next_guess, c(list(df = df, guesses = guesses, ks = ks, seed = 4), params))
print_options("Round 4 — COMPUTE", res)

# ===== YOU PLAY Round 4 =====
 guesses[[4]] <- res$next_guess
# guesses[[4]] <- as.integer(res$top[1, 1:4])
# if (!is.null(res$isolator_tests)) guesses[[4]] <- as.integer(res$isolator_tests$testA)
# if (!is.null(res$isolator_tests)) guesses[[4]] <- as.integer(res$isolator_tests$testB)
# guesses[[4]] <- c(1,2,3,4)
 ks[4]        <- 0

# ============================ ROUND 5 — COMPUTE ============================
set.seed(5)
res <- do.call(propose_next_guess, c(list(df = df, guesses = guesses, ks = ks, seed = 5), params))
print_options("Round 5 — COMPUTE", res)

# ===== YOU PLAY Round 5 =====
# guesses[[5]] <- res$next_guess
# guesses[[5]] <- as.integer(res$top[1, 1:4])
# if (!is.null(res$isolator_tests)) guesses[[5]] <- as.integer(res$isolator_tests$testA)
# if (!is.null(res$isolator_tests)) guesses[[5]] <- as.integer(res$isolator_tests$testB)
# guesses[[5]] <- c(1,2,3,4)
# ks[5]        <- 4
