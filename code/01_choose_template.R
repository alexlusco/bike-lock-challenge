## ---------------- choose hyperparameter columns ----------------
hyper_cols <- c(
  "param_laplace",
  "param_alpha_pmi",
  "param_w12",
  "param_w34",
  "param_w_cross",
  "param_socs_bonus_log",
  "param_leak_on",
  "param_beta_leak",
  "param_use_cross_pool",
  "param_K_prefix",
  "param_K_suffix",
  "param_use_mass_split",
  "param_sample_limit",
  "param_respect_locks",
  "param_lock_threshold",
  "param_always_offer_isolators",
  "param_top_n"
)

## ---------------- choose feature columns ----------------
# use only stat_* columns that appear in BOTH files
common_cols <- intersect(names(templates), names(game))
feature_cols <- common_cols[grepl("^stat_", common_cols)]

if (length(feature_cols) == 0) {
  stop("no stat_* feature columns found in both templates and game.")
}

cat("using feature columns:\n  ",
    paste(feature_cols, collapse = ", "), "\n\n")

X_train <- as.matrix(templates[, feature_cols, drop = FALSE])
X_game  <- as.matrix(game[,      feature_cols, drop = FALSE])

## ---------------- standardize on training stats ----------------
mu <- colMeans(X_train, na.rm = TRUE)
sd <- apply(X_train, 2, sd, na.rm = TRUE)
sd[sd == 0] <- 1

X_train_sc <- scale(X_train, center = mu, scale = sd)
X_game_sc  <- scale(X_game,  center = mu, scale = sd)

## ---------------- helper: 1-nn prediction ----------------
nn_predict_one <- function(x_row, X_train_sc, y_train) {
  diff <- sweep(X_train_sc, 2, x_row, FUN = "-")
  d2   <- rowSums(diff * diff)
  idx  <- which.min(d2)
  list(
    template_id = y_train[idx],
    dist        = sqrt(d2[idx])
  )
}

y_train <- templates$template_id

## ---------------- optional: "weirdness" threshold ----------------
compute_loocv_nn_dist <- function(X) {
  n <- nrow(X)
  d_min <- numeric(n)
  for (i in seq_len(n)) {
    diff <- sweep(X, 2, X[i, ], FUN = "-")
    d2   <- rowSums(diff * diff)
    d2[i] <- Inf
    d_min[i] <- sqrt(min(d2))
  }
  d_min
}

train_nn_dist <- compute_loocv_nn_dist(X_train_sc)
thresh_hi     <- quantile(train_nn_dist, 0.99, na.rm = TRUE)

cat("training NN distance summary:\n")
print(summary(train_nn_dist))
cat("99th percentile distance:", as.numeric(thresh_hi), "\n\n")

## ---------------- build template -> param lookup ----------------
param_cols <- intersect(hyper_cols, names(templates))
if (length(param_cols) == 0) {
  stop("no hyperparameter columns found in templates.")
}

template_params <- unique(templates[, c("template_id", param_cols), drop = FALSE])
template_params <- template_params[!duplicated(template_params$template_id), ]

## ---------------- apply to game data ----------------
game$template_id_pred   <- NA_character_
game$nn_dist            <- NA_real_
game$nn_is_weird_99pct  <- NA

for (i in seq_len(nrow(X_game_sc))) {
  res <- nn_predict_one(X_game_sc[i, ], X_train_sc, y_train)
  game$template_id_pred[i]  <- res$template_id
  game$nn_dist[i]           <- res$dist
  game$nn_is_weird_99pct[i] <- res$dist > thresh_hi
}

## merge params onto game rows
game_with_params <- merge(
  game,
  template_params,
  by.x = "template_id_pred",
  by.y = "template_id",
  all.x = TRUE,
  sort  = FALSE
)

## ---------------- summarize + cat() results to console ----------------

# How many rows got mapped to each template
cat("Summary of NN template assignments for game data:\n")
print(table(game_with_params$template_id_pred))
cat("\n")

# Optionally summarize distances
cat("NN distance summary over all game rows:\n")
print(summary(game_with_params$nn_dist))
cat("99th percentile training distance threshold:",
    sprintf("%.4f", thresh_hi), "\n")
cat("Fraction of game rows flagged 'weird' (dist > threshold):",
    mean(game_with_params$nn_is_weird_99pct, na.rm = TRUE), "\n\n")

# Take the most common predicted template and show its hyperparameters once
template_counts <- table(game_with_params$template_id_pred)
best_template   <- names(template_counts)[which.max(template_counts)]

cat("Chosen template for this game (most frequent 1-NN):", best_template, "\n")
cat("Median NN distance for rows with this template:",
    sprintf("%.4f",
            median(game_with_params$nn_dist[
              game_with_params$template_id_pred == best_template
            ], na.rm = TRUE)),
    "\n\n")

cat("Hyperparameters for template", best_template, ":\n")

# pick one representative row for that template
row_best <- game_with_params[
  match(best_template, game_with_params$template_id_pred),
]

for (p in param_cols) {
  val <- row_best[[p]]
  if (is.logical(val)) {
    cat("  ", p, "=", if (isTRUE(val)) "TRUE" else "FALSE", "\n", sep = "")
  } else {
    cat("  ", p, "=", as.character(val), "\n", sep = "")
  }
}
cat("\n")
