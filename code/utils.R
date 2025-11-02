# ============================== utils.R ==============================

# -------- basics --------
.standardize_guess <- function(g) {
  if (is.character(g) && length(g) == 1L) as.integer(strsplit(g, "")[[1L]]) else as.integer(g)
}
.code_str <- function(v) paste0(v, collapse = "")
.circ_mod <- function(x) ((x %% 10) + 10) %% 10
.enumerate_all_codes <- function() {
  expand.grid(digit_1=0:9, digit_2=0:9, digit_3=0:9, digit_4=0:9,
              KEEP.OUT.ATTRS=FALSE, stringsAsFactors=FALSE)
}
.matches_vec <- function(codes_mat, guess_vec) rowSums(sweep(codes_mat, 2, guess_vec, `==`))
.safe_log <- function(x) log(pmax(x, 1e-12))

# ======================= Leak Insights Builder =======================
.zpad4 <- function(x) sprintf("%04d", as.integer(gsub("\\D", "", as.character(x))))
.to_digits <- function(pin_chr) {
  v <- as.integer(strsplit(pin_chr, "")[[1]])
  names(v) <- paste0("digit_", 1:4); v
}

# Build a 10k prior table: digit_1..digit_4 + logp
build_leak_prior <- function(leak_input, temperature = 600) {
  leak_df <- if (is.character(leak_input)) read.csv(leak_input, stringsAsFactors = FALSE) else leak_input
  stopifnot("pins" %in% names(leak_df))
  pins <- .zpad4(leak_df$pins)
  rank <- seq_along(pins)                 # 1 = most popular
  r0   <- rank - min(rank)
  w    <- exp(-r0 / max(temperature, 1e-6))
  di   <- do.call(rbind, lapply(pins, .to_digits))
  m    <- cbind(as.data.frame(di), w = w)
  agg  <- aggregate(w ~ digit_1 + digit_2 + digit_3 + digit_4, data = m, FUN = sum)
  p    <- agg$w / sum(agg$w)
  agg$logp <- .safe_log(p)
  agg[, c("digit_1","digit_2","digit_3","digit_4","logp")]
}

# Position marginals from prior (force 0..9 levels)
.marginals_from_prior <- function(prior_df) {
  out <- vector("list", 4)
  for (i in 1:4) {
    col <- paste0("digit_", i)
    p <- tapply(exp(prior_df$logp), factor(prior_df[[col]], levels = 0:9), sum)
    p[is.na(p)] <- 0
    p <- p / sum(p)
    out[[i]] <- p
  }
  names(out) <- paste0("pos_", 1:4)
  out
}

# 10x10 joint + PMI
.joint_and_pmi <- function(prior_df, i, j) {
  ci <- paste0("digit_", i); cj <- paste0("digit_", j)
  df2 <- data.frame(
    w  = exp(prior_df$logp),
    di = factor(prior_df[[ci]], levels = 0:9),
    dj = factor(prior_df[[cj]], levels = 0:9)
  )
  Pij <- xtabs(w ~ di + dj, data = df2)  # 10x10
  Pij <- Pij / sum(Pij)
  Pi  <- rowSums(Pij); Pj <- colSums(Pij)
  ratio <- Pij / pmax(outer(Pi, Pj), 1e-12)
  pmi <- log(pmax(ratio, 1e-12))
  list(Pij = Pij, PMI = pmi, Pi = Pi, Pj = Pj)
}

.pattern_stats <- function(prior_df) {
  df <- prior_df
  w  <- exp(df$logp)
  AAAA <- df$digit_1==df$digit_2 & df$digit_2==df$digit_3 & df$digit_3==df$digit_4
  AABB <- df$digit_1==df$digit_2 & df$digit_3==df$digit_4 & df$digit_1!=df$digit_3
  ABAB <- df$digit_1==df$digit_3 & df$digit_2==df$digit_4 & df$digit_1!=df$digit_2
  ABBA <- df$digit_1==df$digit_4 & df$digit_2==df$digit_3 & df$digit_1!=df$digit_2
  pal  <- df$digit_1==df$digit_4 & df$digit_2==df$digit_3
  inc  <- (df$digit_2== (df$digit_1+1)%%10) & (df$digit_3== (df$digit_2+1)%%10) & (df$digit_4== (df$digit_3+1)%%10)
  dec  <- (df$digit_2== (df$digit_1+9)%%10) & (df$digit_3== (df$digit_2+9)%%10) & (df$digit_4== (df$digit_3+9)%%10)
  year_like <- (df$digit_1==1 & df$digit_2 %in% 8:9) | (df$digit_1==2 & df$digit_2==0)
  pr <- function(mask) sum(w[mask])/sum(w)
  list(
    AAAA      = pr(AAAA),
    AABB      = pr(AABB),
    ABAB      = pr(ABAB),
    ABBA      = pr(ABBA),
    palindrome= pr(pal),
    inc_seq   = pr(inc),
    dec_seq   = pr(dec),
    year_like = pr(year_like)
  )
}

leak_insights <- function(leak_prior_df) {
  pos_m <- .marginals_from_prior(leak_prior_df)
  start_bias <- pos_m[[1]]
  pref <- .joint_and_pmi(leak_prior_df, 1, 2)
  suff <- .joint_and_pmi(leak_prior_df, 3, 4)
  coupling_12 <- mean(abs(pref$PMI))
  coupling_34 <- mean(abs(suff$PMI))
  patt <- .pattern_stats(leak_prior_df)
  w12_suggest <- 1.0 + min(0.5, coupling_12 / 2)
  w34_suggest <- 1.0 + min(0.5, coupling_34 / 2)
  beta_suggest <- 0.0
  list(
    pos_marginals = pos_m,
    start_bias    = start_bias,
    prefix = list(P = pref$Pij, PMI = pref$PMI),
    suffix = list(P = suff$Pij, PMI = suff$PMI),
    coupling_strength = c("1,2" = coupling_12, "3,4" = coupling_34),
    patterns = patt,
    reco = list(w12 = w12_suggest, w34 = w34_suggest, beta_leak = beta_suggest)
  )
}

# ======================= Model scoring pieces =======================
build_priors_nb <- function(df, laplace = 1) {
  req <- c("digit_1","digit_2","digit_3","digit_4")
  stopifnot(all(req %in% names(df)))
  for (c in req) {
    df[[c]] <- as.integer(df[[c]])
    if (any(df[[c]] < 0 | df[[c]] > 9)) stop("digits must be 0..9")
  }
  priors <- lapply(req, function(col) {
    tab <- table(factor(df[[col]], levels = 0:9))
    p <- (as.numeric(tab) + laplace) / (sum(tab) + laplace*10)
    names(p) <- as.character(0:9); p
  })
  names(priors) <- paste0("p", 1:4)
  priors
}

build_pairwise_logphi <- function(df, alpha = 0.6, laplace = 1) {
  cols <- c("digit_1","digit_2","digit_3","digit_4")
  Pi <- lapply(cols, function(cn) {
    tab <- table(factor(df[[cn]], levels=0:9))
    (as.numeric(tab) + laplace) / (sum(tab) + 10*laplace)
  })
  names(Pi) <- cols
  pairs <- list(c(1,2), c(1,3), c(1,4), c(2,3), c(2,4), c(3,4))
  logphi <- vector("list", length(pairs))
  names(logphi) <- sapply(pairs, \(p) paste0("p",p[1],"p",p[2]))
  for (k in seq_along(pairs)) {
    i <- pairs[[k]][1]; j <- pairs[[k]][2]
    tab <- table(factor(df[[cols[i]]], levels=0:9),
                 factor(df[[cols[j]]], levels=0:9))
    Pij <- (tab + laplace) / (sum(tab) + 100*laplace)
    outer_PiPj <- outer(Pi[[i]], Pi[[j]])
    ratio <- Pij / pmax(outer_PiPj, 1e-12)
    logphi[[k]] <- alpha * log(pmax(ratio, 1e-12))
  }
  list(logphi = logphi, pairs = pairs)
}

modes_per_wheel <- function(df) {
  sapply(c("digit_1","digit_2","digit_3","digit_4"), function(col) {
    vc <- table(df[[col]]); as.integer(names(which.max(vc)))
  })
}
.socs_suffix_bonus <- function(code_vec, df, S = -4:4, U = -3:3, bonus_log = log(1.5)) {
  m <- modes_per_wheel(df)
  ok3 <- any(((m[3] - S) %% 10) == code_vec[3])
  ok4 <- any(((m[4] - U) %% 10) == code_vec[4])
  if (ok3 && ok4) bonus_log else 0
}

.merge_leak_logp <- function(pool_df, leak_prior_df, beta_leak = 0) {
  if (is.null(leak_prior_df) || beta_leak <= 0) {
    pool_df$logp_leak <- 0
    return(pool_df)
  }
  req <- c("digit_1","digit_2","digit_3","digit_4","logp")
  stopifnot(all(req %in% names(leak_prior_df)))
  merged <- merge(pool_df, leak_prior_df[, req],
                  by = c("digit_1","digit_2","digit_3","digit_4"),
                  all.x = TRUE, sort = FALSE)
  default_logp <- log(1/10000)
  merged$logp_leak[is.na(merged$logp)] <- default_logp
  merged$logp_leak[!is.na(merged$logp)] <- merged$logp[!is.na(merged$logp)]
  merged$logp_leak <- beta_leak * merged$logp_leak
  merged
}

.blockwise_logscore <- function(code_vec, priors, logphi_bundle,
                                df,
                                w12 = 1.0, w34 = 1.0, w_cross = 0.25,
                                socs_bonus_log = log(1.5)) {
  d1 <- code_vec[1]+1L; d2 <- code_vec[2]+1L; d3 <- code_vec[3]+1L; d4 <- code_vec[4]+1L
  L <- logphi_bundle$logphi
  s12_marg <- log(priors$p1[as.character(code_vec[1])]) + log(priors$p2[as.character(code_vec[2])])
  s34_marg <- log(priors$p3[as.character(code_vec[3])]) + log(priors$p4[as.character(code_vec[4])])
  s12_pair <- L[[1]][d1,d2]
  s34_pair <- L[[6]][d3,d4]
  cross_lp <- L[[2]][d1,d3] + L[[3]][d1,d4] + L[[4]][d2,d3] + L[[5]][d2,d4]
  socs_bump <- .socs_suffix_bonus(code_vec, df, bonus_log = socs_bonus_log)
  w12 * (s12_marg + s12_pair) +
    w34 * (s34_marg + s34_pair + socs_bump) +
    w_cross * cross_lp
}

score_all_codes <- function(df,
                            priors, logphi_bundle,
                            w12 = 1.0, w34 = 1.0, w_cross = 0.25,
                            socs_bonus_log = log(1.5),
                            leak_prior = NULL, beta_leak = 0) {
  pool <- .enumerate_all_codes()
  M <- as.matrix(pool[,1:4])
  pool$logscore <- apply(M, 1, .blockwise_logscore,
                         priors = priors,
                         logphi_bundle = logphi_bundle,
                         df = df,
                         w12 = w12, w34 = w34, w_cross = w_cross,
                         socs_bonus_log = socs_bonus_log)
  pool <- .merge_leak_logp(pool, leak_prior, beta_leak = beta_leak)
  pool$logscore <- pool$logscore + pool$logp_leak
  pool
}

apply_feedback <- function(cand_df, guesses = list(), ks = integer()) {
  if (length(guesses) == 0L) return(cand_df)
  if (length(guesses) != length(ks)) stop("guesses and ks must match lengths")
  codes_mat <- as.matrix(cand_df[, c("digit_1","digit_2","digit_3","digit_4")])
  keep <- rep(TRUE, nrow(cand_df))
  exclude <- rep(FALSE, nrow(cand_df))
  for (i in seq_along(guesses)) {
    g <- .standardize_guess(guesses[[i]])
    k <- as.integer(ks[[i]])
    if (length(g) != 4L || k < 0L || k > 4L) stop("invalid guess/k")
    m <- .matches_vec(codes_mat, g)
    keep <- keep & (m == k)
    ex <- (codes_mat[,1]==g[1]) & (codes_mat[,2]==g[2]) &
      (codes_mat[,3]==g[3]) & (codes_mat[,4]==g[4])
    exclude <- exclude | ex
  }
  out <- cand_df[keep & !exclude, , drop = FALSE]
  rownames(out) <- NULL
  out
}

.posterior_weights <- function(filt) {
  lw <- filt$logscore; lw <- lw - max(lw); exp(lw)
}
.posterior_marginals <- function(filt, w = NULL) {
  if (is.null(w)) w <- .posterior_weights(filt)
  M <- as.matrix(filt[, c("digit_1","digit_2","digit_3","digit_4")])
  out <- matrix(0, nrow = 4, ncol = 10)
  for (p in 1:4) {
    fac <- factor(M[,p], levels = 0:9)
    agg <- tapply(w, fac, sum)
    agg[is.na(agg)] <- 0
    out[p,] <- as.numeric(agg) / sum(agg)
  }
  out
}
detect_locks <- function(filt, threshold = 0.90) {
  w <- .posterior_weights(filt)
  marg <- .posterior_marginals(filt, w)
  locks <- rep(NA_integer_, 4)
  for (p in 1:4) {
    dstar <- which.max(marg[p,]) - 1L
    if (marg[p, dstar + 1L] >= threshold) locks[p] <- dstar
  }
  locks
}
restrict_pool_to_locks <- function(pool_df, locks) {
  if (all(is.na(locks))) return(pool_df)
  keep <- rep(TRUE, nrow(pool_df))
  for (p in 1:4) if (!is.na(locks[p])) {
    col <- paste0("digit_", p)
    keep <- keep & (pool_df[[col]] == locks[p])
  }
  out <- pool_df[keep, , drop = FALSE]
  if (!nrow(out)) return(pool_df)
  out
}

.build_cross_pool <- function(scored, K_prefix = 60, K_suffix = 60) {
  df <- scored[order(-scored$logscore), c("digit_1","digit_2","digit_3","digit_4","logscore")]
  agg12 <- aggregate(logscore ~ digit_1 + digit_2, data = df, FUN = max)
  agg34 <- aggregate(logscore ~ digit_3 + digit_4, data = df, FUN = max)
  top12 <- head(agg12[order(-agg12$logscore), c("digit_1","digit_2")], K_prefix)
  top34 <- head(agg34[order(-agg34$logscore), c("digit_3","digit_4")], K_suffix)
  cross <- merge(top12, top34, by = NULL)
  cross[, c("digit_1","digit_2","digit_3","digit_4")]
}

.weighted_bucket <- function(m, w, nbins = 5L) {
  fac <- factor(m, levels = 0:(nbins-1))
  agg <- tapply(w, fac, sum)
  out <- ifelse(is.na(agg), 0, agg)
  as.numeric(out)
}
best_split_guess <- function(cand_df,
                             probe_pool = NULL,
                             sample_limit = 2500,
                             use_mass_split = TRUE,
                             seed = 1L) {
  stopifnot(nrow(cand_df) > 0L)
  logw <- cand_df$logscore; logw <- logw - max(logw)
  w <- exp(logw); Wtot <- sum(w)
  codes_mat <- as.matrix(cand_df[,1:4])
  if (is.null(probe_pool)) probe_pool <- .enumerate_all_codes()[,1:4]
  if (!is.na(sample_limit) && nrow(probe_pool) > sample_limit) {
    set.seed(seed); probe_pool <- probe_pool[sample(seq_len(nrow(probe_pool)), sample_limit), , drop = FALSE]
  }
  pool_mat <- as.matrix(probe_pool)
  best_cost <- Inf; best_g <- NULL; best_tie <- -Inf
  for (i in seq_len(nrow(pool_mat))) {
    g <- as.integer(pool_mat[i, ])
    m <- .matches_vec(codes_mat, g)
    if (use_mass_split) {
      Wk <- .weighted_bucket(m, w, nbins = 5L)
      pk <- Wk / max(Wtot, 1e-12)
      cost <- sum(Wk * Wk) / max(Wtot, 1e-12)
      tie  <- -sum(pk * log(pmax(pk, 1e-12)))
    } else {
      sizes <- tabulate(m + 1L, nbins = 5L)
      pk <- sizes / nrow(codes_mat)
      cost <- sum(sizes * sizes) / nrow(codes_mat)
      tie  <- -sum(pk * log(pmax(pk, 1e-12)))
    }
    if (cost < best_cost - 1e-12 || (abs(cost - best_cost) <= 1e-12 && tie > best_tie)) {
      best_cost <- cost; best_g <- g; best_tie <- tie
    }
  }
  list(guess = best_g, expected_cost = best_cost)
}

generate_isolators <- function(filt, base_guess = NULL, lock_threshold = 0.90) {
  w <- .posterior_weights(filt)
  marg <- .posterior_marginals(filt, w)
  if (is.null(base_guess)) {
    base_guess <- as.integer(filt[1, 1:4])
  } else {
    base_guess <- .standardize_guess(base_guess)
  }
  locks <- detect_locks(filt, threshold = lock_threshold)
  uncert <- which(is.na(locks))
  if (length(uncert) == 0L) return(list(testA = base_guess, testB = base_guess))
  least_digit <- function(p, avoid) {
    ord <- order(marg[p, ], decreasing = FALSE)
    d <- ord[1]-1L
    if (d == avoid) d <- ord[2]-1L
    d
  }
  repl <- sapply(1:4, function(p) least_digit(p, base_guess[p]))
  u_scores <- sapply(uncert, function(p) 1 - max(marg[p, ]))
  ord_u <- uncert[order(u_scores, decreasing = TRUE)]
  half <- ceiling(length(ord_u)/2)
  A_flip <- ord_u[seq_len(half)]
  B_flip <- setdiff(ord_u, A_flip)
  testA <- base_guess; if (length(A_flip)) for (p in A_flip) testA[p] <- repl[p]
  testB <- base_guess; if (length(B_flip)) for (p in B_flip) testB[p] <- repl[p]
  list(testA = testA, testB = testB)
}

propose_next_guess <- function(df,
                               guesses = list(), ks = integer(),
                               laplace = 1, alpha_pmi = 0.6,
                               w12 = 1.0, w34 = 1.0, w_cross = 0.25,
                               socs_bonus_log = log(1.5),
                               leak_prior = NULL, beta_leak = 0,
                               use_cross_pool = TRUE, K_prefix = 60, K_suffix = 60,
                               use_mass_split = TRUE,
                               sample_limit = 2500,
                               seed = 1L,
                               top_n = 20,
                               respect_locks = TRUE, lock_threshold = 0.90,
                               always_offer_isolators = TRUE) {
  pri <- build_priors_nb(df, laplace = laplace)
  logphi <- build_pairwise_logphi(df, alpha = alpha_pmi, laplace = laplace)
  scored <- score_all_codes(df, priors = pri, logphi_bundle = logphi,
                            w12 = w12, w34 = w34, w_cross = w_cross,
                            socs_bonus_log = socs_bonus_log,
                            leak_prior = leak_prior, beta_leak = beta_leak)
  filt <- apply_feedback(scored, guesses = guesses, ks = ks)
  if (!nrow(filt)) stop("No candidates after filtering. Check guesses/k.")
  filt <- filt[order(-filt$logscore), ]; rownames(filt) <- NULL
  locks <- if (respect_locks) detect_locks(filt, threshold = lock_threshold) else rep(NA_integer_, 4)
  probe_pool <- NULL
  if (use_cross_pool) {
    cross_pool <- .build_cross_pool(filt, K_prefix = K_prefix, K_suffix = K_suffix)
    if (nrow(cross_pool)) probe_pool <- cross_pool
  }
  if (!is.null(probe_pool)) {
    probe_pool <- restrict_pool_to_locks(probe_pool, locks)
  } else if (any(!is.na(locks))) {
    probe_pool <- restrict_pool_to_locks(.enumerate_all_codes(), locks)
  }
  bs <- best_split_guess(cand_df = filt,
                         probe_pool = probe_pool,
                         sample_limit = sample_limit,
                         use_mass_split = use_mass_split,
                         seed = seed)
  isolators <- if (always_offer_isolators) generate_isolators(filt, base_guess = bs$guess,
                                                              lock_threshold = lock_threshold) else NULL
  list(
    remaining     = nrow(filt),
    next_guess    = as.integer(bs$guess),
    expected_cost = bs$expected_cost,
    top           = head(filt[, c("digit_1","digit_2","digit_3","digit_4","logscore")], top_n),
    locks         = locks,
    isolator_tests= isolators
  )
}
