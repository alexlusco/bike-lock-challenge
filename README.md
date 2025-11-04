# bike lock challenge

**authors:** alex luscombe & m.p.

**disclaimer:** for education purposes only! this algorithm was created for a mastermind-style guessing game (5 guesses max), and *not* for the purpose of stealing bikes. don't do that.

## overview

this repo holds a small model for guessing a hidden 4-digit code (0000–9999) using only **black-peg feedback**:

- you play a 4-digit guess  
- you are told `k` = how many digits are correct **in the correct position**  
- you get no information about “right digit, wrong place”

the goal is to get strong partial hits (k ≥ 2, ideally k ≥ 3) in as few guesses as possible, using a small “jumble” dataset of past combos plus an optional leak prior.

---

## high-level idea

the model builds a probabilistic picture of 4-digit codes from:

1. your **jumble data** (`real_data` / `exp_data`): a few dozen historical 4-digit combos, treated as samples from your personal distribution.
2. an optional **leak prior** (`pin_data`): a ranked list of popular 4-digit codes (e.g., from a larger leak).

it then:

- scores all 10,000 possible codes using that distribution;
- filters codes to respect all previous `k` values; and
- proposes guesses that either:
  - **shrink** the remaining space as efficiently as possible, or  
  - jump to the current **map** (maximum a-posteriori) candidate.

---

## model assumptions

plain english assumptions baked into the scoring:

- **digits are not uniform.** some digits show up more often at each wheel (e.g., leading “1” or “0”). the model learns these per-position frequencies from your jumbles.
- **positions are weakly coupled.** certain pairs like (1,2) or (3,4) occur together more often than chance (e.g., 19xx, xx00). we capture this with light pairwise dependencies.
- **suffixes may “share offsets.”** people often reuse simple patterns at the end (e.g., shifting a favorite number a few steps). the model adds a small bonus when the last two digits look like a simple back-off from the dominant suffix in your data.
- **feedback is exact on k only.** the hidden code is assumed to be one of the 10,000 possibilities, and any candidate that does not match *all* past `k` values is thrown out.
- **a leak prior is soft information.** if you enable it, global popularity nudges the scores but does not override what your jumbles and feedback say.

none of these are treated as hard rules; they just tilt probabilities.

---

## how the scoring works

for each 4-digit code, the model computes a log-score as:

1. **per-wheel priors (naive bayes).**  
   from your jumbles, estimate `p(digit_d at position_p)` with laplace smoothing. this captures marginal preferences like “0 is common in the last slot.”

2. **pairwise coupling (pmi).**  
   for all position pairs ((1,2), (1,3), (1,4), (2,3), (2,4), (3,4)) we compute **pointwise mutual information (pmi)** from your jumbles:
   - pmi > 0 means two digits co-occur together more often than expected under independence;
   - pmi < 0 means they avoid each other.

   we then add scaled pmi terms to the score with weight `alpha_pmi`, and emphasize (1,2) and (3,4) via `w12` and `w34`.

3. **shared-offset suffix bump.**  
   using the most common digits on wheels 3 and 4, we give a small bonus (`socs_bonus_log`) to candidates where the last two digits look like a simple offset from those modes (within a small range). this is a heuristic to capture “people often shift their favorite suffix by a few steps.”

4. **optional leak prior.**  
   if you build `leak_prior_table` from `pin_data` and set:
   ```r
   params$leak_prior <- leak_prior_table
   params$beta_leak  <- 0.10  # for example
   ```
   then each code gets an extra term `beta_leak * logp_leak(code)`. small `beta_leak` = gentle nudge; large `beta_leak` = trust leak more.

5. **total log-score.**  
   all components are summed into a single `logscore` over the full 10k code space. higher `logscore` means more plausible under the current model.

---

## using feedback: pruning the space

after each guess, you enter `k` (0–4). the algorithm then:

1. **filters candidates.**  
   for every past guess `g_i` with feedback `k_i`, we only keep codes that have **exactly** `k_i` matching positions with `g_i`. anything inconsistent is removed.

2. **re-weights candidates.**  
   among the survivors, `logscore` is treated as a posterior log-probability. codes are sorted from most to least likely.

3. **detects locks.**  
   for each position, it accumulates posterior mass over digits 0–9. if one digit holds more than `lock_threshold` (e.g., 0.88) of the mass at that position, the position is marked as “locked” and future probes prefer to keep it fixed.

---

## shrink vs map vs isolators

each round, `propose_next_guess()` returns:

- `next_guess`: the **shrink guess**. this is the probe that best *splits* the remaining posterior mass across the five possible `k` buckets (0–4). it is chosen by:
  - optionally restricting to a **cross-pool** of strong (prefix, suffix) combinations;  
  - optionally sampling a subset (`sample_limit`) of that pool;  
  - simulating the distribution of `k` for each candidate and picking the one that minimizes expected “mass concentration” (i.e., keeps the posterior as spread out as possible).

- `top`: a table of the top `top_n` candidates by log-score. the first row is the **map guess** (most probable code under the current model).

- `isolator_tests`: two “isolator” guesses `testA` and `testB`:
  - they start from a strong base guess (usually the shrink guess);
  - they deliberately **flip digits** in the most uncertain positions to probe which slots are correct.

your script prints all these options each round:

- `shrink` = best information-gain probe  
- `map`    = highest-probability code  
- `isol a` / `isol b` = diagnostics for which digit slots are right

you then manually pick which one to play and enter the observed `k`.

---

## main knobs and how to tune them

all tuning is through the `params` list in `01_build_guesses.R`.

### structure / prior

- `laplace`  
  smoothing for per-wheel frequencies.  
  - larger ⇒ more uniform digits, less extreme priors.  
  - smaller ⇒ priors follow your jumbles more closely.

- `alpha_pmi`  
  how much to trust pairwise dependencies.  
  - 0.0 ⇒ treat positions as independent.  
  - 0.7–0.9 ⇒ strong coupling; use if your jumbles clearly prefer certain pairs (e.g., 19, 00, 23).

- `w12`, `w34`  
  weights on the (1,2) and (3,4) blocks.  
  - >1.0 ⇒ emphasize evidence from that pair;  
  - <1.0 ⇒ down-weight it.

- `w_cross`  
  weight on cross-pairs (1,3), (1,4), (2,3), (2,4).  
  - small (≈ 0.2) ⇒ treat them as weak hints;  
  - larger ⇒ cross-pair structure matters more.

- `socs_bonus_log`  
  log-bonus for suffixes that look like simple offsets around the suffix modes.  
  - 0 ⇒ no suffix heuristic;  
  - `log(1.3)–log(1.7)` ⇒ gentle bump.

### leak prior

- `leak_prior`  
  `NULL` = leak off.  
  `leak_prior_table` = use `build_leak_prior()` over a ranked `pins` list.

- `beta_leak`  
  mixing weight for the leak log-probabilities.  
  - start at 0 (off);  
  - increase to ~0.10–0.15 if, after a couple of rounds, the shrink guesses are not moving `k` or the remaining space.

### shrink engine / search

- `use_cross_pool`, `K_prefix`, `K_suffix`  
  if `TRUE`, build a focused candidate pool from the best `K_prefix` (digit_1, digit_2) pairs and `K_suffix` (digit_3, digit_4) pairs.

- `use_mass_split`  
  `TRUE` ⇒ pick shrink guesses that split posterior **mass** evenly.  
  `FALSE` ⇒ split by raw candidate **counts** only.

- `sample_limit`  
  cap on how many probes to test when searching for the best shrink guess.

### locks and isolators

- `respect_locks`, `lock_threshold`  
  if `respect_locks = TRUE`, positions whose posterior mass exceeds `lock_threshold` are treated as fixed in future shrink probes.

- `always_offer_isolators`  
  if `TRUE`, compute `isol a` and `isol b` every round.

---

## practical recipe

a simple way to use and tune this:

1. **start with mp-style defaults**, e.g.
   ```r
   laplace        = 1
   alpha_pmi      = 0.7
   w12            = 1.1
   w34            = 1.2
   w_cross        = 0.22
   socs_bonus_log = log(1.5)
   leak_prior     = NULL
   beta_leak      = 0
   ```

2. **play shrink first.** use the shrink guess for round 1, and usually round 2, until either:
   - `k ≥ 2`, or  
   - remaining candidates ≤ ~3000.

3. **then switch to map or an isolator.** once the posterior is tight, using `map` or `isol a / b` is often the fastest way to reach k ≥ 3 and then k = 4.

4. **only enable leak if stuck.** if after two shrink guesses you are still at k ≤ 1 and the candidate space is large, set:
   ```r
   params$leak_prior <- leak_prior_table
   params$beta_leak  <- 0.1
   ```
   and recompute from the current `guesses` and `ks`.

5. **change one knob at a time.** watch how it affects:
   - how fast you reach k ≥ 2;  
   - what the top `map` candidates look like after each round.

this keeps the procedure simple: let the jumbles define the structure, use shrink to carve down the space, then lean on map and isolators once the model is confident.
