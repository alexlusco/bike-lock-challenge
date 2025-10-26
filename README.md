# bike-lock-challenge

## overview
this repo models how a person’s bike-lock habits leak information about their true 4-digit code. from observed “jumble” stops (the digits visible after a scramble), we build wheel-wise priors and combine them with behavioral priors from real-world pin use and simple human patterns. the result is a posterior over all 10,000 codes that updates after each guess + feedback (`k` correct in position, or “anywhere” if you choose that mode). a two-signal policy then proposes the next guess.

developed by alex luscombe & m.p. — october 2025

---

## how it works

### 1) wheel priors from jumbles (kernel model)
- estimate **spin** per wheel by the entropy of observed stops.
- map spin → smoothing **λ** (sharper when spin is low, flatter when spin is high).
- apply a symmetric **circular kernel** on digits 0–9 for each wheel → `p1..p4`.

### 2) behavioral priors (can toggle off)
- **full-code pin prior (csv):** read a ranked list of 4-digit pins (col `pins`) and convert rank to probability via `r^{-alpha}`. this gives a popularity prior over all 10k codes.
- **pattern prior (bike-lock tuned):** simple features (aaaa, aabb, abab, abba, increasing/decreasing sequences with wrap, mmdd/ddmm dates, 19xx/20xx years). convert `w·x` to a probability with a softmax.
- **blend:** multiply kernel posterior by `(pin_prior^eta_pin)` and `(pattern_prior^eta_patt)`, then renormalize.

### 3) calibration (rolling cv by week)
- using week information (in this version of the game, 20 new rows are released/week), the code runs **rolling cv** to choose `(ent_min, ent_max, lambda_sharp, lambda_flat)` that maximize out-of-sample log-lik on wheel digits.
- In absence of `week` column, falls back to defaults.

### 4) ranking + filtering after feedback
- compute weights for all 10,000 codes; optionally **ban first digits** (e.g., 7–9).
- after each guess, filter survivors by the feedback rule:
  - `match_mode="pos"` → exact slot matches (bike-lock default).
  - `match_mode="any"` → mastermind-style anywhere matches.
- renormalize weights and re-rank.

### 5) next-guess policy (two signals)
- **signal a (map):** the most probable surviving code.
- **signal b (shrink):** the code that maximizes expected information gain (or expected remaining mass) over a candidate pool (top-n/top-k/all).
- **decision:** if `Pr(MAP) >= tau_map`, policy recommends to shoot map; otherwise take the shrink pick.

### 6) bootstrap (work in progress)
- use bootstrapping to produce a **consensus next guess** (potentially handy early game when priors dominate?).

---
