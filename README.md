# bike-lock challenge — 4-digit code cracking from jumbles

**authors:** alex luscombe & m.p.

**disclaimer:** for education purposes only! this algorithm was created for a mastermind-style guessing game (5 guesses max), and *not* for the purpose of stealing bikes. don't do that.

## overview
a small, reproducible framework to guess a 4-digit combination lock code using only:
- a few dozen rows of jumble data (the digits you tend to land on when you spin each wheel), and
- the black-peg feedback k you get after each guess (k = number of digits correct *and* in the correct position).

---

## approach (technical)
1. **per-wheel priors (naive bayes):** estimate p(digit at each position) from jumbles with laplace smoothing.
2. **pairwise coupling (pmi):** learn light dependencies between positions with pointwise mutual information (pmi), especially (1,2) and (3,4).
3. **shared-offset suffix bump (optional):** a small log-prob bonus if the suffix (digits 3 and 4) matches a simple **shared back-off** pattern seen in jumbles. 
4. **optional leak prior:** add a soft log-probability over the full 10k code space (weighted by a mixing knob).
5. **candidate scoring:** sum all components to get a log-score per code on the full 10k space.
6. **feedback filtering:** enforce every past guess’s `k` (keep only codes with exactly that many positional matches).
7. **map vs shrink:**  
   - **map** = top posterior code (highest current log-score).  
   - **shrink** = next probe that **maximizes expected reduction** of the remaining mass (best split on k=0..4 buckets), sampled from a focused pool.
8. **isolator tests:** build 2 orthogonal tests that flip digits at the most uncertain positions to identify which slots are correct.
9. **locks:** when a position’s posterior mass ≥ threshold (e.g., 0.9), treat it as “locked” and focus probes around remaining positions.

---

