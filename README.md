# bike-lock-challenge

## overview

this repo contains a small R project built to model how a person’s bike lock habits might leak information about their true 4-digit code.  

the core script builds probabilistic priors for each wheel of a lock based on observed "jumble" data — the random numbers visible after someone scrambles the lock. it uses entropy to estimate how much each wheel tends to spin and applies a circular kernel to weight digits near frequently stopped positions.

from these priors, the script scores all 10,000 possible combinations and can optionally ban unlikely starting digits (e.g., 7–9).  

this is built for a small game where we’re allowed up to five guesses. after each guess, we’re told how many digits were in the correct position (the number of `k` spots). 

each time we get that feedback, the model filters the candidate pool to keep only codes consistent with the result, renormalizes the probabilities, and re-ranks what’s left. multiple guesses can be chained together to iteratively shrink the search space.  
the output is a ranked list of the most likely remaining combinations.

---

developed by alex luscombe & m.p.

october 2025
