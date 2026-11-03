# Moser–Tardos for Ramsey colourings (CS4234 mini-project)

Finds 2-colourings of the edges of K_n with no monochromatic K_k using the
Moser–Tardos resampling algorithm. A success certifies R(k,k) > n.

## Setup
    pip install -r requirements.txt

## Use
    from ramsey_mt.core import solve, RamseyInstance, n_lll
    r = solve(n=12, k=4, budget=10**6, seed=0, rule="random")
    # r: success, steps, min_bad_seen, step_of_min, time_ms, final_bad, colours

Selection rules: "random", "last" (most recently violated), "min_index".
Reuse `RamseyInstance(n, k)` across seeds to avoid rebuilding the event structure.

## Experiments
    python3 run_experiment.py --k 4 --n 5 17 --seeds 50 --budget 1000000 --out results/k4.csv

CSV schema (shared with the baselines):
    method,k,n,seed,budget,success,steps,min_bad_seen,step_of_min,time_ms

## Tests
    python3 -m pytest -q tests
Includes a check that the Paley colouring of K_17 has no monochromatic K_4.

## Interface for baselines
Baselines should expose the same signature and return the same keys:
    solve(n, k, budget, seed) -> {success, steps, min_bad_seen, step_of_min, time_ms}
`RamseyInstance.count_monochromatic(colours)` can be used to score any colouring.
