# Baseline results (k = 4)

Searching for a 2-colouring of K_n with no one-colour K_4, for n = 5–17, with
50 seeds (0–49) and a budget of 1,000,000 steps per run. Hitting the budget
counts as a failure. The LLL guarantees success up to n = 5, and valid
colourings exist up to n = 17 because R(4,4) = 18.

## Main findings

- **Greedy solved every run at every n**, including n = 17 in a median of 1,404
  flips and 0.61 ms. MT solved 0/50 at that size.
- **Annealing reached n = 15 in every run**, then 47/50 at n = 16 and 14/50 at
  n = 17, using the untuned starting values.
- **Restart solved every run up to n = 12**, then 1/50 at n = 13.
- **MT got only one n further than restart:** it solved every run up to n = 13,
  then 1/50 at n = 14. It was far cheaper at n = 12, with 1,172 steps against
  65,551 and 0.22 ms against 51.09 ms.
- **MT's failures from n = 14 are not because solutions are out of reach.** Local
  search finds them at those sizes in milliseconds.

## Runs solved (out of 50)

| n | MT | Restart | Greedy | Annealing |
|---|---:|---:|---:|---:|
| 5 | 50 | 50 | 50 | 50 |
| 6 | 50 | 50 | 50 | 50 |
| 7 | 50 | 50 | 50 | 50 |
| 8 | 50 | 50 | 50 | 50 |
| 9 | 50 | 50 | 50 | 50 |
| 10 | 50 | 50 | 50 | 50 |
| 11 | 50 | 50 | 50 | 50 |
| 12 | 50 | 50 | 50 | 50 |
| 13 | 50 | 1 | 50 | 50 |
| 14 | 1 | 0 | 50 | 50 |
| 15 | 0 | 0 | 50 | 50 |
| 16 | 0 | 0 | 50 | 47 |
| 17 | 0 | 0 | 50 | 14 |

## Median steps

Taken over all 50 runs. "—" means at most half were solved, so the median would
just be the budget.

| n | MT | Restart | Greedy | Annealing |
|---|---:|---:|---:|---:|
| 5 | 0 | 1 | 0 | 0 |
| 6 | 0 | 1 | 0 | 0 |
| 7 | 1 | 2 | 1 | 2 |
| 8 | 2 | 4 | 1 | 6 |
| 9 | 5 | 16 | 2 | 24 |
| 10 | 14 | 72 | 3 | 66 |
| 11 | 62 | 654 | 5 | 341 |
| 12 | 1,172 | 65,551 | 8 | 2,800 |
| 13 | 145,200 | — | 13 | 46,788 |
| 14 | — | — | 22 | 117,560 |
| 15 | — | — | 76 | 163,137 |
| 16 | — | — | 396 | 197,410 |
| 17 | — | — | 1,404 | — |

## Median time per run (ms)

Same rule as above. Sizes below n = 10 all took under 0.01 ms.

| n | MT | Restart | Greedy | Annealing |
|---|---:|---:|---:|---:|
| 10 | <0.01 | 0.03 | <0.01 | 0.01 |
| 11 | 0.01 | 0.33 | <0.01 | 0.03 |
| 12 | 0.22 | 51.09 | <0.01 | 0.25 |
| 13 | 33.20 | — | 0.01 | 4.40 |
| 14 | — | — | 0.01 | 11.47 |
| 15 | — | — | 0.03 | 16.16 |
| 16 | — | — | 0.16 | 20.30 |
| 17 | — | — | 0.61 | — |

## Fewest bad sets on failed runs (min / median / max)

Greedy had no failed runs, and no method failed below n = 13.

| n | MT | Restart | Annealing |
|---|---:|---:|---:|
| 13 | — | 1 / 2 / 2 | — |
| 14 | 1 / 2 / 3 | 2 / 5 / 5 | — |
| 15 | 5 / 7 / 8 | 8 / 10 / 11 | — |
| 16 | 11 / 13 / 14 | 14 / 16 / 17 | 1 / 2 / 2 |
| 17 | 19 / 21 / 23 | 23 / 24 / 26 | 4 / 4 / 5 |

## Methods

| Method | One step | Settings | Total sweep time |
|---|---|---|---:|
| MT | resample one bad set (6 edges) | `random` rule | 68.6 s |
| Restart | one fresh colouring, recounting every set | none | 606.3 s |
| Greedy | one edge flip | noise 0.2 | 0.1 s |
| Annealing | one proposed flip | t0 2.0, cooling 0.99999, t_min 0.05 | 6.4 s |

## Caveats

- **Steps are not equal work across methods.** Compare time as well as steps.
  Restart's total time is summed across 4 processes running at once, while the
  other sweeps ran one at a time, so its timings are approximate.
- **A failure only means the budget ran out.** It doesn't show a method can never
  succeed, and the 10× budget check covers this for MT.
- **Annealing could do better with tuning.** The temperature reaches its floor
  after about 369,000 steps. Every win came by step 268,683, and no failed run
  improved after the floor, so about 63% of each failed run's budget did nothing.
- **This section covers k = 4 only**, with one greedy noise setting.

## Checks

- Replayed all 2,113 solved runs from the CSVs: every colouring was valid and
  every step count reproduced.
- My MT sweep matched the team's `k4_T1e6.csv` on all 650 rows for `success`,
  `steps`, `min_bad_seen` and `step_of_min`.
- At n = 17, all 64 greedy and annealing solutions (50 + 14) passed a separate
  brute-force recount and were 8-regular, as the Paley graph is.
- At n = 18, where no valid colouring exists, greedy and annealing each solved
  0/20 runs, and the fewest bad sets reached was 9.

## Reproduce

CSVs are in `results/`. Results are appended, so delete the old file first.

    python3 run_experiment.py --method mt --k 4 --n 5 17 --seeds 50 --budget 1000000 --out results/k4_mt.csv
    python3 run_experiment.py --method restart --k 4 --n 5 17 --seeds 50 --budget 1000000 --out results/k4_restart.csv
    python3 run_experiment.py --method greedy --k 4 --n 5 17 --seeds 50 --budget 1000000 --out results/k4_greedy.csv
    python3 run_experiment.py --method annealing --k 4 --n 5 17 --seeds 50 --budget 1000000 --out results/k4_annealing.csv

# Baseline results (k = 5)

Searching for a 2-colouring of K_n with no one-colour K_5, with 50 seeds
(0–49) and a budget of 1,000,000 steps per run. Greedy and annealing were
scouted at 10 seeds over n = 20–42; restart was scouted to n = 24 because
it fails from 22. The 50-seed runs cover n = 20–24 and the frontier around
each method's first failures.

## Main findings

- **Greedy is the strongest by a wide margin.** It solved every 50-seed run
  at n = 20–24 and n = 35–36, and the 10-seed scout was clean through
  n = 36. It then dropped to 16/50 at n = 37 and 8/50 at n = 38; the scout
  found one success at n = 38 and none at n = 39–42.
- **Annealing reaches n = 35.** It solved every 50-seed run at n = 20–24
  and n = 34–35, and the 10-seed scout was clean through n = 35. It then
  dropped to 5/50 at n = 36 and 0/50 at n = 37.
- **Restart collapses at n = 22.** It solved every seed at n = 20, 19/50
  at n = 21, 1/50 at n = 22, and 0/50 at n = 23–24.
- The frontier is much higher than at k = 4: greedy n ≈ 37, annealing
  n ≈ 35, restart n ≈ 21.

## Runs solved (out of 50)

| n | Restart | Greedy | Annealing |
|---|---:|---:|---:|
| 20 | 50 | 50 | 50 |
| 21 | 19 | 50 | 50 |
| 22 | 1 | 50 | 50 |
| 23 | 0 | 50 | 50 |
| 24 | 0 | 50 | 50 |
| 34 | — | — | 50 |
| 35 | — | 50 | 50 |
| 36 | — | 50 | 5 |
| 37 | — | 16 | 0 |
| 38 | — | 8 | — |

The frontier rows are not a full sweep: greedy was run at n = 35–38 and
annealing at n = 34–37. A dash means the size was not tested at 50 seeds.

## 10-seed scout

Greedy solved 10/10 at every n from 20 to 36, then 3/10 at n = 37, 1/10
at n = 38, and 0/10 at n = 39–42. Annealing solved 10/10 through n = 35
and 0/10 from n = 36 to 42. Restart solved 10/10 at n = 20, 4/10 at
n = 21, 1/10 at n = 22, and 0/10 at n = 23–24.

## Median steps

Taken over solved runs only; a dash means no run solved.

| n | Restart | Greedy | Annealing |
|---|---:|---:|---:|
| 20 | 49,982 | 8 | 932 |
| 21 | 489,503 | 11 | 3,082 |
| 22 | 580,966 | 12 | 10,059 |
| 23 | — | 16 | 37,230 |
| 24 | — | 20 | 65,208 |
| 34 | — | — | 207,126 |
| 35 | — | 5,026 | 247,900 |
| 36 | — | 55,964 | 755,897 |
| 37 | — | 475,908 | — |
| 38 | — | 471,934 | — |

At n = 37 and 38 greedy's median is over the 16 and 8 successful runs;
the rest hit the budget.

## Timing

Median time per run: restart 2.2 s at n = 20, 61 s at n = 21, 81 s at
n = 22, 204 s at n = 23 and 252 s at n = 24. Greedy takes 6 s at n = 36
and 64–72 s at n = 37–38. Annealing takes 3.4 s at n = 35 and 12–13 s at
n = 36–37. Restart's failures are the expensive ones; each spends the
full budget.

## Caveats

- **Steps are not equal work across methods.** Greedy's steps are cheap
  edge flips; restart's steps are whole fresh colourings.
- **Annealing uses the untuned k = 4 settings.** t0 = 2.0, cooling =
  0.99999, t_min = 0.05.
- **The scout is only 10 seeds.** The 50-seed runs bracket the frontier;
  the scout gives the shape from n = 20 to 42.
- **Restart n = 23–24 are all failures.** They confirm the frontier is
  below 23, at the cost of a full budget per seed.

## Reproduce

CSVs are in `results/`. Results are appended, so delete the old file first.

    python3 run_experiment.py --method restart --k 5 --n 20 24 --seeds 50 --budget 1000000 --out results/k5_restart.csv
    python3 run_experiment.py --method greedy --k 5 --n 20 24 --seeds 50 --budget 1000000 --out results/k5_greedy.csv
    python3 run_experiment.py --method annealing --k 5 --n 20 24 --seeds 50 --budget 1000000 --out results/k5_annealing.csv
    python3 run_experiment.py --method greedy --k 5 --n 35 38 --seeds 50 --budget 1000000 --out results/k5_greedy_frontier.csv
    python3 run_experiment.py --method annealing --k 5 --n 34 37 --seeds 50 --budget 1000000 --out results/k5_annealing_frontier.csv
    python3 run_experiment.py --method greedy --k 5 --n 20 42 --seeds 10 --budget 1000000 --out results/k5_greedy_scout.csv
    python3 run_experiment.py --method annealing --k 5 --n 20 42 --seeds 10 --budget 1000000 --out results/k5_annealing_scout.csv
    python3 run_experiment.py --method restart --k 5 --n 20 24 --seeds 10 --budget 1000000 --out results/k5_restart_scout.csv
