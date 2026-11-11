"""Sweep a search method over n for a fixed k and append results to a CSV.

Example:
  python3 run_experiment.py --method mt --k 4 --n 5 17 --seeds 50 --budget 1000000 --out results/k4_mt.csv

CSV schema (shared with the baselines):
  method,k,n,seed,budget,success,steps,min_bad_seen,step_of_min,time_ms
"""
import argparse
import csv
import os
from functools import partial

import baselines
from core import RamseyInstance, solve

FIELDS = ["method", "k", "n", "seed", "budget", "success", "steps",
          "min_bad_seen", "step_of_min", "time_ms"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--k", type=int, required=True)
    ap.add_argument("--n", type=int, nargs=2, metavar=("N_MIN", "N_MAX"), required=True)
    ap.add_argument("--seeds", type=int, default=20, help="seeds 0..SEEDS-1")
    ap.add_argument("--budget", type=int, default=10**6, help="max steps per run")
    ap.add_argument("--method", default="mt", choices=["mt", "restart", "greedy", "annealing"])
    ap.add_argument("--rule", default="random", choices=["random", "last", "min_index"],
                    help="selection rule (mt only)")
    ap.add_argument("--noise", type=float, default=0.2, help="random-flip probability (greedy only)")
    ap.add_argument("--t0", type=float, default=2.0, help="start temperature (annealing only)")
    ap.add_argument("--cooling", type=float, default=0.99999, help="per-step factor (annealing only)")
    ap.add_argument("--t-min", type=float, default=0.05, help="min temperature (annealing only)")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    if args.method == "mt":
        run, label = partial(solve, rule=args.rule), f"mt-{args.rule}"
    elif args.method == "restart":
        run, label = baselines.random_restart, "restart"
    elif args.method == "greedy":
        run, label = partial(baselines.greedy_local_search, noise=args.noise), f"greedy-p{args.noise}"
    else:
        run = partial(baselines.simulated_annealing, t0=args.t0,
                      cooling=args.cooling, t_min=args.t_min)
        label = f"annealing-t{args.t0}-c{args.cooling}"

    # Compile now so Numba doesn't inflate the first timed run.
    run(args.k + 1, args.k, 10, 0)

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    new_file = not os.path.exists(args.out)
    with open(args.out, "a", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        if new_file:
            w.writeheader()
        for n in range(args.n[0], args.n[1] + 1):
            inst = RamseyInstance(n, args.k)
            wins = 0
            for seed in range(args.seeds):
                r = run(n, args.k, args.budget, seed, instance=inst)
                wins += r["success"]
                w.writerow({"method": label, "k": args.k, "n": n,
                            "seed": seed, "budget": args.budget,
                            **{key: r[key] for key in FIELDS[5:]}})
            f.flush()
            print(f"k={args.k} n={n:3d}: {wins}/{args.seeds} solved", flush=True)


if __name__ == "__main__":
    main()
