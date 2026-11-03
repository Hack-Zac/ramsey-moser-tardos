"""Sweep Moser–Tardos over n for a fixed k and append results to a CSV.

Example:
  python3 run_experiment.py --k 4 --n 5 17 --seeds 50 --budget 1000000 --rule random --out results/k4.csv

CSV schema (shared with the baselines):
  method,k,n,seed,budget,success,steps,min_bad_seen,step_of_min,time_ms
"""
import argparse
import csv
import os

from ramsey_mt.core import RamseyInstance, solve

FIELDS = ["method", "k", "n", "seed", "budget", "success", "steps",
          "min_bad_seen", "step_of_min", "time_ms"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--k", type=int, required=True)
    ap.add_argument("--n", type=int, nargs=2, metavar=("N_MIN", "N_MAX"), required=True)
    ap.add_argument("--seeds", type=int, default=20, help="seeds 0..SEEDS-1")
    ap.add_argument("--budget", type=int, default=10**6, help="max resamplings T")
    ap.add_argument("--rule", default="random", choices=["random", "last", "min_index"])
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

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
                r = solve(n, args.k, args.budget, seed, args.rule, instance=inst)
                wins += r["success"]
                w.writerow({"method": f"mt-{args.rule}", "k": args.k, "n": n,
                            "seed": seed, "budget": args.budget,
                            **{key: r[key] for key in FIELDS[5:]}})
            f.flush()
            print(f"k={args.k} n={n:3d}: {wins}/{args.seeds} solved", flush=True)


if __name__ == "__main__":
    main()
