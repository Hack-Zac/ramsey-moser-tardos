"""k=5 comparison: success rate and median time vs n. Merges 50-seed runs (main+frontier) over 10-seed scouts.
Usage: python compare_k5.py DATA_DIR OUT_PREFIX"""
import csv, sys, os, collections
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#9a9893"
STYLE = {"mt": ("Moser–Tardos", "#2a78d6", "o"), "restart": ("Random restart", "#eb6834", "s"),
         "annealing": ("Simulated annealing", "#1baf7a", "^"), "greedy": ("Greedy local search", "#eda100", "D")}
N_LLL, N_BEST = 9, 42
D, out = sys.argv[1], sys.argv[2]
def rows(f):
    p = os.path.join(D, f); return list(csv.DictReader(open(p))) if os.path.exists(p) else []
data = {}
for m in STYLE:
    by = collections.defaultdict(dict)
    for r in rows(f"k5_{m}_scout.csv"): by[int(r["n"])][r["seed"]] = r
    full = rows(f"k5_{m}.csv") + rows(f"k5_{m}_frontier.csv")
    for n in {int(r["n"]) for r in full}: by[n] = {}
    for r in full: by[int(r["n"])][r["seed"]] = r
    ns = sorted(by); rate, t, cnt = [], [], []
    for n in ns:
        R = list(by[n].values()); ok = np.array([r["success"] == "True" for r in R])
        rate.append(ok.mean()); cnt.append(len(R))
        t.append(np.median([float(r["time_ms"]) for r in R]) if ok.mean() > 0.5 else np.nan)
    data[m] = tuple(map(np.array, (ns, rate, t, cnt)))
plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": INK2,
                     "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False})
fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.8), constrained_layout=True)
for a in ax:
    for x, lab in [(N_LLL, "LLL\nguarantee"), (N_BEST, "best known\n(Exoo)")]:
        a.axvline(x, color=MUTED, lw=1, ls=(0, (3, 3)), zorder=0)
        a.text(x, 1.02, lab, transform=a.get_xaxis_transform(), ha="center", va="bottom", fontsize=7, color=INK2)
    a.set_xlabel("n"); a.set_xlim(8, 43.5); a.set_xticks(range(10, 43, 4))
LABEL_POS = {"restart": (0, 9, "center"), "mt": (5, 0, "left"), "annealing": (6, -9, "left"), "greedy": (6, 7, "left")}
for m, (lab, col, mk) in STYLE.items():
    ns, rate, t, cnt = data[m]
    ax[0].plot(ns, rate, color=col, lw=2, marker=mk, ms=4, label=lab, markeredgecolor="white", markeredgewidth=0.8)
    k = ~np.isnan(t)
    ax[1].plot(ns[k], t[k], color=col, lw=2, marker=mk, ms=4, markeredgecolor="white", markeredgewidth=0.8)
    dx, dy, ha = LABEL_POS[m]
    ax[1].annotate(lab, (ns[k][-1], t[k][-1]), xytext=(dx, dy), textcoords="offset points", va="center", ha=ha, fontsize=7, color=INK)
ax[0].set_ylim(-0.03, 1.03); ax[0].set_ylabel("success rate (budget $10^6$ steps)")
ax[0].set_title("(a) Runs solved", fontsize=9, color=INK, loc="left", pad=22)
ax[0].legend(frameon=False, fontsize=7, loc="lower left")
ax[1].set_yscale("log"); ax[1].set_ylabel("median time per run (ms)")
ax[1].set_title("(b) Median time, where most runs succeed", fontsize=9, color=INK, loc="left", pad=22)
ax[1].set_xlim(8, 50)
fig.savefig(out + ".pdf"); fig.savefig(out + ".png", dpi=200); print("saved", out)
