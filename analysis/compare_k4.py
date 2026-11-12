"""Compare methods on one k: success rate and median time vs n.
Usage: python analysis/compare_k4.py OUT_PREFIX results/k4_mt.csv results/k4_restart.csv ...
"""
import csv, sys, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

INK, INK2, MUTED = "#0b0b0b", "#52514e", "#9a9893"
# fixed categorical order (validated: adjacent CVD dE >= 9.1, normal >= 22.9 on light);
# low-contrast slots are relieved by direct labels + distinct markers.
STYLE = {  # method prefix -> (label, colour, marker)
    "mt":        ("Moser–Tardos",      "#2a78d6", "o"),
    "restart":   ("Random restart",    "#eb6834", "s"),
    "annealing": ("Simulated annealing", "#1baf7a", "^"),
    "greedy":    ("Greedy local search", "#eda100", "D"),
}
N_LLL, N_BEST = 5, 17
LABEL_POS = {"restart": (-6, "right")}  # avoid collision with the MT end label

out, files = sys.argv[1], sys.argv[2:]
data = {}
for f in files:
    rows = list(csv.DictReader(open(f)))
    key = rows[0]["method"].split("-")[0]
    ns = sorted({int(r["n"]) for r in rows})
    rate, med_t = [], []
    for n in ns:
        R = [r for r in rows if int(r["n"]) == n]
        ok = np.array([r["success"] == "True" for r in R])
        t = np.array([float(r["time_ms"]) for r in R])
        rate.append(ok.mean())
        med_t.append(np.median(t) if ok.mean() > 0.5 else np.nan)
    data[key] = (np.array(ns), np.array(rate), np.array(med_t))

order = [k for k in ["mt", "restart", "annealing", "greedy"] if k in data]
plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
                     "xtick.color": INK2, "ytick.color": INK2,
                     "axes.spines.top": False, "axes.spines.right": False})
fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.8), constrained_layout=True)
for a in ax:
    for x, lab in [(N_LLL, "LLL\nguarantee"), (N_BEST, "best known\n(Paley)")]:
        a.axvline(x, color=MUTED, lw=1, ls=(0, (3, 3)), zorder=0)
        a.text(x, 1.02, lab, transform=a.get_xaxis_transform(), ha="center",
               va="bottom", fontsize=7, color=INK2)
    a.set_xlabel("n"); a.set_xlim(4.5, 18.5)

for k in order:
    lab, col, mk = STYLE[k]
    ns, rate, med_t = data[k]
    ax[0].plot(ns, rate, color=col, lw=2, marker=mk, ms=4.5, label=lab,
               markeredgecolor="white", markeredgewidth=0.8)
    m = ~np.isnan(med_t)
    ax[1].plot(ns[m], med_t[m], color=col, lw=2, marker=mk, ms=4.5,
               markeredgecolor="white", markeredgewidth=0.8)
    xe, ye = ns[m][-1], med_t[m][-1]
    dx, ha = LABEL_POS.get(k, (5, "left"))
    ax[1].annotate(lab, (xe, ye), xytext=(dx, 0), textcoords="offset points",
                   va="center", ha=ha, fontsize=7, color=INK)

ax[0].set_ylim(-0.03, 1.03)
ax[0].set_ylabel("success rate (budget $10^6$ steps)")
ax[0].set_title("(a) Runs solved", fontsize=9, color=INK, loc="left", pad=22)
ax[0].legend(frameon=False, fontsize=7, loc="lower left")
ax[1].set_yscale("log")
ax[1].set_ylabel("median time per run (ms)")
ax[1].set_title("(b) Median time, where most runs succeed", fontsize=9, color=INK, loc="left", pad=22)
ax[1].set_xlim(4.5, 20.5)
for a in ax:
    a.set_xticks(range(5, 18, 2))

os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
fig.savefig(out + ".pdf"); fig.savefig(out + ".png", dpi=200)
print("saved", out)
