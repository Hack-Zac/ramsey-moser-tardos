"""Summarise a k=4 sweep: table + 3-panel figure for the report.
Usage: python analysis/k4_summary.py results/k4_T1e6.csv figures/k4
"""
import csv, sys, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

INK, INK2, MUTED, BLUE = "#0b0b0b", "#52514e", "#9a9893", "#2a78d6"
N_LLL, N_BEST = 5, 17   # LLL guarantee and largest n with a valid colouring (R(4,4)=18)

path, out = sys.argv[1], sys.argv[2]
rows = list(csv.DictReader(open(path)))
ns = sorted({int(r["n"]) for r in rows})
T = int(rows[0]["budget"])

summary = []
for n in ns:
    R = [r for r in rows if int(r["n"]) == n]
    ok = np.array([r["success"] == "True" for r in R])
    steps = np.array([int(r["steps"]) for r in R])
    fails_bad = np.array([int(r["min_bad_seen"]) for r, s in zip(R, ok) if not s])
    rate = ok.mean()
    # median over all runs is defined only if more than half succeeded (failures are >= T)
    med = float(np.median(steps)) if rate > 0.5 else None
    q1, q3 = (np.percentile(steps, [25, 75]) if rate > 0.75 else (None, None))
    summary.append(dict(n=n, runs=len(R), rate=rate, median=med, q1=q1, q3=q3,
                        fail_min=fails_bad.min() if len(fails_bad) else None,
                        fail_med=float(np.median(fails_bad)) if len(fails_bad) else None,
                        fail_max=fails_bad.max() if len(fails_bad) else None))

print(f"budget T = {T:,}")
print(f"{'n':>3} {'solved':>8} {'median steps':>13} {'x prev':>7} {'fewest K4s on fails (min/med/max)':>34}")
prev = None
for s in summary:
    med = s["median"]
    ratio = f"{med/prev:6.1f}x" if (med and prev) else ""
    fb = "" if s["fail_med"] is None else f"{s['fail_min']}/{s['fail_med']:.0f}/{s['fail_max']}"
    print(f"{s['n']:>3} {int(s['rate']*s['runs']):>4}/{s['runs']:<3} "
          f"{'' if med is None else f'{med:,.0f}':>13} {ratio:>7} {fb:>34}")
    prev = med if med else prev

# ---------------- figure ----------------
plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
                     "xtick.color": INK2, "ytick.color": INK2, "axes.spines.top": False,
                     "axes.spines.right": False})
fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.5), constrained_layout=True)

def refs(a, top_label=True):
    for x, lab in [(N_LLL, "LLL\nguarantee"), (N_BEST, "best known\n(Paley)")]:
        a.axvline(x, color=MUTED, lw=1, ls=(0, (3, 3)), zorder=0)
        if top_label:
            a.text(x, 1.02, lab, transform=a.get_xaxis_transform(), ha="center",
                   va="bottom", fontsize=7, color=INK2)

# (a) success rate
x = [s["n"] for s in summary]; y = [s["rate"] for s in summary]
ax[0].plot(x, y, color=BLUE, lw=2, marker="o", ms=4)
ax[0].set_ylim(-0.03, 1.03); ax[0].set_xlabel("n"); ax[0].set_ylabel(f"success rate ($T = 10^{{{int(round(np.log10(T)))}}}$)")
ax[0].set_title("(a) Runs solved", fontsize=9, color=INK, loc="left", pad=22)
refs(ax[0])

# (b) median resamplings (log), with IQR
pts = [s for s in summary if s["median"] is not None]
xb = np.array([s["n"] for s in pts]); yb = np.array([s["median"] + 1 for s in pts])
lo = np.array([(s["q1"] if s["q1"] is not None else s["median"]) + 1 for s in pts])
hi = np.array([(s["q3"] if s["q3"] is not None else s["median"]) + 1 for s in pts])
ax[1].errorbar(xb, yb, yerr=[yb - lo, hi - yb], color=BLUE, lw=2, marker="o", ms=4,
               elinewidth=1, capsize=0)
ax[1].set_yscale("log"); ax[1].set_xlabel("n"); ax[1].set_ylabel("1 + median resamplings")
ax[1].set_xlim(min(ns) - 0.5, max(ns) + 0.5)
ax[1].axhline(T, color=MUTED, lw=1, ls=(0, (1, 2)))
ax[1].text(min(ns) + 0.6, T * 0.35, "budget $T$", ha="left", va="top", fontsize=7, color=INK2)
ax[1].set_title("(b) Cost (median, IQR)", fontsize=9, color=INK, loc="left", pad=22)
refs(ax[1])

# (c) closest approach on failed runs
pf = [s for s in summary if s["fail_med"] is not None]
xc = np.array([s["n"] for s in pf]); yc = np.array([s["fail_med"] for s in pf])
ax[2].errorbar(xc, yc, yerr=[yc - [s["fail_min"] for s in pf], [s["fail_max"] for s in pf] - yc],
               color=BLUE, lw=2, marker="o", ms=4, elinewidth=1, capsize=0)
ax[2].set_xlabel("n"); ax[2].set_ylabel("fewest one-colour $K_4$'s seen")
ax[2].set_ylim(bottom=0); ax[2].set_xlim(min(ns) - 0.5, max(ns) + 0.5)
ax[2].set_title("(c) Closest approach (failures)", fontsize=9, color=INK, loc="left", pad=22)
refs(ax[2])

os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
fig.savefig(out + ".pdf"); fig.savefig(out + ".png", dpi=200)
print("saved", out + ".pdf/.png")
