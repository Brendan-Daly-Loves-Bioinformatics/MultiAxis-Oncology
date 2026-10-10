"""Side-by-side PDAC vs LUAD validation figure (PNG + SVG)."""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams["font.family"] = ["Liberation Sans", "DejaVu Sans"]
plt.rcParams["svg.fonttype"] = "none"

pdac = pd.read_csv("/mnt/results/01_target_nomination/tables_v2/pdac_target_ranking_v2.csv")
luad = pd.read_csv("/workspace/MultiAxis-Oncology/work/luad/luad_target_ranking.csv")

fig, axes = plt.subplots(1, 3, figsize=(12.5, 3.8))

# Panel A: AUC comparison
ax = axes[0]
vals = [0.753, 0.829]
bars = ax.bar(["PDAC\n(75/91 bench)", "LUAD\n(72/76 bench)"], vals,
              color=["#0279EE", "#FF9400"], width=0.55)
for b, v, p in zip(bars, vals, ["p = 1.9e-14", "p = 2.5e-22"]):
    ax.text(b.get_x() + b.get_width() / 2, v + 0.015, f"{v:.3f}\n{p}",
            ha="center", fontsize=9)
ax.set_ylim(0.5, 1.0)
ax.set_ylabel("benchmark rank-AUC")
ax.axhline(0.5, ls="--", c="#999999", lw=0.8)
ax.set_title("A. Validation, frozen weights", fontsize=10, loc="left")

# Panel B: top-N enrichment curves
ax = axes[1]
for df, label, col in [(pdac, "PDAC", "#0279EE"), (luad, "LUAD", "#FF9400")]:
    n_val, n_all = df.in_benchmark.sum(), len(df)
    Ns = np.arange(10, 500, 10)
    fold = [df.head(N).in_benchmark.sum() / (N * n_val / n_all) for N in Ns]
    ax.plot(Ns, fold, lw=2, label=label, color=col)
ax.axhline(1, ls="--", c="#999999", lw=0.8)
ax.set_xlabel("top-N genes")
ax.set_ylabel("benchmark enrichment (fold over chance)")
ax.legend(frameon=False, fontsize=9)
ax.set_title("B. Enrichment of literature targets", fontsize=10, loc="left")

# Panel C: shared benchmark targets, rank in each cancer
ax = axes[2]
shared = sorted(set(pdac[pdac.in_benchmark].gene) & set(luad[luad.in_benchmark].gene))
r1 = [pdac.loc[pdac.gene == g, "rank"].iloc[0] for g in shared]
r2 = [luad.loc[luad.gene == g, "rank"].iloc[0] for g in shared]
lim = 300
r1c = np.minimum(r1, lim); r2c = np.minimum(r2, lim)
ax.scatter(r1c, r2c, c="#0279EE", s=28, alpha=0.75)
for g, x, y in zip(shared, r1c, r2c):
    if g in ("KRAS", "EGFR", "AURKA", "TOP2A", "PIK3CA", "TYMS", "ERBB2", "BRCA1"):
        ax.annotate(g, (x, y), fontsize=8, xytext=(4, 4), textcoords="offset points")
ax.plot([1, lim], [1, lim], ls="--", c="#999999", lw=0.8)
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlim(lim, 1); ax.set_ylim(lim, 1)
ax.set_xlabel("rank in PDAC")
ax.set_ylabel("rank in LUAD")
ax.set_title("C. Benchmark target ranks across cancers", fontsize=10, loc="left")

fig.tight_layout()
fig.savefig("/mnt/results/oncotarget/fig_validation_pdac_luad.png", dpi=200)
fig.savefig("/mnt/results/oncotarget/fig_validation_pdac_luad.svg")
print("figure saved")
