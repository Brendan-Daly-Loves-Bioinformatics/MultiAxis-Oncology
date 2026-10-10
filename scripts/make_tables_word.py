"""Render Table 1 (data sources) and Table 2 (weights) as high-res PNGs for Word,
plus standalone LaTeX sources for Overleaf compilation."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams["font.family"] = ["Liberation Sans", "DejaVu Sans"]

OUT = "/mnt/results/oncotarget"

# ---------- Table 1: data sources ----------
cols = ["Dataset", "Source", "Version", "Notes"]
rows = [
    ["Tumor expression", "Toil TCGA matrix (UCSC Xena)", "Toil recompute", "log2(TPM+0.001)"],
    ["Normal expression", "Toil GTEx matrix (UCSC Xena)", "Toil recompute", "pancreas n=167; lung n=288"],
    ["Gene annotation", "GENCODE probemap", "v23", "version-stripped Ensembl IDs"],
    ["Somatic mutations", "GDC Xena hub", "GDC", "xenaPython sparse API"],
    ["Copy number (ASCAT2)", "GDC Xena hub", "GDC", "amp: CN>=4; del: CN<=1"],
    ["Survival / clinical", "GDC Xena hub", "GDC", "xenaPython dataset_fetch"],
    ["CRISPR dependency", "DepMap (Bioconductor depmap)", "22Q2", "Chronos scores"],
    ["Constraint (LOEUF)", "gnomAD constraint metrics", "v4.0", "MANE Select preferred"],
    ["Tractability", "Open Targets Platform GraphQL", "2025 access", "target-level, disease-agnostic"],
]
def table_fig(cols, rows, title, path, col_widths=None):
    fig, ax = plt.subplots(figsize=(9.5, 0.42 * (len(rows) + 1) + 0.6))
    ax.axis("off")
    tbl = ax.table(cellText=rows, colLabels=cols, loc="upper center",
                   cellLoc="left", colWidths=col_widths)
    tbl.auto_set_font_size(False); tbl.set_fontsize(9)
    for (r, c), cell in tbl.get_celld().items():
        cell.set_edgecolor("#bbbbbb")
        if r == 0:
            cell.set_facecolor("#e8eef7"); cell.set_text_props(weight="bold")
        elif r % 2 == 0:
            cell.set_facecolor("#f5f7fa")
    ax.set_title(title, fontsize=11, loc="left", pad=12)
    fig.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)

table_fig(cols, rows, "Table 1. Data sources and versions.",
          f"{OUT}/table1_data_sources.png")

# ---------- Table 2: weights ----------
cols2 = ["Evidence axis", "Weight", "Direction"]
rows2 = [
    ["Somatic mutation frequency", "0.20", "higher = stronger candidate"],
    ["CRISPR dependency (lineage mean Chronos)", "0.20", "more dependent = stronger"],
    ["Amplification frequency", "0.10", "higher = stronger"],
    ["Tumor logFC (tumor vs normal)", "0.10", "higher = stronger"],
    ["Survival association (-log10 q)", "0.10", "more significant = stronger"],
    ["LOEUF (gnomAD constraint)", "0.10", "less constrained = safer"],
    ["Low normal-tissue expression", "0.10", "lower expression = safer"],
    ["Tractability (SM flag + known drugs)", "0.10", "more tractable = stronger"],
]
table_fig(cols2, rows2, "Table 2. Pre-registered axis weights (frozen across diseases).",
          f"{OUT}/table2_weights.png")

# ---------- LaTeX standalone sources ----------
tex1 = r"""\documentclass{standalone}
\usepackage{booktabs}
\begin{document}
\begin{tabular}{p{3.2cm}p{4.6cm}p{2.2cm}p{4.2cm}}
\toprule
Dataset & Source & Version & Notes \\
\midrule
Tumor expression & Toil TCGA matrix (UCSC Xena) & Toil recompute & log2(TPM+0.001) \\
Normal expression & Toil GTEx matrix (UCSC Xena) & Toil recompute & pancreas $n{=}167$; lung $n{=}288$ \\
Gene annotation & GENCODE probemap & v23 & version-stripped Ensembl IDs \\
Somatic mutations & GDC Xena hub & GDC & xenaPython sparse API \\
Copy number (ASCAT2) & GDC Xena hub & GDC & amp: CN$\geq$4; del: CN$\leq$1 \\
Survival / clinical & GDC Xena hub & GDC & xenaPython dataset\_fetch \\
CRISPR dependency & DepMap (Bioconductor depmap) & 22Q2 & Chronos scores \\
Constraint (LOEUF) & gnomAD constraint metrics & v4.0 & MANE Select preferred \\
Tractability & Open Targets Platform GraphQL & 2025 access & target-level, disease-agnostic \\
\bottomrule
\end{tabular}
\end{document}
"""
tex2 = r"""\documentclass{standalone}
\usepackage{booktabs}
\begin{document}
\begin{tabular}{p{7.2cm}p{1.6cm}p{5.4cm}}
\toprule
Evidence axis & Weight & Direction \\
\midrule
Somatic mutation frequency & 0.20 & higher = stronger candidate \\
CRISPR dependency (lineage mean Chronos) & 0.20 & more dependent = stronger \\
Amplification frequency & 0.10 & higher = stronger \\
Tumor logFC (tumor vs normal) & 0.10 & higher = stronger \\
Survival association ($-\log_{10} q$) & 0.10 & more significant = stronger \\
LOEUF (gnomAD constraint) & 0.10 & less constrained = safer \\
Low normal-tissue expression & 0.10 & lower expression = safer \\
Tractability (SM flag + known drugs) & 0.10 & more tractable = stronger \\
\bottomrule
\end{tabular}
\end{document}
"""
with open(f"{OUT}/table1_data_sources.tex", "w") as f: f.write(tex1)
with open(f"{OUT}/table2_weights.tex", "w") as f: f.write(tex2)
print("tables written: PNG (Word-ready, 300dpi) + .tex (Overleaf-ready)")
