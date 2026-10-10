"""Integrated multi-axis target scoring (config-driven).

Pre-registered weights come from the disease config and are FROZEN: never tune
them to improve benchmark concordance. The literature benchmark is validation
ONLY and never enters scoring. Missing axis -> neutral 0.5 percentile.
"""
import numpy as np
import pandas as pd
from scipy import stats

from .utils import load_probemap, ens2sym_map, neutral


def integrate(cfg, tables):
    """tables: dict with keys de, surv, dep, ot, mut, cna, gtex_stats, gnomad, bench.

    gtex_stats must have a column '<tissue>_log2tpm' for the config tissue.
    Returns the ranked DataFrame.
    """
    de, sur, dep, ot = tables["de"], tables["surv"], tables["dep"], tables["ot"]
    mut, cna, gtex, gnomad, bench = (tables["mut"], tables["cna"], tables["gtex_stats"],
                                     tables["gnomad"], tables["bench"])
    tissue_col = f"{cfg['gtex_tissue'].lower()}_log2tpm"

    df = pd.DataFrame({"gene": de.gene.tolist()})
    de_i = de.set_index("gene")
    df["logFC"] = df.gene.map(de_i.logFC)

    df["surv_q"] = df.gene.map(sur.set_index("gene").q)

    dep_i = dep.set_index("gene")
    df["dep_mean"] = df.gene.map(dep_i.mean_lineage)
    df["dep_selective"] = df.gene.map(dep_i.selective).fillna(False).astype(bool)

    ot_i = ot.set_index("gene")
    df["sm_tractable"] = df.gene.map(ot_i.tract_small_molecule).fillna(False).astype(float)
    df["known_drugs"] = df.gene.map(ot_i.known_drugs).fillna(0)

    df["mut_freq"] = df.gene.map(mut)
    df["amp_freq"] = df.gene.map(cna.amp_freq)
    df["del_freq"] = df.gene.map(cna.del_freq)

    gtex_i = gtex.drop_duplicates()
    gtex_i.index = gtex_i.index.str.replace(r"\.\d+$", "", regex=True)
    pm = load_probemap(f"{cfg['paths']['data_dir']}/gencode.v23.probemap.tsv")
    gtex_i.index = gtex_i.index.map(ens2sym_map(pm))
    gtex_i = gtex_i.groupby(level=0).max()
    df["tissue_expr"] = df.gene.map(gtex_i[tissue_col])

    gn = gnomad.copy()
    gn["lof.oe_ci.upper"] = pd.to_numeric(gn["lof.oe_ci.upper"], errors="coerce")
    gn = gn[gn["lof.oe_ci.upper"].notna()].sort_values("mane_select", ascending=False)
    gn = gn.drop_duplicates("gene").set_index("gene")
    df["LOEUF"] = df.gene.map(gn["lof.oe_ci.upper"])

    W = cfg["weights"]
    df["score"] = (W["mut"] * neutral(df.mut_freq.fillna(0.0))
                   + W["amp"] * neutral(df.amp_freq.fillna(0.0))
                   + W["dep"] * neutral(-df.dep_mean.fillna(-0.1))
                   + W["logFC"] * neutral(df.logFC)
                   + W["surv"] * neutral(-np.log10(df.surv_q))
                   + W["loeuf"] * neutral(df.LOEUF)
                   + W["tissue"] * (1 - neutral(df.tissue_expr))
                   + W["tract"] * (df.sm_tractable + np.minimum(df.known_drugs, 10) / 10) / 1.5)

    df = df.sort_values("score", ascending=False).reset_index(drop=True)
    df["rank"] = df.index + 1
    df["in_benchmark"] = df.gene.isin(set(bench.gene))
    return df


def rank_auc(ranks_val, n_all):
    n_val = len(ranks_val)
    return ((np.sum(n_all - ranks_val + 1) - n_val * (n_val + 1) / 2)
            / (n_val * (n_all - n_val)))


def validate(df, cfg=None):
    """Benchmark AUC, enrichment, ablation, weight sensitivity. Returns dict."""
    ranks_val = df.loc[df.in_benchmark, "rank"].to_numpy()
    n_val, n_all = len(ranks_val), len(df)
    auc = rank_auc(ranks_val, n_all)
    u, p = stats.mannwhitneyu(ranks_val, df.loc[~df.in_benchmark, "rank"],
                              alternative="less")
    out = {"n_benchmark_scored": int(n_val), "auc": float(auc), "p": float(p)}
    for N in (50, 100, 250):
        hit = int(df.head(N).in_benchmark.sum())
        exp = N * n_val / n_all
        out[f"top{N}_hits"] = hit
        out[f"top{N}_fold"] = float(hit / max(exp, 0.01))

    W = cfg["weights"] if cfg else dict(mut=0.20, amp=0.10, dep=0.20, logFC=0.10,
                                        surv=0.10, loeuf=0.10, tissue=0.10, tract=0.10)
    base_axes = {
        "mut": neutral(df.mut_freq.fillna(0.0)),
        "amp": neutral(df.amp_freq.fillna(0.0)),
        "dep": neutral(-df.dep_mean.fillna(-0.1)),
        "logFC": neutral(df.logFC),
        "surv": neutral(-np.log10(df.surv_q)),
        "loeuf": neutral(df.LOEUF),
        "tissue": 1 - neutral(df.tissue_expr),
        "tract": (df.sm_tractable + np.minimum(df.known_drugs, 10) / 10) / 1.5,
    }

    def auc_from_weights(wts):
        s = sum(w * base_axes[k] for k, w in wts.items())
        r = (-s).rank(method="first").to_numpy()
        rv = r[df.in_benchmark.to_numpy()]
        return rank_auc(rv, n_all)

    out["ablation"] = {f"drop_{k}": float(auc_from_weights(
        {a: w / (1 - W[a]) for a, w in W.items() if a != k})) for k in base_axes}
    out["ablation"]["full"] = float(auc_from_weights(W))

    rng = np.random.default_rng(42)
    keys = list(W.keys())
    dirich = rng.dirichlet([2, 1, 2, 1, 1, 1, 1, 1], size=200)
    sens = np.array([auc_from_weights(dict(zip(keys, row))) for row in dirich])
    out["sens_median"] = float(np.median(sens))
    out["sens_ci"] = [float(np.percentile(sens, 5)), float(np.percentile(sens, 95))]
    return out
