"""Regression test: refactored integration must reproduce the published v2 PDAC
ranking exactly (same scores, same ranks, AUC 0.753) from the same axis tables.

Run: python tests/test_regression_pdac.py
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))
from multiaxisoncology.config import load_config
from multiaxisoncology.integrate import integrate, validate

ROOT = Path(__file__).parent.parent
V2 = Path("/mnt/results/01_target_nomination/tables_v2")
DATA = Path("/workspace/pdac_ai/data")

cfg = load_config(ROOT / "configs/pdac.yaml")

# axis tables: v2 outputs (dep table uses old column name mean_pdac)
dep = pd.read_csv(V2 / "depmap_pdac_selective_dependencies_v2.csv")
dep = dep.rename(columns={"mean_pdac": "mean_lineage"})

tables = dict(
    de=pd.read_csv("/mnt/results/01_target_nomination/tables/de_paad_vs_normal_limma.csv"),
    surv=pd.read_csv(V2 / "survival_cox_paad.csv"),
    dep=dep,
    ot=pd.read_csv(V2 / "opentargets_tractability_v2.csv"),
    mut=pd.read_csv(V2 / "mutation_frequency.csv", index_col=0).iloc[:, 0],
    cna=pd.read_csv(V2 / "cna_frequencies.csv", index_col=0),
    gtex_stats=pd.read_csv(DATA / "gtex_tissue_expression_stats.csv", index_col=0)
        .rename(columns={"pancreas_log2tpm": "pancreas_log2tpm"}),  # tissue col name matches
    gnomad=pd.read_csv(DATA / "gnomad_constraint.tsv", sep="\t", low_memory=False),
    bench=pd.read_csv(ROOT / "configs/benchmarks/pdac_benchmark.csv"),
)

new = integrate(cfg, tables)
old = pd.read_csv(V2 / "pdac_target_ranking_v2.csv")

assert len(new) == len(old), f"universe size changed: {len(new)} vs {len(old)}"
m = old[["gene", "rank", "score"]].merge(new[["gene", "rank", "score"]], on="gene",
                                         suffixes=("_old", "_new"))
assert len(m) == len(old), "gene sets differ"
assert (m["rank_old"] == m["rank_new"]).all(), "ranks differ"
assert np.allclose(m["score_old"], m["score_new"], atol=1e-10), "scores differ"

metrics = validate(new, cfg)
assert abs(metrics["auc"] - 0.753) < 0.001, f"AUC changed: {metrics['auc']}"
assert metrics["n_benchmark_scored"] == 75

print("REGRESSION PASS: ranking reproduced exactly "
      f"({len(new)} genes, AUC {metrics['auc']:.3f}, p = {metrics['p']:.2e})")
print("top-5:", new.head(5).gene.tolist())
