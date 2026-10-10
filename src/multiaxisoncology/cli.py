"""MultiAxis Oncology CLI.

Usage:
  python -m multiaxisoncology.cli run <config.yaml> [--skip-fetch] [--skip-axes]
  python -m multiaxisoncology.cli validate <config.yaml> <ranking.csv>
"""
import argparse
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from .config import load_config
from .integrate import integrate, validate
from .fetch import expression as fx
from .fetch import gtex_stats as fg
from .fetch import mutations_cna as fm
from .fetch import survival as fs
from .fetch import opentargets as fo
from .utils import load_probemap, sym2ens_map

AXES_DIR = Path(__file__).parent / "axes"


def rscript(script, *args):
    cmd = ["Rscript", str(AXES_DIR / script)] + [str(a) for a in args]
    print("+", " ".join(cmd))
    subprocess.run(cmd, check=True)


def load_tables(cfg):
    work = cfg["paths"]["work_dir"]
    data = cfg["paths"]["data_dir"]
    code = cfg["disease_code"].lower()
    tissue = cfg["gtex_tissue"].lower()
    ot_path = cfg.get("opentargets_table") or f"{data}/opentargets_tractability_genome.csv"
    return dict(
        de=pd.read_csv(f"{work}/de_{code}_vs_normal_limma.csv"),
        surv=pd.read_csv(f"{work}/survival_cox_{code}.csv"),
        dep=pd.read_csv(f"{work}/depmap_{code}_selective_dependencies.csv"),
        ot=pd.read_csv(ot_path),
        mut=pd.read_csv(f"{work}/mutation_frequency.csv", index_col=0).iloc[:, 0],
        cna=pd.read_csv(f"{work}/cna_frequencies.csv", index_col=0),
        gtex_stats=pd.read_csv(f"{data}/gtex_tissue_expression_stats_{tissue}.csv",
                               index_col=0),
        gnomad=pd.read_csv(f"{data}/gnomad_constraint.tsv", sep="\t", low_memory=False),
        bench=pd.read_csv(cfg["benchmark"]),
    )


def cmd_run(cfg, skip_fetch=False, skip_axes=False):
    work = Path(cfg["paths"]["work_dir"]); work.mkdir(parents=True, exist_ok=True)
    out = Path(cfg["paths"]["out_dir"]); out.mkdir(parents=True, exist_ok=True)
    data = cfg["paths"]["data_dir"]
    code = cfg["disease_code"]; lc = code.lower()
    tissue = cfg["gtex_tissue"]

    # shared, disease-agnostic inputs
    gnomad_path = f"{data}/gnomad_constraint.tsv"
    ot_path = cfg.get("opentargets_table") or f"{data}/opentargets_tractability_genome.csv"
    if not Path(ot_path).exists():
        print("Open Targets genome-wide table missing; querying (once, reused later)...")
        pm = load_probemap(f"{data}/gencode.v23.probemap.tsv")
        s2e = sym2ens_map(pm)
        de_universe = pd.read_csv(f"{work}/de_{lc}_vs_normal_limma.csv")
        pairs = [(g, s2e.get(g)) for g in de_universe.gene]
        pairs = [(g, e) for g, e in pairs if e]
        fo.query_opentargets(pairs, ot_path)

    if not skip_fetch:
        expr_base = fx.slice_expression(cfg)
        if not Path(f"{data}/gtex_tissue_expression_stats_{tissue.lower()}.csv").exists():
            fg.gtex_tissue_stats(cfg)
        fs.fetch_survival_clinical(cfg)

    if not skip_axes:
        expr_csv = f"{work}/{lc}_gtex_{tissue.lower()}_tpm.csv.gz"
        if not Path(expr_csv).exists():
            pd.read_parquet(f"{work}/{lc}_gtex_{tissue.lower()}_tpm.parquet") \
              .to_csv(expr_csv)
        rscript("de_limma.R", expr_csv, f"{data}/gencode.v23.probemap.tsv",
                f"{work}/de_{lc}_vs_normal_limma.csv")
        rscript("survival_cox.R", expr_csv, f"{work}/{code}_survival.csv",
                f"{work}/{code}_clinical.csv", f"{data}/gencode.v23.probemap.tsv",
                f"{work}/survival_cox_{lc}.csv")
        rscript("depmap_dependency.R", f"{data}/depmap_crispr_22Q2.rds",
                f"{data}/depmap_metadata_22Q2.rds", cfg["depmap_lineage"],
                cfg.get("depmap_lineage_subtype", "NA"),
                cfg.get("depmap_lineage_sub_subtype", "NA"),
                f"{work}/depmap_{lc}_selective_dependencies.csv")
        de = pd.read_csv(f"{work}/de_{lc}_vs_normal_limma.csv")
        fm.run(cfg, de.gene.tolist())

    tables = load_tables(cfg)
    ranked = integrate(cfg, tables)
    ranked.to_csv(f"{out}/{lc}_target_ranking.csv", index=False)

    metrics = validate(ranked, cfg)
    import json
    with open(f"{out}/validation_metrics.json", "w") as f:
        json.dump(metrics, f, indent=2)

    print(f"\n=== {cfg['disease_name']} ===")
    print(f"AUC {metrics['auc']:.3f} | p = {metrics['p']:.2e} "
          f"| {metrics['n_benchmark_scored']} benchmark targets scored")
    for N in (50, 100, 250):
        print(f"top-{N}: {metrics[f'top{N}_hits']} hits ({metrics[f'top{N}_fold']:.1f}x)")
    print("TOP 20:")
    print(ranked.head(20)[["rank", "gene", "score", "in_benchmark"]].to_string(index=False))
    return ranked, metrics


def main():
    ap = argparse.ArgumentParser(prog="multiaxisoncology")
    sub = ap.add_subparsers(dest="cmd", required=True)
    run = sub.add_parser("run")
    run.add_argument("config")
    run.add_argument("--skip-fetch", action="store_true")
    run.add_argument("--skip-axes", action="store_true")
    args = ap.parse_args()
    cfg = load_config(args.config)
    if args.cmd == "run":
        cmd_run(cfg, skip_fetch=args.skip_fetch, skip_axes=args.skip_axes)


if __name__ == "__main__":
    main()
