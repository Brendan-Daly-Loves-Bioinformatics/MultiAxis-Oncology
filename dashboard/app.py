"""Oncotarget dashboard — Streamlit app over the validated ranking engine.

Reads precomputed pipeline outputs (no live fetching). Disease-agnostic:
add a cancer type by dropping its ranking CSV + metrics JSON into data/.
"""
from pathlib import Path
import json

import numpy as np
import pandas as pd
import streamlit as st

st.set_page_config(page_title="MultiAxis Oncology — cancer target prioritization",
                   page_icon="🧬", layout="wide")

DATA = Path(__file__).parent / "data"

DISEASES = {
    "PDAC — Pancreatic ductal adenocarcinoma": {
        "ranking": "pdac_target_ranking.csv", "metrics": "pdac_metrics.json",
        "short": "PDAC",
    },
    "LUAD — Lung adenocarcinoma": {
        "ranking": "luad_target_ranking.csv", "metrics": "luad_metrics.json",
        "short": "LUAD",
    },
}

AXIS_LABELS = {
    "mut_freq": ("Somatic mutation frequency", "fraction of patients with a damaging mutation"),
    "amp_freq": ("Amplification frequency", "fraction of tumors with total CN >= 4"),
    "dep_mean": ("CRISPR dependency (lineage mean Chronos)",
                 "more negative = stronger dependency"),
    "logFC": ("Tumor vs normal logFC", "expression fold-change, log2"),
    "surv_q": ("Survival association (q)", "genome-wide Cox, BH-adjusted"),
    "LOEUF": ("LOEUF (gnomAD constraint)", "lower = more LOF-constrained = safety concern"),
    "tissue_expr": ("Normal-tissue expression (log2 TPM)",
                    "lower = better on-target safety window"),
    "sm_tractable": ("Small-molecule tractable", "Open Targets tractability"),
    "known_drugs": ("Known drugs / clinical candidates", "Open Targets count"),
}


@st.cache_data
def load_ranking(path):
    df = pd.read_csv(path)
    # harmonize column names across disease versions
    df = df.rename(columns={"dep_mean_pdac": "dep_mean", "pancreas_expr": "tissue_expr"})
    return df


@st.cache_data
def load_metrics(path):
    return json.loads(Path(path).read_text())


def fmt(v, kind):
    if pd.isna(v):
        return "—"
    if kind == "pct":
        return f"{v:.3f}"
    if kind == "q":
        return f"{v:.2e}" if v > 0 else "—"
    if kind == "num":
        return f"{v:,.0f}"
    return f"{v:.2f}"


def evidence_card(gene, df):
    row = df[df.gene == gene]
    if row.empty:
        st.warning(f"{gene} not in the scored universe.")
        return
    r = row.iloc[0]
    rank, n = int(r["rank"]), len(df)
    st.subheader(f"{gene}  ·  rank {rank} of {n:,}")
    if r["in_benchmark"]:
        st.success("Literature benchmark target — established evidence for this disease.")
    pct = {k: (df[k].rank(pct=True)[row.index[0]]) for k in AXIS_LABELS if k in df.columns}

    c1, c2, c3 = st.columns(3)
    with c1:
        st.metric("Score", f"{r['score']:.3f}")
        st.metric("Percentile", f"{(n - rank + 1) / n * 100:.1f}")
    with c2:
        st.metric("Mutation freq", fmt(r.get("mut_freq"), "pct"))
        st.metric("Amplification freq", fmt(r.get("amp_freq"), "pct"))
    with c3:
        st.metric("Dependency (Chronos)", fmt(r.get("dep_mean"), "num"))
        st.metric("Tumor logFC", fmt(r.get("logFC"), "num"))

    st.markdown("**Axis detail**")
    detail = []
    for k, (label, note) in AXIS_LABELS.items():
        if k not in r.index:
            continue
        v = r[k]
        p = pct.get(k)
        ptxt = f"{p * 100:.0f}th pct" if p is not None and not pd.isna(p) else "—"
        detail.append({"Axis": label, "Value": fmt(v, "q" if k == "surv_q" else "pct"),
                       "Universe percentile": ptxt, "Note": note})
    st.dataframe(pd.DataFrame(detail), hide_index=True, use_container_width=True)

    # deterministic, template-based interpretation (no LLM involved)
    st.markdown("**Automated summary** (template-generated from the values above)")
    bits = []
    if not pd.isna(r.get("mut_freq")) and r["mut_freq"] > 0.10:
        bits.append(f"frequently mutated ({r['mut_freq'] * 100:.0f}% of patients)")
    if not pd.isna(r.get("amp_freq")) and r["amp_freq"] > 0.10:
        bits.append(f"frequently amplified ({r['amp_freq'] * 100:.0f}%)")
    if not pd.isna(r.get("dep_mean")) and r["dep_mean"] < -0.5:
        bits.append("strong lineage CRISPR dependency")
    if not pd.isna(r.get("logFC")) and r["logFC"] > 1:
        bits.append(f"overexpressed in tumors (logFC {r['logFC']:.1f})")
    if r.get("sm_tractable"):
        kd = int(r.get("known_drugs", 0))
        bits.append(f"small-molecule tractable" + (f" with {kd} known drugs/candidates" if kd else ""))
    if not pd.isna(r.get("tissue_expr")) and r["tissue_expr"] > 3:
        bits.append("expressed in normal tissue (safety penalty applied)")
    if not bits:
        bits.append("no dominant single-axis signal; ranked on combined evidence")
    st.write(f"**{gene}** ranks {rank} of {n:,}: " + ", ".join(bits) + ".")


def validation_panel(metrics):
    st.subheader("Benchmark validation")
    c1, c2, c3 = st.columns(3)
    c1.metric("Rank-AUC", f"{metrics['auc']:.3f}")
    c2.metric("Mann-Whitney p", f"{metrics['p']:.2e}")
    c3.metric("Benchmark targets scored", metrics["n_benchmark_scored"])
    en = pd.DataFrame({
        "Top-N": [50, 100, 250],
        "Hits": [metrics.get(f"top{N}_hits") for N in (50, 100, 250)],
        "Fold over chance": [f"{metrics.get(f'top{N}_fold'):.1f}x" for N in (50, 100, 250)],
    })
    st.dataframe(en, hide_index=True, use_container_width=True)
    st.markdown("**Axis ablation** (AUC when each axis is dropped)")
    abl = pd.DataFrame({"AUC when dropped": metrics["ablation"]}).T
    st.dataframe(abl, use_container_width=True)
    st.caption(f"Weight sensitivity (200 Dirichlet draws): median {metrics['sens_median']:.3f}, "
               f"5–95% [{metrics['sens_ci'][0]:.3f}, {metrics['sens_ci'][1]:.3f}]")


# ---------------- layout ----------------
st.title("MultiAxis Oncology")
st.caption("Transparent, benchmark-validated multi-axis cancer target prioritization. "
           "Pre-registered weights; literature benchmarks used for validation only.")

disease = st.sidebar.selectbox("Cancer type", list(DISEASES))
cfg = DISEASES[disease]
df = load_ranking(DATA / cfg["ranking"])
metrics = load_metrics(DATA / cfg["metrics"])

tab_rank, tab_gene, tab_val, tab_about = st.tabs(
    ["Ranked targets", "Per-gene evidence", "Validation", "About"])

with tab_rank:
    c1, c2, c3 = st.columns([2, 1, 1])
    q = c1.text_input("Search gene", "")
    top_n = c2.selectbox("Show top", [25, 50, 100, 250, "all"], index=1)
    bench_only = c3.checkbox("Benchmark targets only")
    view = df
    if q:
        view = view[view.gene.str.upper().str.contains(q.upper())]
    if bench_only:
        view = view[view.in_benchmark]
    if top_n != "all":
        view = view.head(int(top_n))
    show_cols = ["rank", "gene", "score", "mut_freq", "amp_freq", "dep_mean", "logFC",
                 "surv_q", "LOEUF", "tissue_expr", "sm_tractable", "known_drugs",
                 "in_benchmark"]
    show_cols = [c for c in show_cols if c in view.columns]
    st.dataframe(view[show_cols], hide_index=True, use_container_width=True, height=520)
    st.download_button("Download full ranking (CSV)",
                       df.to_csv(index=False).encode(),
                       file_name=f"multiaxisoncology_{cfg['short'].lower()}_ranking.csv")

with tab_gene:
    gene = st.selectbox(
        "Gene",
        df["gene"].tolist(),
        index=0,
        format_func=lambda g: (
            f"{g} (#{int(df.loc[df['gene'] == g, 'rank'].iloc[0])})"
        )
    )
    evidence_card(gene, df)

with tab_val:
    validation_panel(metrics)

with tab_about:
    st.markdown("""
**What this is.** A transparent, multi-axis target-prioritization framework for cancer
research. Genes are scored on eight evidence axes — somatic mutation frequency, copy-number
amplification, lineage CRISPR dependency, tumor-vs-normal expression, survival association,
genetic constraint (LOEUF), normal-tissue expression, and druggability/tractability — with
**pre-registered weights frozen across cancer types**. Literature benchmarks (approved drugs,
clinical-trial targets, validated preclinical targets) are used **for validation only** and
never enter scoring.

**What this is not.** Not clinical decision support. Outputs are computational hypotheses
from public data, intended to prioritize research targets.

**Interpretation.** Every rank is decomposable: the per-gene evidence tab shows each axis
value, its percentile in the scored universe, and a template-generated summary. No language
model participates in scoring.
""")
    st.info("Optional LLM narration layer: planned, disabled by default, and — when enabled — "
            "grounded strictly in the pipeline's own outputs. It will never alter rankings.")
