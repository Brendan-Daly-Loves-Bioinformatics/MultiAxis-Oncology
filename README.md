# MultiAxis Oncology

An open-source, benchmark-validated framework for multi-axis cancer target
prioritization. It scores every protein-coding gene for a given cancer type
across eight evidence axes; somatic mutation frequency, CRISPR dependency,
amplification, differential expression, survival association, germline
constraint (LOEUF), tissue expression, and druggability. It then combines them
into a single transparent, pre-registered ranking.

**What this is not:** This is not machine learning or artificial intelligence of 
any kind. The engine is a deterministic, weighted-percentile scoring function 
with frozen weights and a pre-registered validation protocol. Every intermediate
value is inspectable in the output tables. This is for maximum reproducibility
and accuracy.

## Validation

Benchmarks are literature-curated gene sets (91 PDAC, 75 LUAD) used only for
validation, not for tuning. Weights were frozen before any ranking was
computed and are shared across cancer types.

| Cancer | Genes scored | Benchmark AUC | p | Top-50 enrichment |
|---|---|---|---|---|
| PDAC (TCGA + GTEx) | 16,598 | 0.753 | 1.9e-14 | 7 of 50 (31x) |
| LUAD (TCGA + GTEx) | 18,240 | 0.829 | 2.5e-22 | 6 of 50 (30x) |

## Quickstart

Install Python 3.10+, then:

    pip install -r requirements.txt

Run a cancer type end to end (fetches TCGA/GTEx/DepMap/Open Targets data,
computes all axes, ranks, validates):

    python -m multiaxisoncology.cli run configs/pdac.yaml

Validate a ranking against its benchmark:

    python -m multiaxisoncology.cli validate configs/pdac.yaml <ranking.csv>

Reproduce the published PDAC result exactly:

    python tests/test_regression_pdac.py

Launch the dashboard:

    streamlit run dashboard/app.py

Adding a new cancer type is a config file: copy `configs/luad.yaml`, change
the TCGA study, GTEx tissue, DepMap lineage, and benchmark CSV.

## Repository structure

    configs/          per-cancer settings (frozen weights, data sources)
    src/multiaxisoncology/   fetch, axis computation, integration, validation
    tests/            regression test (must reproduce published ranking)
    scripts/          one-off analysis and figure scripts
    dashboard/        Streamlit app (disease-agnostic; drop in a CSV + JSON)

## Data sources

TCGA (via UCSC Xena GDC hub), GTEx v8, DepMap 22Q2, Open Targets Platform,
gnomAD v4.0. All data are public; no patient-identifiable data are used.

## License

MIT — see LICENSE.
