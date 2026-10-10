"""Fetch TCGA survival + clinical matrices from the GDC Xena hub via xenaPython
dataset_fetch (direct download URLs are 403-gated). dataset_fetch returns a flat
list of per-probe rows (row i = probes[i] across samples). Writes sample-x-field
CSVs in the transposed format the Cox script expects.
"""
import csv
import xenaPython as xena


def fetch_survival_clinical(cfg):
    hub = cfg["xena_hub"]
    work = cfg["paths"]["work_dir"]
    code = cfg["disease_code"]
    for ds, out in [(cfg["xena_survival_dataset"], "survival"),
                    (cfg["xena_clinical_dataset"], "clinical")]:
        samples = xena.dataset_samples(hub, ds, None)
        probes = xena.dataset_field(hub, ds)
        rows = xena.dataset_fetch(hub, ds, samples, probes)
        path = f"{work}/{code}_{out}.csv"
        with open(path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["sample"] + list(samples))
            for p, row in zip(probes, rows):
                w.writerow([p] + list(row))
        print(out, len(samples), "samples ->", path)
