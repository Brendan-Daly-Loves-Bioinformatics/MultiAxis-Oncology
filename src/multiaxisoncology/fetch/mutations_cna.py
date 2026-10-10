"""Genome-wide per-gene somatic mutation frequency (GDC Xena hub, sparse API)
+ ASCAT2 gene-level CNA amp/del frequencies (Ensembl-keyed, symbol-mapped).

Mutation frequency = fraction of patients with >=1 non-silent (damaging-class)
mutation. CNA: amp = total copy number >= 4, del = <= 1.
"""
import warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import xenaPython as xena

from ..utils import load_probemap, ens2sym_map

NON_SILENT_EXCLUDE = {"synonymous_variant", "3_prime_UTR_variant", "5_prime_UTR_variant",
                      "intron_variant", "non_coding_transcript_exon_variant",
                      "splice_region_variant;synonymous_variant"}


def fetch_mutations(cfg, genes):
    hub = cfg["xena_hub"]
    ds = cfg["xena_mutation_dataset"]
    samples = xena.dataset_samples(hub, ds, None)
    patients = set(s[:12] for s in samples)  # TCGA-XX-XXXX
    print("samples:", len(samples), "| patients:", len(patients))

    recs = []
    CHUNK = 1500
    for i in range(0, len(genes), CHUNK):
        chunk = genes[i:i + CHUNK]
        out = xena.sparse_data(hub, ds, samples, chunk)
        rows = out["rows"]
        for j in range(len(rows["sampleID"])):
            glist = rows["genes"][j] or []
            eff = rows["effect"][j]
            sid = rows["sampleID"][j]
            for g in glist:
                recs.append((g, sid, eff))
        if (i // CHUNK) % 5 == 0:
            print(f"chunk {i // CHUNK + 1}/{(len(genes) - 1) // CHUNK + 1}, "
                  f"records so far: {len(recs)}")

    mdf = pd.DataFrame(recs, columns=["gene", "sample", "effect"])
    mdf["patient"] = mdf["sample"].str[:12]
    dam = mdf[~mdf.effect.isin(NON_SILENT_EXCLUDE)]
    mut_freq = dam.groupby("gene")["patient"].nunique() / len(patients)
    return mut_freq


def fetch_cna(cfg):
    hub = cfg["xena_hub"]
    ds = cfg["xena_cna_dataset"]
    data = cfg["paths"]["data_dir"]
    samples2 = xena.dataset_samples(hub, ds, None)
    fields2 = xena.dataset_field(hub, ds)
    codes, rows2 = xena.dataset_probe_values(hub, ds, samples2, fields2)
    cna = pd.DataFrame(np.array(rows2), index=fields2, columns=samples2).apply(
        pd.to_numeric, errors="coerce")
    amp = (cna >= 4).mean(axis=1)
    dele = (cna <= 1).mean(axis=1)
    out = pd.DataFrame({"amp_freq": amp, "del_freq": dele})
    out.index = out.index.str.split(":").str[0].str.replace(r"\.\d+$", "", regex=True)
    pm = load_probemap(f"{data}/gencode.v23.probemap.tsv")
    out.index = out.index.map(ens2sym_map(pm))
    out = out.groupby(level=0).max()
    return out


def run(cfg, genes):
    work = cfg["paths"]["work_dir"]
    mut_freq = fetch_mutations(cfg, genes)
    mut_freq.to_csv(f"{work}/mutation_frequency.csv")
    print("genes with >=1 damaging mutation:", len(mut_freq))
    cna_out = fetch_cna(cfg)
    cna_out.to_csv(f"{work}/cna_frequencies.csv")
    print("genes with CNA:", len(cna_out))
    return mut_freq, cna_out
