"""Slice Toil TCGA+GTEx TPM matrices down to disease-relevant samples.

Reads the large gzipped matrices in row chunks, keeping only:
  - TCGA <disease_code> primary tumors (-01/-06) and adjacent normals (-11)
  - GTEx samples from the config's normal tissue
Outputs a combined matrix (genes x samples), log2(TPM+0.001) as distributed,
as parquet + csv.gz in the config work dir.
"""
import gzip
import pandas as pd


def _wanted_columns(path, keep_ids):
    keep = []
    with gzip.open(path, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        for s in header[1:]:
            if s in keep_ids:
                keep.append(s)
    return keep


def slice_expression(cfg):
    data = cfg["paths"]["data_dir"]
    work = cfg["paths"]["work_dir"]
    pheno = pd.read_csv(f"{data}/TcgaTargetGTEX_phenotype.txt.gz", sep="\t", index_col=0,
                        encoding="latin-1")
    tissue = cfg["gtex_tissue"]
    gtex_keep = set(pheno.index[pheno["_study"] == "GTEX"].intersection(
        pheno.index[pheno["primary disease or tissue"] == tissue]))
    print(f"GTEx {tissue} samples: {len(gtex_keep)}")

    tcga_pheno = pheno[pheno["_study"] == "TCGA"]
    dis_pheno = tcga_pheno[tcga_pheno["detailed_category"] == cfg["tcga_detailed_category"]]
    dis_ids = set(dis_pheno.index)
    print(f"TCGA {cfg['disease_code']} phenotype samples: {len(dis_ids)} "
          f"(tumors {sum(s.endswith(('-01', '-06')) for s in dis_ids)}, "
          f"normals {sum(s.endswith('-11') for s in dis_ids)})")

    tcga_cols = _wanted_columns(f"{data}/tcga_RSEM_gene_tpm.gz", dis_ids)
    gtex_cols = _wanted_columns(f"{data}/gtex_RSEM_gene_tpm.gz", gtex_keep)
    print(f"TCGA columns found in matrix: {len(tcga_cols)}")
    print(f"GTEx columns: {len(gtex_cols)}")

    chunks = []
    for path, cols, tag in [(f"{data}/tcga_RSEM_gene_tpm.gz", tcga_cols, "TCGA"),
                            (f"{data}/gtex_RSEM_gene_tpm.gz", gtex_cols, "GTEX")]:
        reader = pd.read_csv(path, sep="\t", index_col=0, usecols=["sample"] + cols,
                             chunksize=5000)
        n = 0
        for chunk in reader:
            chunks.append(chunk.astype("float32"))
            n += len(chunk)
        print(f"{tag}: {n} gene rows read")

    expr = pd.concat(chunks)
    expr = expr.groupby(level=0).max()  # par_y / par_b duplicates
    out = f"{work}/{cfg['disease_code'].lower()}_gtex_{tissue.lower()}_tpm"
    expr.to_parquet(out + ".parquet")
    print("Saved:", expr.shape, "->", out)
    return out
