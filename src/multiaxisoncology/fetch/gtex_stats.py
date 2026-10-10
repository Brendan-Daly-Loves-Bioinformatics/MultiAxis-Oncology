"""Per-gene GTEx tissue expression stats: target-tissue mean TPM + max mean across
other tissues. Reads the full GTEx matrix in row chunks; accumulates per-tissue sums.
Disease-agnostic except for the target tissue column name.
"""
import numpy as np
import pandas as pd


def gtex_tissue_stats(cfg, target_tissue=None):
    data = cfg["paths"]["data_dir"]
    target_tissue = target_tissue or cfg["gtex_tissue"]
    pheno = pd.read_csv(f"{data}/TcgaTargetGTEX_phenotype.txt.gz", sep="\t", index_col=0,
                        encoding="latin-1")
    tissue = pheno["_primary_site"].to_dict()

    parts = []
    reader = pd.read_csv(f"{data}/gtex_RSEM_gene_tpm.gz", sep="\t", index_col=0,
                         chunksize=4000)
    for n_chunks, chunk in enumerate(reader, 1):
        cols_t = [tissue.get(s) for s in chunk.columns]
        cols_t = [t if isinstance(t, str) else "Unknown" for t in cols_t]
        chunk = chunk.clip(lower=-10)  # Toil floor is log2(TPM+0.001) ~ -9.97
        tissues = sorted(set(cols_t))
        acc = {t: np.zeros(chunk.shape[0]) for t in tissues}
        cnt = {t: 0 for t in tissues}
        for t in tissues:
            mask = np.array([c == t for c in cols_t])
            if mask.any():
                acc[t] += chunk.loc[:, np.array(chunk.columns)[mask]].sum(axis=1).to_numpy()
                cnt[t] = int(mask.sum())
        parts.append(pd.DataFrame({t: acc[t] / max(cnt[t], 1) for t in tissues},
                                  index=chunk.index))
        if n_chunks % 4 == 0:
            print(f"chunk {n_chunks}")

    ml_all = pd.concat(parts)
    ml = ml_all.groupby(level=0).max()  # par_y/par_b duplicates
    target = ml[target_tissue]
    other = ml.drop(columns=[target_tissue])
    out = pd.DataFrame({
        f"{target_tissue.lower()}_log2tpm": target,
        "max_other_tissue_log2tpm": other.max(axis=1),
        "max_other_tissue": other.idxmax(axis=1),
    })
    out_path = f"{data}/gtex_tissue_expression_stats_{target_tissue.lower()}.csv"
    out.to_csv(out_path)
    print("saved", out.shape, "->", out_path)
    return out_path
