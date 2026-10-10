"""Shared helpers: Ensembl <-> symbol mapping, percentile scoring."""
import pandas as pd


def load_probemap(path):
    """GENCODE probemap -> DataFrame with versionless ensembl id + symbol."""
    pm = pd.read_csv(path, sep="\t")
    pm["ens"] = pm["id"].str.replace(r"\.\d+$", "", regex=True)
    return pm.drop_duplicates("ens")


def ens2sym_map(pm):
    return pm.drop_duplicates("ens").set_index("ens")["gene"]


def sym2ens_map(pm):
    return pm.drop_duplicates("gene").set_index("gene")["ens"]


def pct(s):
    return s.rank(pct=True)


def neutral(s):
    """Percentile with missing values -> neutral 0.5."""
    return pct(s).fillna(0.5)
