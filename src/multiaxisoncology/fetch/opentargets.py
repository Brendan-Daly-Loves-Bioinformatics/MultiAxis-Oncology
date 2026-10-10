"""Query Open Targets Platform (v4 schema) for target-level tractability + known drugs.

Disease-agnostic: tractability and known-drug counts are target properties, so the
genome-wide table is computed once and reused across disease configs.
"""
import json
import time
import urllib.request
import pandas as pd

API = "https://api.platform.opentargets.org/api/v4/graphql"

QUERY = """
query targets($ids: [String!]!) {
  targets(ensemblIds: $ids) {
    id
    approvedSymbol
    tractability {
      label
      modality
      value
    }
    drugAndClinicalCandidates {
      count
    }
  }
}
"""


def _batch(ids):
    payload = json.dumps({"query": QUERY, "variables": {"ids": ids}}).encode()
    req = urllib.request.Request(API, data=payload,
                                 headers={"Content-Type": "application/json"})
    for attempt in range(4):
        try:
            with urllib.request.urlopen(req, timeout=60) as r:
                return json.loads(r.read())
        except Exception:
            if attempt == 3:
                raise
            time.sleep(5 * (attempt + 1))


def query_opentargets(pairs, out_path):
    """pairs: list of (gene_symbol, ensembl_id). Saves genome-wide table."""
    results = []
    for i in range(0, len(pairs), 50):
        chunk = pairs[i:i + 50]
        data = _batch([e for _, e in chunk])
        tg = {t["id"]: t for t in (data.get("data", {}).get("targets") or [])}
        for g, e in chunk:
            t = tg.get(e)
            if t is None:
                results.append({"gene": g, "ensembl": e, "found": False,
                                "tract_small_molecule": False, "tract_antibody": False,
                                "known_drugs": 0})
                continue
            tr = t.get("tractability") or []
            sm = any(x["value"] for x in tr if x["modality"] == "SM")
            ab = any(x["value"] for x in tr if x["modality"] == "AB")
            kd = ((t.get("drugAndClinicalCandidates") or {}).get("count") or 0)
            results.append({"gene": g, "ensembl": e, "found": True,
                            "tract_small_molecule": sm, "tract_antibody": ab,
                            "known_drugs": kd})
        if (i // 50) % 4 == 0:
            print(f"batch {i // 50 + 1}/{(len(pairs) - 1) // 50 + 1} done")
        time.sleep(0.5)

    df = pd.DataFrame(results)
    df.to_csv(out_path, index=False)
    print("saved", len(df), "rows;", int(df.found.sum()), "matched")
    return df
