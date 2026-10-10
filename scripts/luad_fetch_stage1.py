"""LUAD fetch stage 1: slice expression matrices + GTEx lung tissue stats."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from multiaxisoncology.config import load_config
from multiaxisoncology.fetch import expression as fx
from multiaxisoncology.fetch import gtex_stats as fg
from multiaxisoncology.fetch import survival as fs

cfg = load_config(Path(__file__).parent.parent / "configs" / "luad.yaml")
Path(cfg["paths"]["work_dir"]).mkdir(parents=True, exist_ok=True)

print("=== slicing expression ===", flush=True)
fx.slice_expression(cfg)
print("=== gtex lung tissue stats ===", flush=True)
fg.gtex_tissue_stats(cfg)
print("=== survival + clinical ===", flush=True)
fs.fetch_survival_clinical(cfg)
print("STAGE1 DONE", flush=True)
