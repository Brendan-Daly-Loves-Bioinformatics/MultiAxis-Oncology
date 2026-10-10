"""Config loading for MultiAxis Oncology."""
from pathlib import Path
import yaml


def load_config(path):
    path = Path(path)
    with open(path) as f:
        cfg = yaml.safe_load(f)
    cfg["_config_path"] = str(path)
    cfg["_name"] = path.stem
    # resolve relative paths against the config file's parent's parent (repo root)
    root = path.parent.parent
    for key in ("benchmark",):
        if cfg.get(key):
            p = Path(cfg[key])
            if not p.is_absolute():
                cfg[key] = str(root / p)
    for section in ("paths",):
        for k, v in cfg.get(section, {}).items():
            p = Path(v)
            if not p.is_absolute():
                cfg[section][k] = str(root / p)
    w = cfg["weights"]
    assert abs(sum(w.values()) - 1.0) < 1e-9, "weights must sum to 1"
    return cfg
