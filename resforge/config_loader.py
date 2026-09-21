"""Config yükleme: ${ENV} değişkenlerini çözer. KURAL: kodda mutlak yol YOK.
BacForge deseninden uyarlandı (ayrı proje, bağımsız kopya)."""
from __future__ import annotations

import os
import re
from pathlib import Path

try:
    import yaml
except ImportError:
    yaml = None

_ENV_PATTERN = re.compile(r"\$\{([A-Z_][A-Z0-9_]*)\}")


def ensure_env_defaults() -> str:
    """RESFORGE_HOME / _DB / _WORK ayarlanmadıysa makul varsayılan ver."""
    home = os.environ.get("RESFORGE_HOME") or str(Path(__file__).resolve().parents[1])
    os.environ.setdefault("RESFORGE_HOME", home)
    os.environ.setdefault("RESFORGE_DB", str(Path(home) / "databases"))
    os.environ.setdefault("RESFORGE_WORK", str(Path(home) / "runs"))
    return home


def _expand(value):
    if isinstance(value, str):
        prev, cur = None, value
        while cur != prev:
            prev = cur
            cur = _ENV_PATTERN.sub(lambda m: os.environ.get(m.group(1), m.group(0)), cur)
        return cur
    if isinstance(value, dict):
        return {k: _expand(v) for k, v in value.items()}
    if isinstance(value, list):
        return [_expand(v) for v in value]
    return value


def load_config(path: str | os.PathLike | None = None) -> dict:
    ensure_env_defaults()
    if path is None:
        path = Path(os.environ["RESFORGE_HOME"]) / "config" / "config.yaml"
    path = Path(path)
    if yaml is None:
        raise RuntimeError("pyyaml kurulu değil. 'conda env create -f environment.yml' çalıştır.")
    with open(path) as fh:
        raw = yaml.safe_load(fh) or {}
    return _expand(raw)
