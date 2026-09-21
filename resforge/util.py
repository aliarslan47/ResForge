"""Ortak yardımcılar: tool→conda-env haritası, thread hesabı, kohort IO.
KURAL: env adları sistemde kurulu izole ortamlarla birebir eşleşir (yeni kurulum gerekmez)."""
from __future__ import annotations

import os
from glob import glob
from pathlib import Path

# Mantıksal araç anahtarı -> kurulu izole conda ortamı adı.
# (Keşifle doğrulandı: bunlar sistemde mevcut.)
ENV = {
    "amrfinder": "ali-amrfinder",   # NCBI AMRFinderPlus
    "rgi": "ali-rgi",               # CARD/RGI
    "abricate": "ali-virulence",    # ABRicate (card + resfinder DB'leri kurulu)
    "hamronize": "hamr",            # hAMRonization (test edildi)
}


def threads(cfg: dict) -> int:
    total = os.cpu_count() or 4
    res = (cfg or {}).get("resources", {})
    reserve = int(res.get("reserve_cores", 1))
    frac = float(res.get("thread_fraction", 0.8))
    return max(1, int((total - reserve) * frac))


def genome_id(path: str | Path) -> str:
    """Dosya adından kararlı genom kimliği (uzantısız stem)."""
    return Path(path).stem


def list_genomes(cfg: dict) -> list[Path]:
    pattern = cfg["input"]["genomes_glob"]
    return sorted(Path(p) for p in glob(pattern))
