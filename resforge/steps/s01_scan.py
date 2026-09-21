"""s01 — Tarama: her genomu birden çok AMR tarayıcısıyla tara.
Çıktı: runs/<ts>/scan/<genome_id>/<tool>.tsv (+ <tool>.provenance.json) ve scan_manifest.json
KURAL (Ali): kurulu olmayan/başarısız araç sessizce yutulmaz → exit_code kaydedilir, manifest'e yazılır."""
from __future__ import annotations

import json
from pathlib import Path

from .. import util


def _amrfinder(runner, genome, out_dir, cfg, t):
    out = out_dir / "amrfinderplus.tsv"
    tp = cfg.get("tools", {}).get("amrfinderplus", {})
    cmd = ["amrfinder", "-n", str(genome), "--threads", str(t)]
    cmd += list(tp.get("extra", ["--plus"]))
    if tp.get("organism"):
        cmd += ["--organism", tp["organism"]]
    db = Path(cfg["paths"]["db"]) / "amrfinderplus" / "latest"
    if db.exists():
        cmd += ["--database", str(db)]   # yoksa araç kendi kurulu DB'sini kullanır
    return runner.run("amrfinderplus", cmd, conda_env=util.ENV["amrfinder"],
                      version_cmd=["amrfinder", "--version"], stdout_path=str(out), check=False), out


def _rgi(runner, genome, out_dir, t):
    prefix = out_dir / "rgi"
    prov = runner.run("rgi", ["rgi", "main", "-i", str(genome), "-o", str(prefix),
                              "-t", "contig", "-n", str(t), "--clean"],
                      conda_env=util.ENV["rgi"], version_cmd=["rgi", "main", "--version"], check=False)
    return prov, Path(str(prefix) + ".txt")


def _abricate(runner, genome, out_dir, db, cfg):
    out = out_dir / f"abricate_{db}.tsv"
    ab = cfg.get("tools", {}).get("abricate", {})
    cmd = ["abricate", "--db", db, "--nopath",
           "--minid", str(ab.get("min_identity", 80)),
           "--mincov", str(ab.get("min_coverage", 80)), str(genome)]
    return runner.run(f"abricate_{db}", cmd, conda_env=util.ENV["abricate"],
                      version_cmd=["abricate", "--version"], stdout_path=str(out), check=False), out


# hAMRonization alt komutu + hangi native dosyaya karşılık geldiği
TOOL_MAP = {
    "amrfinderplus": {"hamronize": "amrfinderplus"},
    "rgi":           {"hamronize": "rgi"},
    "abricate_card": {"hamronize": "abricate"},
    "abricate_resfinder": {"hamronize": "abricate"},
}


def run(cfg: dict, run_dir: Path, runner_factory) -> dict:
    genomes = util.list_genomes(cfg)
    scan_root = run_dir / "scan"
    scan_root.mkdir(parents=True, exist_ok=True)
    which = cfg.get("scanners", {})

    manifest = {"genomes": [], "scanners": which, "records": []}
    if not genomes:
        raise FileNotFoundError(
            f"Genom bulunamadı: {cfg['input']['genomes_glob']} — data/genomes/ altına *.fasta koyun.")

    for g in genomes:
        gid = util.genome_id(g)
        manifest["genomes"].append(gid)
        gdir = scan_root / gid
        gdir.mkdir(parents=True, exist_ok=True)
        runner = runner_factory(gdir)
        t = util.threads(cfg)

        results = {}
        if which.get("amrfinderplus", True):
            prov, out = _amrfinder(runner, g, gdir, cfg, t)
            results["amrfinderplus"] = (prov, out)
        if which.get("rgi", True):
            prov, out = _rgi(runner, g, gdir, t)
            results["rgi"] = (prov, out)
        if which.get("abricate_card", True):
            prov, out = _abricate(runner, g, gdir, "card", cfg)
            results["abricate_card"] = (prov, out)
        if which.get("abricate_resfinder", True):
            prov, out = _abricate(runner, g, gdir, "resfinder", cfg)
            results["abricate_resfinder"] = (prov, out)

        for tool, (prov, out) in results.items():
            ok = prov.get("exit_code") == 0
            has_out = out.exists() and out.stat().st_size > 0
            manifest["records"].append({
                "genome_id": gid, "tool": tool,
                "hamronize_parser": TOOL_MAP[tool]["hamronize"],
                "native_output": str(out), "exit_code": prov.get("exit_code"),
                "version": prov.get("version"), "status": "OK" if (ok and has_out) else "WARNING",
            })

    (run_dir / "scan_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
    warned = [r for r in manifest["records"] if r["status"] != "OK"]
    print(f"[s01] {len(genomes)} genom × {len([k for k,v in which.items() if v])} tarayıcı; "
          f"{len(warned)} WARNING (araç yok/başarısız/boş).")
    return manifest
