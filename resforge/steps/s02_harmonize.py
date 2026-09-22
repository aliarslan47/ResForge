"""s02 — Harmonizasyon: her native çıktıyı hAMRonization ile ortak şemaya çevir, sonra birleştir.
Çıktı: runs/<ts>/harmonized/<genome_id>__<tool>.tsv ve harmonized/combined.tsv
Not: hAMRonization her araç için farklı zorunlu metadata ister; per-tool doğru bayrak seti verilir."""
from __future__ import annotations

import json
from pathlib import Path

from .. import util

# hAMRonization parser'ının zorunlu metadata bayrakları (parser'a göre değişir).
# amrfinderplus/rgi input_file_name ister; abricate #FILE sütunundan alır.
NEEDS_INPUT_NAME = {"amrfinderplus", "rgi"}


def _harmonize_one(runner, parser, native, out_tsv, gid, sw_version):
    args = ["hamronize", parser,
            "--analysis_software_version", sw_version or "unknown",
            "--reference_database_version", "unknown",
            "--format", "tsv", "--output", str(out_tsv)]
    if parser in NEEDS_INPUT_NAME:
        args += ["--input_file_name", gid]
    args += [str(native)]
    return runner.run(f"hamronize_{gid}_{parser}", args, conda_env=util.ENV["hamronize"],
                      version_cmd=["hamronize", "--version"], check=False)


def run(cfg: dict, run_dir: Path, runner_factory) -> dict:
    manifest = json.loads((run_dir / "scan_manifest.json").read_text(encoding="utf-8"))
    hdir = run_dir / "harmonized"
    hdir.mkdir(parents=True, exist_ok=True)
    runner = runner_factory(hdir)

    harmonized_files = []
    per = {"harmonized": [], "skipped": []}
    for rec in manifest["records"]:
        native = Path(rec["native_output"])
        gid, tool, parser = rec["genome_id"], rec["tool"], rec["hamronize_parser"]
        if rec["status"] != "OK" or not (native.exists() and native.stat().st_size > 0):
            per["skipped"].append({"genome_id": gid, "tool": tool, "reason": "native yok/boş/başarısız"})
            continue
        out_tsv = hdir / f"{gid}__{tool}.tsv"
        prov = _harmonize_one(runner, parser, native, out_tsv, gid, rec.get("version"))
        if prov.get("exit_code") == 0 and out_tsv.exists() and out_tsv.stat().st_size > 0:
            harmonized_files.append(out_tsv)
            per["harmonized"].append({"genome_id": gid, "tool": tool, "file": str(out_tsv)})
        else:
            per["skipped"].append({"genome_id": gid, "tool": tool, "reason": f"hamronize exit {prov.get('exit_code')}"})

    combined = hdir / "combined.tsv"
    if harmonized_files:
        runner.run("hamronize_summarize",
                   ["hamronize", "summarize", "--summary_type", "tsv", "--output", str(combined)]
                   + [str(p) for p in harmonized_files],
                   conda_env=util.ENV["hamronize"], check=False)
    else:
        print("[s02] WARNING: harmonize edilecek çıktı yok — combined.tsv üretilmedi.")

    per["combined"] = str(combined) if combined.exists() else None
    (run_dir / "harmonize_manifest.json").write_text(
        json.dumps(per, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"[s02] {len(per['harmonized'])} harmonize, {len(per['skipped'])} atlandı; "
          f"combined={'var' if per['combined'] else 'YOK'}")
    return per
