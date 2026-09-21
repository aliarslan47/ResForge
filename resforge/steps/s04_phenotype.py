"""s04 — Fenotip merceği (İNCE): harmonize genotip ile AST fenotibini karşılaştır.
AMRForge'la çakışmaz — rakip sınıflandırıcı EĞİTİLMEZ. Yalnız basit örtüşme:
  "gen var → dirençli beklenir" kuralı AST ile ne kadar uyuyor, ve
  uyumsuzluklar araçların AYRIŞTIĞI yerlerde mi yoğunlaşıyor?

Girdi:
  - harmonized/combined.tsv (genotip)
  - data/ast.tsv           : genome_id, antibiotic, phenotype (R/S/I)
  - data/antibiotic_class.tsv : antibiotic, drug_class_keyword  (ilaç→sınıf eşlemesi; yoksa WARNING)
Çıktı: phenotype/genotype_vs_phenotype.tsv, phenotype/summary.json
"""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from .s03_concord import _norm_gid


def run(cfg: dict, run_dir: Path, runner_factory=None) -> dict:
    pdir = run_dir / "phenotype"
    pdir.mkdir(parents=True, exist_ok=True)
    combined = run_dir / "harmonized" / "combined.tsv"
    ast_path = Path(cfg["input"]["ast_table"])
    class_map = Path(cfg["paths"]["data"]) / "antibiotic_class.tsv"

    def _warn(reason):
        print(f"[s04] WARNING: {reason} — fenotip merceği atlandı.")
        (pdir / "summary.json").write_text(json.dumps({"status": "SKIPPED", "reason": reason},
                                                       indent=2, ensure_ascii=False), encoding="utf-8")
        return {"status": "SKIPPED", "reason": reason}

    if not (combined.exists() and combined.stat().st_size > 0):
        return _warn("combined.tsv yok")
    if not ast_path.exists():
        return _warn(f"AST tablosu yok: {ast_path}")
    if not class_map.exists():
        return _warn(f"antibiyotik→sınıf eşlemesi yok: {class_map}")

    df = pd.read_csv(combined, sep="\t", dtype=str).fillna("")
    manifest = json.loads((run_dir / "scan_manifest.json").read_text(encoding="utf-8"))
    known = set(manifest.get("genomes", []))
    df["genome_id"] = df["input_file_name"].map(lambda v: _norm_gid(v, known))
    df["drug_class"] = df["drug_class"].str.upper()
    df["gene"] = df["gene_symbol"].str.strip()
    df["tool"] = df["analysis_software_name"].str.strip()

    ast = pd.read_csv(ast_path, sep="\t", dtype=str).fillna("")
    ast.columns = [c.strip().lower() for c in ast.columns]
    cmap = pd.read_csv(class_map, sep="\t", dtype=str).fillna("")
    cmap.columns = [c.strip().lower() for c in cmap.columns]
    ab2class = {r["antibiotic"].strip().lower(): r["drug_class_keyword"].strip().upper()
                for _, r in cmap.iterrows()}

    n_tools = df["tool"].nunique()
    # gen bazında araç desteği (uyuşmazlık bağı)
    support = (df.groupby(["genome_id", "gene"])["tool"].nunique()).to_dict()

    rows = []
    for _, a in ast.iterrows():
        gid, ab, pheno = a.get("genome_id", "").strip(), a.get("antibiotic", "").strip(), a.get("phenotype", "").strip().upper()
        if pheno not in ("R", "S", "I") or gid == "" or ab == "":
            continue
        cls = ab2class.get(ab.lower())
        if not cls:
            continue
        hits = df[(df["genome_id"] == gid) & (df["drug_class"].str.contains(cls, na=False))]
        genotype_R = not hits.empty
        # bu antibiyotiği destekleyen genlerin araç desteği (düşük destek = uyuşmazlıklı)
        supports = [support.get((gid, gene), 0) for gene in hits["gene"].unique()]
        min_support = min(supports) if supports else 0
        max_support = max(supports) if supports else 0
        agree = (genotype_R and pheno == "R") or (not genotype_R and pheno == "S")
        rows.append({
            "genome_id": gid, "antibiotic": ab, "drug_class": cls, "phenotype": pheno,
            "genotype_call": "R" if genotype_R else "S", "agreement": agree,
            "n_support_genes": len(hits["gene"].unique()),
            "max_tool_support": max_support, "min_tool_support": min_support,
            "discordant_support": max_support < n_tools if genotype_R else "",
        })

    out = pd.DataFrame(rows)
    out.to_csv(pdir / "genotype_vs_phenotype.tsv", sep="\t", index=False)

    summary = {"status": "OK", "n_pairs": int(len(out))}
    if len(out):
        graded = out[out["phenotype"].isin(["R", "S"])]
        summary["overall_agreement"] = round(float(graded["agreement"].mean()), 4) if len(graded) else None
        # uyuşmazlık bağı: destek TAM (tüm araçlar) olan vs olmayan çağrılarda uyum
        rcall = out[out["genotype_call"] == "R"]
        full = rcall[rcall["discordant_support"] == False]
        disc = rcall[rcall["discordant_support"] == True]
        summary["agreement_when_tools_concordant"] = round(float(full["agreement"].mean()), 4) if len(full) else None
        summary["agreement_when_tools_discordant"] = round(float(disc["agreement"].mean()), 4) if len(disc) else None
        summary["n_R_calls_tool_concordant"] = int(len(full))
        summary["n_R_calls_tool_discordant"] = int(len(disc))
    (pdir / "summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"[s04] {summary.get('n_pairs',0)} genom-antibiyotik çifti; "
          f"genel uyum={summary.get('overall_agreement')}")
    return summary
