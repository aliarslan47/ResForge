"""s04 — Fenotip merceği (İNCE): harmonize genotip ile AST fenotibini karşılaştır.
AMRForge'la çakışmaz — rakip sınıflandırıcı EĞİTİLMEZ. İki katmanlı basit örtüşme:

  Katman 1 (SINIF-kanıtı): antibiyotiğin sınıfından (ör. amikacin→AMINOGLYCOSIDE)
      HERHANGİ bir gen varsa "R beklenir". Kaba; sınıf-içi ilaç ayrımını gözden kaçırır
      (ör. aac(3)-Ia gentamisin verir ama amikacin vermez → sahte amikacin-R).

  Katman 2 (İLAÇ-özgü): sınıf-hitleri arasında, aracın `antimicrobial_agent`
      etiketinde TAM O antibiyotiği (ör. 'amikacin') anan bir gen varsa "R beklenir".
      Bu, sınıf-kabalığını çözer. Etiketi AMRFinderPlus + RGI taşır; ABRicate BOŞ
      (efektif 2 araç — bilinen kısıt, rapora yazılır).

Çıktı iki katmanın uyumunu yan yana + ilaç-başına kırılım verir → drug-özgü katmanın
sahte-R'leri (özellikle amikacin) düzelttiğini gösterir.

Girdi:
  - harmonized/combined.tsv (genotip; sütunlar: drug_class, antimicrobial_agent, ...)
  - data/ast.tsv              : genome_id, antibiotic, phenotype (R/S/I)
  - data/antibiotic_class.tsv : antibiotic, drug_class_keyword
Çıktı: phenotype/genotype_vs_phenotype.tsv, phenotype/per_antibiotic.tsv, phenotype/summary.json
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import pandas as pd

from .s03_concord import _norm_gid


def _agent_components(antibiotic: str) -> list[str]:
    """AST antibiyotik adını eşleşme bileşenlerine indir (bileşik ilaçları böl).
    'trimethoprim/sulfamethoxazole' -> ['trimethoprim','sulfamethoxazole']."""
    ab = antibiotic.strip().lower()
    parts = re.split(r"[\/\+;,]", ab)
    return [p.strip() for p in parts if p.strip()]


def _agent_names(agent_field: str) -> str:
    """`antimicrobial_agent` alanını normalize et (AMRFinder 'A/B', RGI 'a; b')."""
    return re.sub(r"[\/;,]", " ", str(agent_field)).lower()


def _drug_specific_hit(hits: pd.DataFrame, antibiotic: str) -> tuple[bool, bool]:
    """(ilaç-özgü R kanıtı var mı, herhangi bir araç ilaç-etiketi sağladı mı).
    İkinci bayrak: sınıf-hitlerinin en az birinde dolu antimicrobial_agent varsa True
    (yoksa ilaç-özgü çözünürlük UYGULANAMAZ → çağrıyı sınıf katmanına bırak)."""
    comps = _agent_components(antibiotic)
    any_annotation = False
    specific = False
    for agent in hits["antimicrobial_agent"].tolist():
        norm = _agent_names(agent)
        if norm.strip():
            any_annotation = True
            if any(c in norm for c in comps):
                specific = True
    return specific, any_annotation


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
    if "antimicrobial_agent" not in df.columns:
        df["antimicrobial_agent"] = ""

    ast = pd.read_csv(ast_path, sep="\t", dtype=str).fillna("")
    ast.columns = [c.strip().lower() for c in ast.columns]
    cmap = pd.read_csv(class_map, sep="\t", dtype=str).fillna("")
    cmap.columns = [c.strip().lower() for c in cmap.columns]
    ab2class = {r["antibiotic"].strip().lower(): r["drug_class_keyword"].strip().upper()
                for _, r in cmap.iterrows()}

    n_tools = df["tool"].nunique()
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

        # Katman 1: sınıf-kanıtı
        genotype_R_class = not hits.empty

        # Katman 2: ilaç-özgü
        specific, any_annotation = _drug_specific_hit(hits, ab) if not hits.empty else (False, False)
        if not genotype_R_class:
            genotype_R_drug = False
            drug_resolution = "no_class_hit"          # sınıf hiti yok → ikisi de S
        elif not any_annotation:
            genotype_R_drug = genotype_R_class         # etiket yok → sınıfa geri düş
            drug_resolution = "unresolved_falls_back"  # ör. yalnız ABRicate çağırmış
        else:
            genotype_R_drug = specific
            drug_resolution = "drug_specific"

        supports = [support.get((gid, gene), 0) for gene in hits["gene"].unique()]
        max_support = max(supports) if supports else 0
        agree_class = (genotype_R_class and pheno == "R") or (not genotype_R_class and pheno == "S")
        agree_drug = (genotype_R_drug and pheno == "R") or (not genotype_R_drug and pheno == "S")
        rows.append({
            "genome_id": gid, "antibiotic": ab, "drug_class": cls, "phenotype": pheno,
            "genotype_call_class": "R" if genotype_R_class else "S",
            "genotype_call_drug": "R" if genotype_R_drug else "S",
            "agreement_class": agree_class, "agreement_drug": agree_drug,
            "drug_resolution": drug_resolution,
            "n_support_genes": len(hits["gene"].unique()),
            "max_tool_support": max_support,
            "discordant_support": (max_support < n_tools) if genotype_R_class else "",
        })

    out = pd.DataFrame(rows)
    out.to_csv(pdir / "genotype_vs_phenotype.tsv", sep="\t", index=False)

    summary = {"status": "OK", "n_pairs": int(len(out)), "n_tools": int(n_tools)}
    per_ab_rows = []
    if len(out):
        graded = out[out["phenotype"].isin(["R", "S"])].copy()
        summary["overall_agreement_class"] = round(float(graded["agreement_class"].mean()), 4)
        summary["overall_agreement_drug"] = round(float(graded["agreement_drug"].mean()), 4)
        summary["n_flipped_to_S"] = int(((graded["genotype_call_class"] == "R") &
                                         (graded["genotype_call_drug"] == "S")).sum())
        # kaç sahte-R düzeldi (sınıf R & pheno S iken drug S'e döndü)
        fixed = graded[(graded["genotype_call_class"] == "R") & (graded["genotype_call_drug"] == "S") &
                       (graded["phenotype"] == "S")]
        broke = graded[(graded["genotype_call_class"] == "R") & (graded["genotype_call_drug"] == "S") &
                       (graded["phenotype"] == "R")]
        summary["false_R_fixed_by_drug_layer"] = int(len(fixed))
        summary["true_R_lost_by_drug_layer"] = int(len(broke))

        # ilaç-başına kırılım (amikacin iyileşmesini görünür kılar)
        for ab, g in graded.groupby("antibiotic"):
            per_ab_rows.append({
                "antibiotic": ab, "n": int(len(g)),
                "R": int((g["phenotype"] == "R").sum()), "S": int((g["phenotype"] == "S").sum()),
                "agreement_class": round(float(g["agreement_class"].mean()), 4),
                "agreement_drug": round(float(g["agreement_drug"].mean()), 4),
                "delta": round(float(g["agreement_drug"].mean() - g["agreement_class"].mean()), 4),
            })
        # araç-uyuşmazlık bağı (sınıf katmanı üzerinden, geriye-uyumlu)
        rcall = out[out["genotype_call_class"] == "R"]
        full = rcall[rcall["discordant_support"] == False]
        disc = rcall[rcall["discordant_support"] == True]
        summary["agreement_when_tools_concordant"] = round(float(full["agreement_class"].mean()), 4) if len(full) else None
        summary["agreement_when_tools_discordant"] = round(float(disc["agreement_class"].mean()), 4) if len(disc) else None

    per_ab = pd.DataFrame(per_ab_rows).sort_values("delta", ascending=False) if per_ab_rows else pd.DataFrame()
    per_ab.to_csv(pdir / "per_antibiotic.tsv", sep="\t", index=False)
    summary["per_antibiotic"] = per_ab_rows
    (pdir / "summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"[s04] {summary.get('n_pairs',0)} çift; sınıf-uyum={summary.get('overall_agreement_class')} "
          f"→ ilaç-özgü={summary.get('overall_agreement_drug')} "
          f"(sahte-R düzelen={summary.get('false_R_fixed_by_drug_layer')}, "
          f"gerçek-R kaybolan={summary.get('true_R_lost_by_drug_layer')})")
    return summary
