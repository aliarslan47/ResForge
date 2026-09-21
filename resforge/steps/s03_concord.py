"""s03 — Uyuşmazlık: harmonize edilmiş combined.tsv'den araç-araç uyuşmazlık metrikleri.
Çıktı: concord/per_genome.tsv, concord/gene_discordance.tsv, concord/tool_profile.tsv, concord/summary.json

Metrikler:
  - genom başına araçlar arası Jaccard (gen kümesi kesişim/birleşim)
  - gen bazında: kaç araç çağırdı / hangi araçlar (uyuşmazlık sıklığı)
  - araç profili: her aracın ne kadarı ortak, ne kadarı yalnız o araçta (fazla/eksik çağırma eğilimi)
Not: genler araç-adlandırmasıyla (hAMRonization gene_symbol) eşleştirilir; isim varyasyonu bilinen bir sınırdır."""
from __future__ import annotations

import json
from itertools import combinations
from pathlib import Path

import pandas as pd


def _norm_gid(value: str, known_ids: set[str]) -> str:
    """input_file_name'i genome_id'ye indir (abricate dosya adı verebilir)."""
    if value in known_ids:
        return value
    stem = Path(str(value)).name
    for ext in (".fasta", ".fa", ".fna", ".fasta.gz", ".fna.gz"):
        if stem.endswith(ext):
            stem = stem[: -len(ext)]
    return stem


def run(cfg: dict, run_dir: Path, runner_factory=None) -> dict:
    combined = run_dir / "harmonized" / "combined.tsv"
    cdir = run_dir / "concord"
    cdir.mkdir(parents=True, exist_ok=True)
    if not (combined.exists() and combined.stat().st_size > 0):
        print("[s03] WARNING: combined.tsv yok — uyuşmazlık hesaplanamadı.")
        (cdir / "summary.json").write_text(json.dumps({"status": "NO_INPUT"}, indent=2), encoding="utf-8")
        return {"status": "NO_INPUT"}

    df = pd.read_csv(combined, sep="\t", dtype=str).fillna("")
    manifest = json.loads((run_dir / "scan_manifest.json").read_text(encoding="utf-8"))
    known = set(manifest.get("genomes", []))
    df["genome_id"] = df["input_file_name"].map(lambda v: _norm_gid(v, known))
    df["tool"] = df["analysis_software_name"].str.strip()
    df["gene"] = df["gene_symbol"].str.strip()
    df = df[df["gene"] != ""]

    tools = sorted(df["tool"].unique())
    genomes = sorted(df["genome_id"].unique())

    # 1) Genom başına Jaccard (araç çiftleri ortalaması)
    per_genome_rows = []
    for gid in genomes:
        sub = df[df["genome_id"] == gid]
        sets = {t: set(sub[sub["tool"] == t]["gene"]) for t in tools}
        union = set().union(*sets.values()) if sets else set()
        inter = set.intersection(*[s for s in sets.values() if s]) if any(sets.values()) else set()
        pair_j = []
        for a, b in combinations([t for t in tools if sets[t]], 2):
            u = sets[a] | sets[b]
            pair_j.append(len(sets[a] & sets[b]) / len(u) if u else 1.0)
        per_genome_rows.append({
            "genome_id": gid,
            "n_tools_with_calls": sum(1 for s in sets.values() if s),
            "union_genes": len(union),
            "core_genes_all_tools": len(inter),
            "mean_pairwise_jaccard": round(sum(pair_j) / len(pair_j), 4) if pair_j else "",
            **{f"n_{t}": len(sets[t]) for t in tools},
        })
    per_genome = pd.DataFrame(per_genome_rows)
    per_genome.to_csv(cdir / "per_genome.tsv", sep="\t", index=False)

    # 2) Gen bazında uyuşmazlık: her (genom,gen) kaç/hangi araç
    gd = (df.groupby(["genome_id", "gene"])["tool"]
            .agg(lambda x: sorted(set(x))).reset_index())
    gd["n_tools"] = gd["tool"].map(len)
    gd["tools"] = gd["tool"].map(lambda x: ",".join(x))
    gd["concordant_all"] = gd["n_tools"] == len(tools)
    gd[["genome_id", "gene", "n_tools", "tools", "concordant_all"]].to_csv(
        cdir / "gene_discordance.tsv", sep="\t", index=False)

    # 3) Araç profili: toplam / yalnız-bu-araç çağrı sayısı
    tool_rows = []
    for t in tools:
        calls = df[df["tool"] == t][["genome_id", "gene"]].drop_duplicates()
        total = len(calls)
        unique = 0
        for _, row in calls.iterrows():
            others = df[(df["genome_id"] == row["genome_id"]) & (df["gene"] == row["gene"]) & (df["tool"] != t)]
            if others.empty:
                unique += 1
        tool_rows.append({"tool": t, "total_calls": total, "unique_to_tool": unique,
                          "shared_calls": total - unique,
                          "unique_fraction": round(unique / total, 4) if total else ""})
    pd.DataFrame(tool_rows).to_csv(cdir / "tool_profile.tsv", sep="\t", index=False)

    jvals = [r["mean_pairwise_jaccard"] for r in per_genome_rows if r["mean_pairwise_jaccard"] != ""]
    summary = {
        "status": "OK", "n_genomes": len(genomes), "tools": tools,
        "overall_mean_pairwise_jaccard": round(sum(jvals) / len(jvals), 4) if jvals else None,
        "total_gene_calls": int(len(df)),
        "fully_concordant_gene_instances": int(gd["concordant_all"].sum()),
        "discordant_gene_instances": int((~gd["concordant_all"]).sum()),
    }
    (cdir / "summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"[s03] {len(genomes)} genom, araçlar={tools}, "
          f"ortalama Jaccard={summary['overall_mean_pairwise_jaccard']}")
    return summary
