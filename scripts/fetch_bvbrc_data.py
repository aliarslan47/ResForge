#!/usr/bin/env python3
"""ResForge veri çekici — BV-BRC (eski PATRIC) A. baumannii genom + AST pilot kohortu.

DÜRÜSTLÜK / İZOLASYON:
  - Yalnız halka açık BV-BRC verisi. AMRForge `TL_AMR_data` KULLANILMAZ (ayrı proje/organizma).
  - Deterministik seçim (sabit sıralama + seed) → tekrar-üretilebilir.
  - Her indirilen FASTA doğrulanır (contig>0, boyut>0); başarısız = atlanır + WARNING (sessiz hata yok).

ÜRETİR:
  data/genomes/<genome_id>.fasta        (seçilen kohort)
  data/ast.tsv                          (genome_id, antibiotic, phenotype[R/S/I])
  data/antibiotic_class.tsv             (antibiotic, drug_class_keyword)
  data/cohort_metadata.tsv              (provenance: kalite, kaynak, PMID, MIC özeti)

KULLANIM:
  python scripts/fetch_bvbrc_data.py --n 40                 # pilot
  python scripts/fetch_bvbrc_data.py --n 400 --pool 3000    # ölçekleme
"""
from __future__ import annotations

import argparse
import json
import sys
import time
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path

API = "https://www.bv-brc.org/api"
TAXON = 470  # Acinetobacter baumannii

# Antibiyotik -> drug_class_keyword (araç `drug_class` alanına SUBSTRING eşleşir; pilotta doğrulanır).
# Yalnız BV-BRC'de bol ve klinik olarak A. baumannii için anlamlı olanlar.
ANTIBIOTIC_CLASS = {
    "imipenem": "CARBAPENEM",
    "meropenem": "CARBAPENEM",
    "ciprofloxacin": "QUINOLONE",
    "levofloxacin": "QUINOLONE",
    "amikacin": "AMINOGLYCOSIDE",
    "tobramycin": "AMINOGLYCOSIDE",
    "gentamicin": "AMINOGLYCOSIDE",
    "tetracycline": "TETRACYCLINE",
    "ceftazidime": "CEPHALOSPORIN",
    "trimethoprim/sulfamethoxazole": "SULFONAMIDE",
    "ampicillin/sulbactam": "BETA-LACTAM",
}
# Denge için "karar antibiyotiği" (kohort R/S dengesi bunlara göre kurulur).
BALANCE_ANTIBIOTICS = ("meropenem", "imipenem")
# Bir genomun kohorta girmesi için gereken en az farklı antibiyotik sayısı.
MIN_ANTIBIOTICS = 3


def _get(path: str, params: str, accept="application/json", tries=4, timeout=90):
    """BV-BRC RQL GET; basit retry. Hata = yüksek sesle."""
    url = f"{API}/{path}/?{params}&http_accept={urllib.parse.quote(accept)}"
    last = None
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={
                "Accept": accept,
                "User-Agent": "ResForge/0.1 (BV-BRC data fetch; +https://github.com/aliarslan47/ResForge)",
            })
            with urllib.request.urlopen(req, timeout=timeout) as r:
                data = r.read()
                if accept == "application/json":
                    return json.loads(data.decode("utf-8"))
                return data.decode("utf-8", "replace")
        except Exception as e:  # noqa: BLE001
            last = e
            time.sleep(2 * (i + 1))
    raise RuntimeError(f"BV-BRC isteği başarısız ({tries} deneme): {url}\n  son hata: {last}")


def fetch_genome_pool(pool_size: int, taxon: int = TAXON,
                      len_min: int = 3_500_000, len_max: int = 4_400_000) -> list[dict]:
    """Yüksek kaliteli aday genomlar: Complete önce, sonra iyi WGS. Deterministik (genome_id sıralı).
    len_min/len_max organizmaya göre değişir (A. baumannii ~4 Mb, P. aeruginosa ~6.3 Mb)."""
    fields = "genome_id,genome_name,genome_status,genome_length,contigs,genome_quality,assembly_accession,contig_n50"
    out: list[dict] = []
    seen: set[str] = set()
    tiers = [
        f"and(eq(taxon_id,{taxon}),eq(public,true),eq(genome_status,Complete),eq(genome_quality,Good),"
        f"gt(genome_length,{len_min}),lt(genome_length,{len_max}))",
        f"and(eq(taxon_id,{taxon}),eq(public,true),eq(genome_status,WGS),eq(genome_quality,Good),"
        f"gt(genome_length,{len_min}),lt(genome_length,{len_max}),lt(contigs,200))",
    ]
    for q in tiers:
        start = 0
        page = 5000
        while len(out) < pool_size:
            params = f"{q}&select({fields})&sort(+genome_id)&limit({page},{start})"
            rows = _get("genome", params)
            if not rows:
                break
            for r in rows:
                gid = str(r.get("genome_id", ""))
                if gid and gid not in seen:
                    seen.add(gid)
                    out.append(r)
            if len(rows) < page:
                break
            start += page
        if len(out) >= pool_size:
            break
    return out[:pool_size]


def fetch_ast(genome_ids: list[str]) -> dict[str, list[dict]]:
    """Seçili genomların AST kayıtları (yalnız curated antibiyotikler, R/S). in() ile batch."""
    ab_list = ",".join(sorted(ANTIBIOTIC_CLASS))
    fields = "genome_id,antibiotic,resistant_phenotype,measurement,measurement_unit,measurement_sign,testing_standard,pmid"
    per_genome: dict[str, list[dict]] = defaultdict(list)
    batch = 120
    for i in range(0, len(genome_ids), batch):
        chunk = genome_ids[i:i + batch]
        gid_in = ",".join(chunk)
        q = (f"and(in(genome_id,({gid_in})),in(antibiotic,({ab_list})),"
             f"or(eq(resistant_phenotype,Resistant),eq(resistant_phenotype,Susceptible)))")
        params = f"{q}&select({fields})&limit(25000)"
        rows = _get("genome_amr", params)
        for r in rows:
            per_genome[str(r.get("genome_id", ""))].append(r)
    return per_genome


def select_cohort(pool: list[dict], ast: dict[str, list[dict]], n: int) -> list[dict]:
    """AST genişliği yeterli + karar antibiyotiğinde R/S dengeli deterministik kohort."""
    cand = []
    for g in pool:
        gid = str(g["genome_id"])
        recs = ast.get(gid, [])
        abs_tested = {r["antibiotic"].lower() for r in recs}
        if len(abs_tested) < MIN_ANTIBIOTICS:
            continue
        bal_pheno = None
        for ab in BALANCE_ANTIBIOTICS:
            m = [r for r in recs if r["antibiotic"].lower() == ab]
            if m:
                bal_pheno = m[0]["resistant_phenotype"]
                break
        if bal_pheno is None:  # karar antibiyotiği yoksa pilotta dışla (mercek anlamsız)
            continue
        g = dict(g, _n_ab=len(abs_tested), _bal=bal_pheno)
        cand.append(g)
    # Deterministik: geniş AST önce, sonra genome_id
    cand.sort(key=lambda g: (-g["_n_ab"], str(g["genome_id"])))
    res = [g for g in cand if g["_bal"] == "Resistant"]
    sus = [g for g in cand if g["_bal"] == "Susceptible"]
    half = n // 2
    picked = res[:half] + sus[:n - half]
    if len(picked) < n:  # bir taraf yetmezse diğerinden doldur
        extra = (res[half:] + sus[n - half:])
        picked += extra[: n - len(picked)]
    picked.sort(key=lambda g: str(g["genome_id"]))
    return picked[:n]


def download_fasta(gid: str, dest: Path) -> tuple[bool, int, int]:
    """genome_sequence FASTA indir + doğrula. -> (ok, boyut, contig)."""
    fasta = _get("genome_sequence", f"eq(genome_id,{gid})&limit(25000)",
                 accept="application/dna+fasta")
    n_contig = fasta.count(">")
    size = len(fasta.encode("utf-8"))
    if n_contig == 0 or size == 0:
        return False, size, n_contig
    dest.write_text(fasta, encoding="utf-8")
    return True, size, n_contig


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=40, help="kohort genom sayısı")
    ap.add_argument("--pool", type=int, default=1500, help="aday genom havuzu boyutu")
    ap.add_argument("--home", default=None, help="RESFORGE_HOME (varsayılan: repo kökü)")
    ap.add_argument("--taxon", type=int, default=TAXON, help="BV-BRC taxon_id (470=A.baumannii, 287=P.aeruginosa)")
    ap.add_argument("--subdir", default="", help="data/ altında organizma alt-klasörü (izolasyon; ör. paeruginosa)")
    ap.add_argument("--len-min", type=int, default=3_500_000, help="min genom uzunluğu (bp)")
    ap.add_argument("--len-max", type=int, default=4_400_000, help="max genom uzunluğu (bp)")
    args = ap.parse_args()

    home = Path(args.home) if args.home else Path(__file__).resolve().parents[1]
    data = home / "data" / args.subdir if args.subdir else home / "data"
    genomes = data / "genomes"
    genomes.mkdir(parents=True, exist_ok=True)

    print(f"[1/4] Aday havuz çekiliyor (taxon={args.taxon}, uzunluk {args.len_min//10**6}-{args.len_max//10**6}Mb, hedef {args.pool}) …")
    pool = fetch_genome_pool(args.pool, taxon=args.taxon, len_min=args.len_min, len_max=args.len_max)
    print(f"      {len(pool)} yüksek kaliteli genom bulundu.")

    print(f"[2/4] AST çekiliyor ({len(pool)} genom, {len(ANTIBIOTIC_CLASS)} antibiyotik) …")
    ast = fetch_ast([str(g["genome_id"]) for g in pool])
    print(f"      {sum(len(v) for v in ast.values())} AST kaydı, {len(ast)} genomda.")

    print(f"[3/4] Kohort seçiliyor (n={args.n}, denge={BALANCE_ANTIBIOTICS[0]}) …")
    cohort = select_cohort(pool, ast, args.n)
    if not cohort:
        print("HATA: uygun kohort bulunamadı (AST genişliği/denge kriteri).", file=sys.stderr)
        return 2
    nR = sum(1 for g in cohort if g["_bal"] == "Resistant")
    print(f"      {len(cohort)} genom seçildi ({nR} R / {len(cohort)-nR} S — {BALANCE_ANTIBIOTICS[0]}).")

    print(f"[4/4] FASTA indiriliyor + AST/metadata yazılıyor …")
    ast_rows, meta_rows, ok = [], [], 0
    for g in cohort:
        gid = str(g["genome_id"])
        dest = genomes / f"{gid}.fasta"
        good, size, ncont = download_fasta(gid, dest)
        if not good:
            print(f"  [WARN] {gid}: FASTA boş/indirilemedi — atlandı.")
            continue
        ok += 1
        recs = ast[gid]
        pmap = {"Resistant": "R", "Susceptible": "S", "Intermediate": "I"}
        mic_bits = []
        for r in recs:
            ab = r["antibiotic"].lower()
            if ab not in ANTIBIOTIC_CLASS:
                continue
            ph = pmap.get(r.get("resistant_phenotype", ""), "")
            if ph not in ("R", "S", "I"):
                continue
            ast_rows.append((gid, ab, ph))
            mic = f'{r.get("measurement_sign","")}{r.get("measurement","")}{r.get("measurement_unit","")}'.strip()
            mic_bits.append(f"{ab}={ph}({mic})")
        meta_rows.append({
            "genome_id": gid, "genome_name": g.get("genome_name", ""),
            "genome_status": g.get("genome_status", ""), "genome_quality": g.get("genome_quality", ""),
            "genome_length": g.get("genome_length", ""), "contigs": g.get("contigs", ""),
            "assembly_accession": g.get("assembly_accession", ""),
            "fasta_bytes": size, "fasta_contigs": ncont,
            "n_antibiotics": g["_n_ab"], "source": "BV-BRC",
            "ast_summary": "; ".join(sorted(mic_bits)),
        })

    # Yaz: ast.tsv
    with (data / "ast.tsv").open("w", encoding="utf-8") as f:
        f.write("genome_id\tantibiotic\tphenotype\n")
        for gid, ab, ph in sorted(set(ast_rows)):
            f.write(f"{gid}\t{ab}\t{ph}\n")
    # antibiotic_class.tsv
    with (data / "antibiotic_class.tsv").open("w", encoding="utf-8") as f:
        f.write("antibiotic\tdrug_class_keyword\n")
        for ab, cls in sorted(ANTIBIOTIC_CLASS.items()):
            f.write(f"{ab}\t{cls}\n")
    # cohort_metadata.tsv (provenance)
    cols = ["genome_id", "genome_name", "genome_status", "genome_quality", "genome_length",
            "contigs", "assembly_accession", "fasta_bytes", "fasta_contigs", "n_antibiotics",
            "source", "ast_summary"]
    with (data / "cohort_metadata.tsv").open("w", encoding="utf-8") as f:
        f.write("\t".join(cols) + "\n")
        for m in sorted(meta_rows, key=lambda r: r["genome_id"]):
            f.write("\t".join(str(m[c]) for c in cols) + "\n")

    print(f"\nTAMAM: {ok}/{len(cohort)} genom indirildi.")
    print(f"  data/genomes/*.fasta   ({ok} dosya)")
    print(f"  data/ast.tsv           ({len(set(ast_rows))} çift)")
    print(f"  data/antibiotic_class.tsv, data/cohort_metadata.tsv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
