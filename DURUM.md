# ResForge — DURUM

**Proje:** AMR genotip→fenotip tahmin güvenilirliği. Özgün çekirdek = **çoklu-araç uyuşmazlığı** (hAMRonization). Bağımsız proje — AMRForge (Transfer Learning for AMR) ile ÇAKIŞMAZ; onun ML tahmin işini tekrar etmez.

**İlk tür:** *Acinetobacter baumannii*
**İşlem:** yalnız CPU, GPU yok, yeni wet-lab yok.

## Yapıldı
- Repo iskeleti: `pyproject.toml`, `environment.yml`, `.gitignore`, `config/config.yaml`, çift dilli README (EN+TR).
- Paket: `resforge/` — `config_loader` (${RESFORGE_*} genişletme), `runner` (izole conda + provenance + dürüst exit kaydı), `util` (tool→env haritası), `cli` (info/run/scan/harmonize/concord/phenotype/report).
- Adımlar: `steps/s01_scan` (4 tarayıcı: AMRFinderPlus/RGI/ABRicate-CARD/ABRicate-ResFinder) → `s02_harmonize` (hAMRonization) → `s03_concord` (Jaccard + gen/araç uyuşmazlığı) → `s04_phenotype` (ince AST merceği) → `s05_report` (çift dilli HTML).
- Hazır conda ortamları kullanılıyor (yeni kurulum yok): `ali-amrfinder`, `ali-rgi`, `ali-virulence` (abricate), `hamr` (hAMRonization, test edildi).

## Yapıldı (devam)
- **Veri ÇEKİLDİ (2026-09-22):** `scripts/fetch_bvbrc_data.py` — BV-BRC (eski PATRIC) A. baumannii genom+AST çekici (deterministik, provenance'lı, FASTA doğrulamalı). Kaynakta 16.811 AST'li genom / 207.797 R-S kaydı mevcut. Pilot kohort: **40 genom** (Good-kalite WGS, <150 contig, 3.5-4.4 Mb), meropenem'de **20R/20S dengeli**, 360 genom-antibiyotik çifti, 11 antibiyotik. `data/genomes/*.fasta` (157 MB, hepsi doğrulandı) + `data/ast.tsv` + `data/antibiotic_class.tsv` + `data/cohort_metadata.tsv`. AMRForge verisinden BAĞIMSIZ. (data/ gitignore'da; script commit'li.) Ölçekleme: `--n 400 --pool 3000`.

## Yapıldı (devam)
- **PİLOT KOŞTU (2026-09-22): duman testi (3 genom) + tam 40 genom.**
  - Duman testinde blocking bug bulundu+düzeltildi: `s02_harmonize` `hamronize <parser>` yerine `<parser>` çağırıyordu (exit 127) → commit `73d2dc5`.
  - Araçlar hazır+taze DB: AMRFinderPlus 4.2.7 (DB 2026-05-15), RGI 6.0.8/CARD 4.0.1, ABRicate 1.4.0 (card+resfinder 2026-Aug), hAMRonization 1.3.1.
  - 40 genom tam akış exit 0; süre ~27 dk (scan ağırlıklı ~40s/genom; rgi en yavaş). 160/160 araç çağrısı exit 0, 0 WARNING.
  - **ÇEKİRDEK SONUÇ (n=40):** araçlar hemfikirken genotip-fenotip uyumu **0.78** (n=23), ayrışırken **0.61** (n=301) → ResForge tezi doğru yönde. Genel uyum 0.61. Araç ayrışması yüksek: ort. pairwise Jaccard 0.22; 3220 gen çağrısından yalnız 103 tam-konkordan, 2100 ayrışık.
  - Antibiyotik-bazlı uyum: cipro/gentamisin 0.85, imipenem 0.62, mero/tetra/tobra ~0.53-0.55, seftazidim 0.57, **amikacin 0.40** (en düşük).
  - Çıktı: `runs/20260922_110252_acinetobacter/` (report.html + concord/ + phenotype/).

## Sırada
1. **İnce (drug-özgü) mercek:** amikacin 0.40 düşüklüğü sınıf-düzeyi keyword'ün kabalığından — aac(6') genleri AMINOGLYCOSIDE sayılıp amikacin-R öngörüyor ama izolat çoğu kez amikacin-S. ResFinder'ın ilaç-adı etiketlerinden (ör. `AMIKACIN;...;TOBRAMYCIN`) drug-özgü çözünürlük ekle → s04 iki-katmanlı: sınıf-kanıtı vs ilaç-özgü-kanıt. (ResForge'un "ince fenotip merceği" amacı bu.)
2. **abricate çakışması:** card+resfinder hAMRonization'da tek `abricate` adına birleşiyor → efektif araç=3. Konkordans paydası/araç-kimliği için DB'yi ayrı tut (ör. analysis_software_name'e DB ekle) düşünülmeli.
3. Ölçekle (`--n 400`), çift dilli rapor cilası.
3. PipelineForge ile `docs/pipeline_architecture.html` (spec: `PipelineForge/specs/resforge.yml`) + GitHub Pages.
4. Ölçekle, commit+push, çift dilli rapor cilası.

## Notlar
- Dürüstlük: araç yok/başarısız/boş = WARNING; uydurma sonuç yok; her araç çağrısı `<tool>.provenance.json`.
- hAMRonization AMRFinderPlus'ta v4 sütun şeması ister (test edildi); native çıktı sürüm uyumu pilotta doğrulanacak.
