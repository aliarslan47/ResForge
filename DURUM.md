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

## Sırada
1. **Küçük pilot:** 40 genomda `resforge run` → süre + doğruluk kontrolü (Ali kuralı: denenmemiş kodu tam sette koşma). Önce 2-3 genomda duman testi.
2. Uyuşmazlık + fenotip merceği doğrulaması (ör. karbapenem geni var ama AST duyarlı → doğru "uyumsuz" mu). `antibiotic_class.tsv` keyword'leri gerçek `combined.tsv` drug_class sözlüğüne göre doğrula.
3. PipelineForge ile `docs/pipeline_architecture.html` (spec: `PipelineForge/specs/resforge.yml`) + GitHub Pages.
4. Ölçekle, commit+push, çift dilli rapor cilası.

## Notlar
- Dürüstlük: araç yok/başarısız/boş = WARNING; uydurma sonuç yok; her araç çağrısı `<tool>.provenance.json`.
- hAMRonization AMRFinderPlus'ta v4 sütun şeması ister (test edildi); native çıktı sürüm uyumu pilotta doğrulanacak.
