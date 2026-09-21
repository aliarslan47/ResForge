# ResForge — DURUM

**Proje:** AMR genotip→fenotip tahmin güvenilirliği. Özgün çekirdek = **çoklu-araç uyuşmazlığı** (hAMRonization). Bağımsız proje — AMRForge (Transfer Learning for AMR) ile ÇAKIŞMAZ; onun ML tahmin işini tekrar etmez.

**İlk tür:** *Acinetobacter baumannii*
**İşlem:** yalnız CPU, GPU yok, yeni wet-lab yok.

## Yapıldı
- Repo iskeleti: `pyproject.toml`, `environment.yml`, `.gitignore`, `config/config.yaml`, çift dilli README (EN+TR).
- Paket: `resforge/` — `config_loader` (${RESFORGE_*} genişletme), `runner` (izole conda + provenance + dürüst exit kaydı), `util` (tool→env haritası), `cli` (info/run/scan/harmonize/concord/phenotype/report).
- Adımlar: `steps/s01_scan` (4 tarayıcı: AMRFinderPlus/RGI/ABRicate-CARD/ABRicate-ResFinder) → `s02_harmonize` (hAMRonization) → `s03_concord` (Jaccard + gen/araç uyuşmazlığı) → `s04_phenotype` (ince AST merceği) → `s05_report` (çift dilli HTML).
- Hazır conda ortamları kullanılıyor (yeni kurulum yok): `ali-amrfinder`, `ali-rgi`, `ali-virulence` (abricate), `hamr` (hAMRonization, test edildi).

## Sırada
1. **Veri (Claude buluyor):** taze halka açık *A. baumannii* genom + AST seti (izolasyon: AMRForge `TL_AMR_data` KULLANILMAZ). `data/genomes/*.fasta` + `data/ast.tsv` (genome_id, antibiotic, phenotype) + `data/antibiotic_class.tsv`.
2. **Küçük pilot:** 20-50 genomda `resforge run` → süre + doğruluk kontrolü (Ali kuralı: denenmemiş kodu tam sette koşma).
3. Uyuşmazlık + fenotip merceği doğrulaması (ör. karbapenem geni var ama AST duyarlı → doğru "uyumsuz" mu).
4. PipelineForge ile `docs/pipeline_architecture.html` (spec: `PipelineForge/specs/resforge.yml`) + GitHub Pages.
5. Ölçekle, commit+push, çift dilli rapor cilası.

## Notlar
- Dürüstlük: araç yok/başarısız/boş = WARNING; uydurma sonuç yok; her araç çağrısı `<tool>.provenance.json`.
- hAMRonization AMRFinderPlus'ta v4 sütun şeması ister (test edildi); native çıktı sürüm uyumu pilotta doğrulanacak.
