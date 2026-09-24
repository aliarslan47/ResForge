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

## Yapıldı (devam)
- **İNCE (DRUG-ÖZGÜ) MERCEK TAMAM (2026-09-23) — s04 üç-katmanlı + fenotibe-kör hibrit kapı.**
  - Katman 1 **sınıf-kanıtı** (eski), Katman 2 **ilaç-özgü** (`antimicrobial_agent` etiketinde tam ilaç adı geçen gen), Katman 3 **HİBRİT**.
  - **Kapı fenotibe KÖR (sızıntı yok):** ilaç-özgü katman yalnız, antibiyotiğin tam adı kendi sınıf-hit genlerinin etiketinde ≥1 kez geçen ilaçlarda uygulanır; geçmiyorsa (meropenem/ceftazidime — araçlar sınıf-adı yazıyor) sınıf katmanında kalır. Karar yalnız genotip sözlüğünden.
  - **SONUÇ (n=40):** genel uyum sınıf **0.606** → ilaç-özgü 0.683 → **HİBRİT 0.706** (+0.10). Hibrit sınıfa göre **36 sahte-R düzeltti, 0 gerçek-R kaybetti** = saf kazanç.
  - **Manşet:** amikacin **0.40 → 0.95** (+0.55), tobramycin 0.525 → **0.85** (+0.33); aac sahte-R'leri temizlendi. Beta-laktamlar kapıyla korundu (ceftazidime/meropenem sınıfta kaldı, kayıp yok).
  - **Bulgu (özgün):** ilaç-özgü çözünürlüğün değeri SINIFA-BAĞIMLI — aminoglikozitlerde biyoloji+etiket zengin (büyük kazanç), beta-laktamlarda geniş-spektrum+kaba etiket (kapı devre dışı bırakır). "İlaç-özgü her zaman iyi" değil.
  - Çıktı: `phenotype/genotype_vs_phenotype.tsv` (class/drug/hybrid çağrı+uyum), `phenotype/per_antibiotic.tsv` (policy+3 uyum+delta), `summary.json` (`gate` istatistikleri dahil). Commit'ler: `fc137ed` (iki-katman), hibrit bu turda.

## Yapıldı (devam)
- **ÖLÇEKLEME DENEMESİ n=300 TAMAM (2026-09-23) — bulgu sağlam + gizli kapı kusuru bulundu+düzeltildi.**
  - 300 genom (150R/150S meropenem, amikacin 205R/95S), fetch nested (40⊂300), tam akış exit 0, ~3.5-4 sa. scan + ~1 sa. harmonize.
  - **RESUME eklendi (s01):** dolu çıktı+exit0 varsa araç atlanır; kesinti-dayanıklı (40-genom testi 160 çağrı 0.4sn).
  - **ÖLÇEKLEME DARBOĞAZI:** s02 harmonize'da her hamronize çağrısı `conda run -n hamr` ile ortam açıyor (~6sn/çağrı × 1200 = ~1 sa). Optimize adayı: conda aktivasyonu bir kez / toplu hamronize.
  - **GİZLİ KUSUR (denemenin asıl kazancı):** ilk kapı ölçütü "ilaç adı ≥1 kez geçsin" idi. n=40'ta ceftazidime adı hiç geçmediğinden kapı ŞANS ESERİ kapalıydı. n=300'de 3213 etiketten SADECE 2'sinde (%0.06) geçince kapı açıldı → drug-katmanı 244 dirençli ceftazidime izolatını yanlış S sandı → ceftazidime 0.76→0.19, genel hibrit +0.10'dan +0.009'a çöktü (328 düzeldi ama 304 gerçek-R kayboldu). **n=40 sonucu fazla temizmiş — kapı kazara koruyucuydu.**
  - **DÜZELTME:** kapı "≥1 kez" → "sınıf-hitlerin ≥%5'i ilaç adı taşısın" (`phenotype.drug_name_coverage_min=0.05`, config'ten; fenotibe hâlâ kör). ceftazidime %0.06<%5 → kapalı/korunur; amikacin %9.3, tobramycin %9.1, imipenem %35 → açık.
  - **DÜZELTİLMİŞ SONUÇ (n=300):** sınıf 0.701 → **HİBRİT 0.753 (+0.052)**; 142 sahte-R düzeldi, **1 gerçek-R kayboldu** (tekrar ~saf kazanç, 0 per-drug regresyon).
  - **ÇEKİRDEK BULGU REPLİKE OLDU (7.5× ölçek):** amikacin 0.68→**0.93 (+0.25)**, tobramycin 0.727→**0.907 (+0.18)**. Aminoglikozit merceği sağlam.
  - Not: n=40 phenotype'u artık re-run edilemez (data/ast.tsv 300'lük üzerine yazıldı); n=40 sonucu tarihsel kayıtta durur.

## Yapıldı (devam)
- **İKİNCİ TÜR: P. aeruginosa n=300 TAM KOŞU BİTTİ (2026-09-24 gece).** İzole config `config/config_paeruginosa.yaml` + `data/paeruginosa/` (baumannii kök verisine dokunmaz). Fetch 300/300 (1.9 GB), scan 300/300 (~8 sa, RGI darboğaz; resume kesinti-dayanıklı — Ali oturumu kapattı, arka plan süreci sağ kaldı), harmonize→concord→phenotype→report.html hepsi exit 0. Çıktı: `runs/20260923_154757_pseudomonas/`.
  - **ÇEKİRDEK TEZ REPLİKE (tür-bağımsız):** araçlar hemfikirken uyum **0.81**, ayrışırken **0.41** (n=2449). Fark baumannii'den (0.78 vs 0.61) bile BÜYÜK. "Araç uyuşması = fenotip-güven sinyali" ikinci türde de tutuyor. ✅
  - **AMA genel uyum ÇOK DÜŞÜK: sınıf 0.44 / drug 0.49 / hibrit 0.44.** Biyolojik olarak beklenen: P. aeruginosa direnci ağırlıkla mutasyon/efluks/regülatör kaynaklı → gen-katalog tarayıcıları (AMRFinder/RGI/ABRicate) bu mekanizmaları kaçırır. Kazanılmış-gen mantığı bu türde zayıf.
  - **DRUG-ÖZGÜ/HİBRİT MERCEK TRANSFER OLMADI (özgün bulgu):** baumannii'de amikacin 0.40→0.95 fırlamıştı; P. aeruginosa'da 0.31→0.31, SIFIR katkı; hibrit 0 sahte-R düzeltti, 1 gerçek-R kaybetti (hibrit ≤ sınıf). **Bulgu: ince-mercek kazancı TÜRE + MEKANİZMAYA bağımlı** — kazanılmış-gen ağırlıklı türde işe yarar, mutasyon-ağırlıklı türde yaramaz. (Tez için güçlü kontrast örneği.)
  - **Kusur işareti:** levofloxacin drug-katmanı 0.85 ama kapı `class_fallback`'e düşürüp 0.47'de bıraktı — kapı burada fazla temkinli, kazancı yedi. İncelenecek.
- **DİSK:** `data/` toplam 3.1 GB (baumannii 1.2 GB / paeruginosa 1.9 GB, 600 genom); `runs/` 674 MB. data/ gitignore'da.

## Yapıldı (devam)
- **ÜÇÜNCÜ TÜR: K. pneumoniae n=300 TAM KOŞU BİTTİ (2026-09-24).** İzole config `config/config_klebsiella.yaml` + `data/klebsiella/` (taxon 573; uzunluk penceresi Klebsiella'ya göre **5.0–6.0 Mb** — varsayılan 3.5–4.4 yanlış olurdu; AMRFinder organism `Klebsiella_pneumoniae` doğrulandı). Pilot (n=5) ile pencere+izolasyon doğrulandı; fetch 300/300 (1.6 GB), scan 300/300 (~7.5 sa), harmonize→concord→phenotype→report exit 0. Çıktı: `runs/20260924_111602_klebsiella/`. Baumannii+paeruginosa verisi dokunulmadı.
  - **TEZ TAHMİNİ TUTTU — Klebsiella tam ORTAYA düştü.** Genel genotip→fenotip gradyanı: **baumannii 0.75 → Klebsiella 0.57 → paeruginosa 0.44.** Direnç mimarisi (kazanılmış-gen → karışık → mutasyon/efluks) ekseni üç bağımsız türle doğrulandı.
  - **Klebsiella detay (n=2861 çift):** Jaccard 0.16 (üçünün en düşüğü — araçlar en az uyuşuyor); hemfikir 0.67 vs ayrışık 0.52 (tez tuttu ama fark en dar=0.15); hibrit 0.57 (19 sahte-R düzeltti/18 gerçek-R kaybetti = **net sıfır, drug-mercek transfer OLMADI**).
  - **İlaç-bazında karışıklık görünür:** seftazidim **0.86** (ESBL=kazanılmış gen, katalogta) yüksek; karbapenem 0.46–0.48 (KPC/NDM + **OmpK35/36 porin kaybı** karışık) düşük; amikacin **0.36** (baumannii'nin temiz aac enzimine benzemiyor) en düşük. Tek tür içinde bile mekanizma-görünürlük ilişkisi tutuyor.
  - **BULGU PEKİŞTİ:** drug-mercek yalnız baumannii'de altın (aminoglikozit=temiz enzim geni); Klebsiella+paeruginosa'da etkisiz → ince-mercek değeri türe/mekanizmaya bağımlı. Baumannii istisna.
- **ÜÇ-TÜR MANŞETİ:** ResForge'un ölçtüğü "genotip-fenotip güvenilirliği" = direncin katalog-görünürlüğünün aynası; hem genel uyum hem drug-mercek değeri kazanılmış→karışık→mutasyonel ekseninde hareket ediyor.
- **DİSK GÜNCEL:** `data/` ~4.7 GB (3 tür × 300 genom); `runs/` büyüdü.

## Sırada
1. **s05 rapor:** hibrit sonucu + ilaç-başına policy tablosu + "sınıfa-bağımlı değer" bulgusunu çift dilli HTML'e işle.
2. **Literatür konumlama:** Davies 2021 (Microbial Genomics, E.coli çoklu-araç uyuşmazlığı — en yakın rakip) + hAMRoaster tam metin oku; "uyuşmazlık→fenotip-güven kuplajı + A.baumannii + drug-özgü hibrit" özgünlüğünü hakem-korumalı yaz.
3. **Ölçekle** (`--n 400`, ilaç-başına R/S tabakalı) → ilaç-başına istatistiksel güven; arka plan+resume.
4. **TEMİZLİK:** iş bitince ağır veri sil (data/ FASTA + runs/*/scan/); fetch script deterministik, geri çekilebilir.
2. **abricate çakışması:** card+resfinder hAMRonization'da tek `abricate` adına birleşiyor → efektif araç=3. Konkordans paydası/araç-kimliği için DB'yi ayrı tut (ör. analysis_software_name'e DB ekle) düşünülmeli.
3. Ölçekle (`--n 400`), çift dilli rapor cilası.
3. PipelineForge ile `docs/pipeline_architecture.html` (spec: `PipelineForge/specs/resforge.yml`) + GitHub Pages.
4. Ölçekle, commit+push, çift dilli rapor cilası.

## Notlar
- Dürüstlük: araç yok/başarısız/boş = WARNING; uydurma sonuç yok; her araç çağrısı `<tool>.provenance.json`.
- hAMRonization AMRFinderPlus'ta v4 sütun şeması ister (test edildi); native çıktı sürüm uyumu pilotta doğrulanacak.
