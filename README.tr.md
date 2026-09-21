# ResForge

**AMR genotip→fenotip tahmin güvenilirliği** üzerine tekrarlanabilir bir çalışma ve araç seti — direnç genleri ölçülen direnci ne kadar doğru tahmin ediyor ve bu cevap hangi aracı kullandığına ne kadar bağlı.

![aşama](https://img.shields.io/badge/a%C5%9Fama-tasar%C4%B1m-c07211)
![işlem](https://img.shields.io/badge/i%C5%9Flem-yaln%C4%B1z%20CPU-2f8f5b)
![tür](https://img.shields.io/badge/t%C3%BCr-in%20silico-0d6b8f)

**Türkçe** · [English](README.md)

> Erken aşama. Bu depo tasarım/planlama adımında — henüz hiçbir şey doğrulanmadı. Gerçek veri koşulana kadar iddialar "planlanan" olarak işaretlidir.

## Nedir?

ResForge, Forge ailesinin antibiyotik-direnci üyesi — BacForge (bakteri WGS) ve VirusForge (virüs/faj) ile aynı kurallar, ama ayrı ve izole bir proje. Ham veri→rapor pipeline'ı değil; genotip-tabanlı AMR tahminindeki **iki belirsizlik katmanını** ölçen bir çalışma kod tabanı:

1. **Araç uyuşmazlığı** — aynı genom, hangi araca sorduğuna göre farklı direnç-geni listesi veriyor (ABRicate, AMRFinderPlus, RGI/CARD, ResFinder). [hAMRonization](https://github.com/pha4ge/hAMRonization) ile hizalanıp varsayılmadan ölçülür.
2. **Genotip–fenotip açığı** — araçlar hemfikir olsa bile, gen varlığı ölçülen fenotiple (AST) her zaman örtüşmüyor. Klasik ML ile modellenir ve asıl önemlisi tahminin *nerede ve neden* başarısız olduğu incelenir.

Yeni laboratuvar deneyi yok; yalnızca halka açık genom + AST verisi; GPU gerekmez.

## Ne yapar?

Antibiyogram (AST) sonuçlarıyla eşleştirilmiş bir bakteri genom seti verildiğinde:

- her genomda birkaç AMR gen tarayıcısı çalıştırıp çıktılarını tek şemaya **harmonize eder**;
- araçlar arası **uyuşmazlık metriklerini** hesaplar (gen ailesi ve ilaç sınıfı bazında);
- genotip öznitelik matrisi kurup antibiyotik başına sınıflandırıcı eğitir (Random Forest / XGBoost / Lojistik Regresyon, CPU modu);
- uyuşmayan izolatlarda **hata analizi** yapar — hangi gen kombinasyonları fenotibi tahmin edemiyor ve bu başarısızlık araç uyuşmazlığıyla örtüşüyor mu;
- hepsini tek bir çift dilli (TR+EN) HTML raporda toplar.

Tasarımı gereği dürüst: eksik değer `WARNING`, uygulanamayan adım `NOT_APPLICABLE`; gerçek araç çıkışı ve gerçek çıktı olmadan PASS yok; tam girdi→araç→veritabanı→komut→çıktı izlenebilirliği.

## Kurulum

```bash
git clone https://github.com/aliarslan47/ResForge.git
cd ResForge

conda env create -f environment.yml
conda activate resforge
pip install -e .
```

> Genom + AST verisi kullanıcı tarafından indirilir (BV-BRC/PATRIC veya NCBI Pathogen Detection gibi halka açık kaynaklar) ve `data/` altına konur; ham veri asla depoya işlenmez.

## Kullanım

```bash
# planlanan CLI — henüz uygulanmadı
python3 -m resforge.cli info                 # algılanan kaynaklar / yollar
python3 -m resforge.cli scan   --data <dizin>  # genomlarda AMR tarayıcıları
python3 -m resforge.cli concord              # araç-uyuşmazlığı metrikleri
python3 -m resforge.cli predict --antibiotic <ad>   # genotip→fenotip ML
```

## Modüller

Planlanan yapı. R00–R02 harmonize genotip tablosunu hazırlar; R03 araç uyuşmazlığını ölçer; R04–R06 genotip→fenotip sorusunu ele alır; R07 raporlar.

| Kod | Modül | Amaç |
|:---:|---|---|
| R00 | Veri | genom + AST alımı, manifest oluşturma |
| R01 | Tarama | AMR tarayıcıları (ABRicate, AMRFinderPlus, RGI/CARD, ResFinder) |
| R02 | Harmonize | tüm çıktıları hAMRonization ile normalize etme |
| R03 | Uyuşmazlık | araç-araç uyuşmazlık metrikleri |
| R04 | Öznitelik | harmonize çağrılardan genotip öznitelik matrisi |
| R05 | Tahmin | antibiyotik başına sınıflandırıcı (RF / XGBoost / LogReg) |
| R06 | Açık analizi | genotibin fenotibi nerede/neden tahmin edemediği |
| R07 | Rapor | çift dilli, kendi kendine yeten HTML |
