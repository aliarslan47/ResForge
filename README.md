# ResForge

A reproducible study and toolkit for **AMR genotype→phenotype prediction reliability** — how well resistance genes predict measured resistance, and how much that answer depends on which tool you ask.

![stage](https://img.shields.io/badge/stage-design-c07211)
![compute](https://img.shields.io/badge/compute-CPU%20only-2f8f5b)
![type](https://img.shields.io/badge/type-in%20silico-0d6b8f)

[Türkçe](README.tr.md) · **English**

> Early stage. This repository is at the design/planning step — nothing here is validated yet. Claims are marked as planned until a real dataset is run.

## What is it?

ResForge is the antimicrobial-resistance member of the Forge family — same conventions as BacForge (bacterial WGS) and VirusForge (virus/phage), but a separate, isolated project. Instead of a raw-reads→report pipeline, ResForge is a study codebase that quantifies **two layers of uncertainty** in genotype-based AMR prediction:

1. **Tool discordance** — the same genome gives different resistance-gene calls depending on the caller (ABRicate, AMRFinderPlus, RGI/CARD, ResFinder). Harmonized with [hAMRonization](https://github.com/pha4ge/hAMRonization) and measured, not assumed.
2. **Genotype–phenotype gap** — even when tools agree, gene presence does not always match the measured phenotype (AST). Modelled with classical ML and, crucially, analysed *where and why* prediction fails.

No new wet-lab experiments; public genome + AST data only; GPU not required.

## What it does

Given a set of bacterial genomes with paired antibiotic-susceptibility (AST) results:

- runs several AMR gene callers on each genome and **harmonizes** their outputs into one schema;
- computes **concordance metrics** across tools (per gene family, per drug class);
- builds a genotype feature matrix and trains per-antibiotic classifiers (Random Forest / XGBoost / Logistic Regression, CPU mode);
- performs **error analysis** on discordant isolates — which gene combinations fail to predict phenotype, and whether the failure tracks tool disagreement;
- gathers everything into one bilingual (TR+EN) HTML report.

Honest by design: a value that is missing is reported as `WARNING`, a step that does not apply as `NOT_APPLICABLE`; no PASS without a real tool exit and real output; full input→tool→database→command→output provenance.

## Installation

```bash
git clone https://github.com/aliarslan47/ResForge.git
cd ResForge

conda env create -f environment.yml
conda activate resforge
pip install -e .
```

> Genome + AST data are downloaded by the user (public sources such as BV-BRC/PATRIC or NCBI Pathogen Detection) and placed under `data/`; raw data is never committed.

## Usage

```bash
# planned CLI — not yet implemented
python3 -m resforge.cli info                 # detected resources / paths
python3 -m resforge.cli scan   --data <dir>  # run AMR callers on genomes
python3 -m resforge.cli concord              # tool-discordance metrics
python3 -m resforge.cli predict --antibiotic <name>   # genotype→phenotype ML
```

## Modules

Planned structure. R00–R02 prepare the harmonized genotype table; R03 measures tool discordance; R04–R06 handle the genotype→phenotype question; R07 reports.

| Code | Module | Purpose |
|:---:|---|---|
| R00 | Data | ingest genomes + AST, build manifest |
| R01 | Scan | run AMR callers (ABRicate, AMRFinderPlus, RGI/CARD, ResFinder) |
| R02 | Harmonize | normalize all outputs with hAMRonization |
| R03 | Concordance | tool-vs-tool discordance metrics |
| R04 | Features | genotype feature matrix from harmonized calls |
| R05 | Predict | per-antibiotic classifiers (RF / XGBoost / LogReg) |
| R06 | Gap analysis | where/why genotype fails to predict phenotype |
| R07 | Report | bilingual self-contained HTML |
