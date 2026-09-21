"""s05 — Çift dilli (TR+EN) kendi kendine yeten HTML rapor.
Girdi: concord/*, phenotype/* özetleri.  Çıktı: runs/<ts>/report.html
Tasarım: tek dosya, gömülü CSS, TR/EN toggle; aile alışkanlığı (BacForge raporu TR-only idi, burada çift dilli)."""
from __future__ import annotations

import html
import json
from pathlib import Path

import pandas as pd


def _load_json(p: Path) -> dict:
    return json.loads(p.read_text(encoding="utf-8")) if p.exists() else {}


def _table(p: Path, max_rows=25) -> str:
    if not (p.exists() and p.stat().st_size > 0):
        return "<p class='muted'>—</p>"
    df = pd.read_csv(p, sep="\t", dtype=str).fillna("")
    return df.head(max_rows).to_html(index=False, border=0, classes="tbl", escape=True)


def _kv(d: dict, keys) -> str:
    out = []
    for k, label_tr, label_en in keys:
        v = d.get(k, "—")
        out.append(f"<div class='kv'><span data-tr='1'>{label_tr}</span><span data-en='1' hidden>{label_en}</span>"
                   f"<b>{html.escape(str(v))}</b></div>")
    return "".join(out)


_CSS = """
:root{--bg:#fff;--fg:#1a2230;--muted:#6b7683;--card:#f6f8fa;--line:#e2e8f0;--accent:#0d6b8f}
:root:not([data-theme=light]) @media (prefers-color-scheme:dark){}
@media (prefers-color-scheme:dark){:root:not([data-theme=light]){--bg:#0f1620;--fg:#e6edf3;--muted:#9aa7b4;--card:#161f2b;--line:#243040;--accent:#4bb3d6}}
:root[data-theme=dark]{--bg:#0f1620;--fg:#e6edf3;--muted:#9aa7b4;--card:#161f2b;--line:#243040;--accent:#4bb3d6}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--fg);font:15px/1.5 system-ui,-apple-system,Segoe UI,Roboto,sans-serif}
.wrap{max-width:960px;margin:0 auto;padding:24px 16px}
h1{font-size:1.5rem;margin:.2em 0}h2{font-size:1.15rem;margin:1.4em 0 .5em;color:var(--accent)}
.sub{color:var(--muted)}
.pills{display:flex;gap:6px;margin:12px 0}
.pills button{border:1px solid var(--line);background:var(--card);color:var(--fg);padding:5px 12px;border-radius:999px;cursor:pointer}
.pills button[aria-pressed=true]{background:var(--accent);color:#fff;border-color:var(--accent)}
.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(180px,1fr));gap:10px;margin:10px 0}
.kv{background:var(--card);border:1px solid var(--line);border-radius:10px;padding:10px 12px;display:flex;flex-direction:column;gap:2px}
.kv span{font-size:.8rem;color:var(--muted)}.kv b{font-size:1.15rem}
table.tbl{border-collapse:collapse;width:100%;font-size:.85rem;margin:.4em 0}
table.tbl th,table.tbl td{border-bottom:1px solid var(--line);padding:5px 8px;text-align:left}
table.tbl th{color:var(--muted);font-weight:600}
.muted{color:var(--muted)}.note{background:var(--card);border:1px solid var(--line);border-radius:10px;padding:12px;font-size:.85rem}
footer{margin-top:2em;color:var(--muted);font-size:.8rem;border-top:1px solid var(--line);padding-top:12px}
"""

_JS = """
function setLang(l){document.querySelectorAll('[data-tr]').forEach(e=>e.hidden=(l!=='tr'));
document.querySelectorAll('[data-en]').forEach(e=>e.hidden=(l!=='en'));
document.getElementById('bt').setAttribute('aria-pressed',l==='tr');
document.getElementById('be').setAttribute('aria-pressed',l==='en');}
"""


def run(cfg: dict, run_dir: Path, runner_factory=None) -> dict:
    concord = _load_json(run_dir / "concord" / "summary.json")
    pheno = _load_json(run_dir / "phenotype" / "summary.json")
    species = cfg.get("project", {}).get("species", "—")

    concord_kv = _kv(concord, [
        ("n_genomes", "Genom sayısı", "Genomes"),
        ("overall_mean_pairwise_jaccard", "Ortalama araç uyumu (Jaccard)", "Mean tool concordance (Jaccard)"),
        ("discordant_gene_instances", "Uyuşmazlık gösteren gen örneği", "Discordant gene instances"),
        ("fully_concordant_gene_instances", "Tam uyumlu gen örneği", "Fully concordant gene instances"),
    ])
    pheno_kv = _kv(pheno, [
        ("overall_agreement", "Genotip-fenotip genel uyumu", "Genotype-phenotype agreement"),
        ("agreement_when_tools_concordant", "Araçlar uyumluyken uyum", "Agreement when tools concordant"),
        ("agreement_when_tools_discordant", "Araçlar uyuşmazken uyum", "Agreement when tools discordant"),
    ])

    body = f"""<!doctype html><html lang="tr"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>ResForge Report</title><style>{_CSS}</style></head><body><div class="wrap">
<div class="pills"><button id="bt" aria-pressed="true" onclick="setLang('tr')">TR</button>
<button id="be" aria-pressed="false" onclick="setLang('en')">EN</button></div>

<h1>ResForge</h1>
<p class="sub" data-tr="1">AMR araç uyuşmazlığı & genotip-fenotip güvenilirliği · Tür: {html.escape(species)}</p>
<p class="sub" data-en="1" hidden>AMR tool discordance & genotype-phenotype reliability · Species: {html.escape(species)}</p>

<h2 data-tr="1">1. Araç Uyuşmazlığı</h2><h2 data-en="1" hidden>1. Tool Discordance</h2>
<div class="grid">{concord_kv}</div>
<h3 data-tr="1" class="sub">Araç profili</h3><h3 data-en="1" hidden class="sub">Tool profile</h3>
{_table(run_dir / 'concord' / 'tool_profile.tsv')}
<h3 data-tr="1" class="sub">Genom başına uyum</h3><h3 data-en="1" hidden class="sub">Per-genome concordance</h3>
{_table(run_dir / 'concord' / 'per_genome.tsv')}

<h2 data-tr="1">2. Fenotip Merceği</h2><h2 data-en="1" hidden>2. Phenotype Lens</h2>
<div class="grid">{pheno_kv}</div>
<div class="note"><span data-tr="1">Not: Bu bölüm rakip bir tahmin modeli değildir; yalnız "gen var → dirençli beklenir" kuralının AST ile örtüşmesini ve uyumsuzlukların araçların ayrıştığı yerlerde yoğunlaşıp yoğunlaşmadığını ölçer.</span>
<span data-en="1" hidden>Note: this is not a competing predictor; it only measures how the "gene present → expected resistant" rule matches AST, and whether mismatches concentrate where tools disagree.</span></div>

<footer>ResForge · <span data-tr="1">bağımsız in-silico çalışma · uydurma sonuç yok (eksik=WARNING)</span>
<span data-en="1" hidden>independent in-silico study · no fabricated results (missing=WARNING)</span></footer>
</div><script>{_JS}</script></body></html>"""

    out = run_dir / "report.html"
    out.write_text(body, encoding="utf-8")
    print(f"[s05] rapor yazıldı: {out}")
    return {"status": "OK", "report": str(out)}
