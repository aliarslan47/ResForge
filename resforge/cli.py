"""ResForge CLI: info / run / scan / harmonize / concord / phenotype / report.
Kohort akışı: s01_scan → s02_harmonize → s03_concord → s04_phenotype → s05_report."""
from __future__ import annotations

import argparse
import datetime as dt
import os
from pathlib import Path

from . import __version__, util
from .config_loader import load_config
from .runner import ToolRunner
from .steps import s01_scan, s02_harmonize, s03_concord, s04_phenotype, s05_report

STEPS = {
    "scan": s01_scan.run,
    "harmonize": s02_harmonize.run,
    "concord": s03_concord.run,
    "phenotype": s04_phenotype.run,
    "report": s05_report.run,
}
ORDER = ["scan", "harmonize", "concord", "phenotype", "report"]


def _runner_factory(logs_dir):
    return ToolRunner(logs_dir)


def _new_run_dir(cfg) -> Path:
    ts = dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    label = cfg.get("project", {}).get("species", "cohort").split()[0].lower()
    d = Path(cfg["paths"]["work"]) / f"{ts}_{label}"
    d.mkdir(parents=True, exist_ok=True)
    return d


def _resolve_run_dir(cfg, arg) -> Path:
    if arg:
        return Path(arg)
    work = Path(cfg["paths"]["work"])
    existing = sorted([p for p in work.glob("*_*") if p.is_dir()]) if work.exists() else []
    if not existing:
        raise SystemExit("Mevcut run yok. Önce 'resforge run' veya 'resforge scan' çalıştırın (--run-dir verin).")
    return existing[-1]


def _cmd_info(cfg, args):
    print(f"ResForge {__version__}")
    print(f"HOME: {os.environ.get('RESFORGE_HOME')}")
    print(f"WORK: {cfg['paths']['work']}")
    print(f"DATA: {cfg['paths']['data']}")
    print(f"Tür : {cfg.get('project', {}).get('species')}")
    print(f"CPU threads (hesap): {util.threads(cfg)}")
    print(f"Tarayıcılar: {[k for k,v in cfg.get('scanners',{}).items() if v]}")
    print(f"conda env haritası: {util.ENV}")
    genomes = util.list_genomes(cfg)
    print(f"Genom (data/genomes): {len(genomes)} dosya"
          + ("" if genomes else f"  [YOK — {cfg['input']['genomes_glob']}]"))
    print(f"AST tablosu: {cfg['input']['ast_table']} "
          f"({'var' if Path(cfg['input']['ast_table']).exists() else 'YOK'})")


def _cmd_run(cfg, args):
    run_dir = _new_run_dir(cfg)
    print(f"[run] run_dir: {run_dir}")
    ctx = (cfg, run_dir, _runner_factory)
    for name in ORDER:
        print(f"--- {name} ---")
        STEPS[name](*ctx)
    print(f"[run] bitti → {run_dir}")


def _cmd_step(cfg, args):
    run_dir = _new_run_dir(cfg) if args.step == "scan" and not args.run_dir else _resolve_run_dir(cfg, args.run_dir)
    print(f"[{args.step}] run_dir: {run_dir}")
    STEPS[args.step](cfg, run_dir, _runner_factory)


def main(argv=None):
    ap = argparse.ArgumentParser(prog="resforge", description="ResForge — AMR araç uyuşmazlığı & genotip-fenotip güvenilirliği")
    ap.add_argument("--config", default=None, help="config.yaml yolu (varsayılan: config/config.yaml)")
    sub = ap.add_subparsers(dest="cmd", required=True)

    sub.add_parser("info", help="ortam / yollar / girdi durumu")
    sub.add_parser("run", help="tüm akışı çalıştır (scan→report)")
    for s in ORDER:
        sp = sub.add_parser(s, help=f"yalnız {s} adımı")
        sp.add_argument("--run-dir", default=None, help="mevcut run dizini (scan hariç zorunlu olabilir)")

    args = ap.parse_args(argv)
    cfg = load_config(args.config)

    if args.cmd == "info":
        _cmd_info(cfg, args)
    elif args.cmd == "run":
        _cmd_run(cfg, args)
    else:
        args.step = args.cmd
        _cmd_step(cfg, args)


if __name__ == "__main__":
    main()
