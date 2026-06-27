from __future__ import annotations

import argparse
from pathlib import Path

from probioscore import __version__
from probioscore.engine.consolidated_engine import score_input


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="probioscore",
        description="ProbioScore portable AHP plus FCE probiogenomic calculator.",
    )
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--v4-tsv-dir", default=None, help="Directory containing V4 TSV outputs with the expected V4 subfolders.")
    source.add_argument("--input-tsv", default=None, help="Prebuilt ProbioScore or V4-derived TSV containing genome_id and ProbioScore_Status.")
    parser.add_argument("--out-dir", required=True, help="Output directory.")
    parser.add_argument("--mode", choices=["prospective_frozen", "benchmark_curated"], default="prospective_frozen", help="Decision mode. Default is prospective_frozen.")
    parser.add_argument("--config", default=None, help="Optional base AHP plus FCE JSON config. Embedded config is used by default.")
    parser.add_argument("--curated-tables-dir", default=None, help="Optional directory containing the five frozen curated TSV tables.")
    parser.add_argument("--allow-missing-curated-tables", action="store_true", help="Allow degraded execution if curated tables are missing.")
    parser.add_argument("--skip-visuals", action="store_true", help="Skip PNG chart generation.")
    parser.add_argument("--skip-xlsx", action="store_true", help="Skip XLSX human report generation.")
    parser.add_argument("--external-panel-name", default=None, help="Name recorded in run_summary.json for prospective external panels.")
    parser.add_argument("--metadata-xlsx", default=None, help="Optional metadata workbook passed to the raw V4 feature builder.")
    parser.add_argument("--version", action="version", version=f"probioscore {__version__}")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    results = score_input(
        out_dir=Path(args.out_dir),
        mode=args.mode,
        v4_tsv_dir=args.v4_tsv_dir,
        input_tsv=args.input_tsv,
        config_path=args.config,
        curated_tables_dir=args.curated_tables_dir,
        allow_missing_curated_tables=args.allow_missing_curated_tables,
        skip_visuals=args.skip_visuals,
        skip_xlsx=args.skip_xlsx,
        external_panel_name=args.external_panel_name,
        metadata_xlsx=args.metadata_xlsx,
    )
    print(f"ProbioScore completed: {len(results)} genomes scored.")
    print(f"Results: {Path(args.out_dir) / '2.probioscore_results' / 'probioscore_calculator_results.tsv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
