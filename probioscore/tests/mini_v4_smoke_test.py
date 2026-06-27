from __future__ import annotations

from pathlib import Path
import sys

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from probioscore.engine.consolidated_engine import score_dataframe


def main() -> int:
    fixture = Path(__file__).resolve().parent / "fixtures_v4_mini" / "mini_v4_input.tsv"
    expected_path = Path(__file__).resolve().parent / "fixtures_v4_mini" / "mini_v4_expected.tsv"
    observed = score_dataframe(pd.read_csv(fixture, sep="\t"), mode="prospective_frozen")
    expected = pd.read_csv(expected_path, sep="\t")
    merged = observed[["genome_id", "ProbioScore_Status"]].merge(expected, on="genome_id", suffixes=("_observed", "_expected"))
    ok = merged["ProbioScore_Status_observed"].eq(merged["ProbioScore_Status_expected"]).all()
    print(f"Mini fixture rows: {len(merged)}. PASS={ok}.")
    if not ok:
        print(merged.to_string(index=False))
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
