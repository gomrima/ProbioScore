from __future__ import annotations

from pathlib import Path

import pandas as pd

from probioscore.engine.consolidated_engine import PRIORITY, INSUFFICIENT, score_dataframe


def test_mini_fixture_reproduces_expected_statuses() -> None:
    root = Path(__file__).resolve().parent / "fixtures_v4_mini"
    observed = score_dataframe(pd.read_csv(root / "mini_v4_input.tsv", sep="\t"), mode="prospective_frozen")
    expected = pd.read_csv(root / "mini_v4_expected.tsv", sep="\t")
    merged = observed[["genome_id", "ProbioScore_Status"]].merge(expected, on="genome_id", suffixes=("_observed", "_expected"))
    assert merged["ProbioScore_Status_observed"].eq(merged["ProbioScore_Status_expected"]).all()


def test_membership_partition_is_exact_after_governance() -> None:
    root = Path(__file__).resolve().parent / "fixtures_v4_mini"
    observed = score_dataframe(pd.read_csv(root / "mini_v4_input.tsv", sep="\t"), mode="prospective_frozen")
    assert float(observed["max_partition_drift"].max()) == 0.0
    assert PRIORITY in set(observed["ProbioScore_Status"])
    assert INSUFFICIENT in set(observed["ProbioScore_Status"])
