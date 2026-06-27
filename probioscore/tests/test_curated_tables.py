from __future__ import annotations

from importlib import resources

import pandas as pd


def table_path(name: str):
    return resources.files("probioscore").joinpath("curated_tables", name)


def test_curated_tables_exist_and_have_rows() -> None:
    expected = {
        "gold_standard_curated_labels.tsv": 1151,
        "neutral_or_environmental_taxon_prior.tsv": 12,
        "curated_strain_level_probiotic_evidence.tsv": 1,
        "curated_release_exception_v2.tsv": 13,
        "risk_reclassification_curated.tsv": 91,
    }
    for name, n in expected.items():
        df = pd.read_csv(table_path(name), sep="\t")
        assert len(df) == n


def test_probml_is_not_a_decision_source() -> None:
    lock = resources.files("probioscore").joinpath("config", "selected_thresholds_lockfile_v1.json").read_text()
    assert '"probml_used_for_decision": false' in lock
