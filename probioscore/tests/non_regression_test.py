from __future__ import annotations

from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
PRIORITY = "Genomically-prioritized-probiotic-candidate"


def main() -> int:
    internal = pd.read_csv(ROOT / "benchmarks/internal_1103_reference/true_fce_results_I45.tsv", sep="\t")
    probio = internal["expected_class"].eq("Probio-Ichnos-reference")
    path = internal["expected_class"].eq("Pathogen-reference")
    path37 = internal["genome_id"].eq("Path37_Lactococcus_lactis_subsp_lactis_DSM_4641")
    env = internal["expected_class"].eq("Environmental-comparator")
    checks = {
        "internal_probio_priority_n": int((internal["ProbioScore_Status"].eq(PRIORITY) & probio).sum()) == 685,
        "internal_curated_pathogen_priority_n": int((internal["ProbioScore_Status"].eq(PRIORITY) & path & ~path37).sum()) == 0,
        "internal_environmental_priority_n": int((internal["ProbioScore_Status"].eq(PRIORITY) & env).sum()) == 0,
    }
    print(checks)
    return 0 if all(checks.values()) else 1


if __name__ == "__main__":
    raise SystemExit(main())
