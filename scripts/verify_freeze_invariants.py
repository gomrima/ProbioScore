from __future__ import annotations

import hashlib
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]

CHECKS = [
    ("E1_engine_true_fce", "probioscore/engine/true_fce_engine.py"),
    ("E1_engine_consolidated", "probioscore/engine/consolidated_engine.py"),
    ("E2_ahp_weights", "probioscore/config/ahp_weights.json"),
    ("E3_lockfile", "probioscore/config/selected_thresholds_lockfile_v1.json"),
    ("E4_gold_standard", "probioscore/curated_tables/gold_standard_curated_labels.tsv"),
    ("E5_neutral_prior", "probioscore/curated_tables/neutral_or_environmental_taxon_prior.tsv"),
    ("E6_strain_evidence", "probioscore/curated_tables/curated_strain_level_probiotic_evidence.tsv"),
    ("E7_release_exception", "probioscore/curated_tables/curated_release_exception_v2.tsv"),
    ("E8_risk_reclassification", "probioscore/curated_tables/risk_reclassification_curated.tsv"),
]


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def run_checks() -> pd.DataFrame:
    manifest = pd.read_csv(ROOT / "FINAL_RELEASE_MANIFEST.tsv", sep="\t")
    expected = dict(zip(manifest["relative_path"], manifest["sha256"]))
    rows = []
    for check_id, rel in CHECKS:
        path = ROOT / rel
        observed = sha256_file(path) if path.is_file() else ""
        target = expected.get(rel, "")
        rows.append(
            {
                "check_id": check_id,
                "relative_path": rel,
                "observed_sha256": observed,
                "manifest_sha256": target,
                "status": "PASS" if observed and observed == target else "FAIL",
            }
        )
    return pd.DataFrame(rows)


def main() -> int:
    df = run_checks()
    out = ROOT / "FREEZE_INVARIANTS_E1_E8.tsv"
    df.to_csv(out, sep="\t", index=False)
    print(df.to_string(index=False))
    print(f"Wrote {out}")
    return 0 if df["status"].eq("PASS").all() else 1


if __name__ == "__main__":
    raise SystemExit(main())
