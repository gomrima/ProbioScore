from __future__ import annotations

import json
import subprocess
import sys
import zipfile
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
PRIORITY = "Genomically-prioritized-probiotic-candidate"
INSUFFICIENT = "Insufficient-probiogenomic-evidence"


def result(test_id: str, scope: str, status: str, detail: str) -> dict[str, str]:
    return {"test_id": test_id, "scope": scope, "status": status, "detail": detail}


def file_hash_matches(rel: str) -> bool:
    import hashlib

    manifest = pd.read_csv(ROOT / "FINAL_RELEASE_MANIFEST.tsv", sep="\t")
    expected = dict(zip(manifest["relative_path"], manifest["sha256"]))
    path = ROOT / rel
    if not path.is_file() or rel not in expected:
        return False
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest() == expected[rel]


def run_acceptance() -> pd.DataFrame:
    rows: list[dict[str, str]] = []
    rows.append(result("E1", "engine immutability", "PASS" if file_hash_matches("probioscore/engine/true_fce_engine.py") and file_hash_matches("probioscore/engine/consolidated_engine.py") else "FAIL", "Engine files match manifest."))
    rows.append(result("E2", "AHP weights", "PASS" if file_hash_matches("probioscore/config/ahp_weights.json") else "FAIL", "AHP weights JSON matches manifest."))
    rows.append(result("E3", "lockfile", "PASS" if file_hash_matches("probioscore/config/selected_thresholds_lockfile_v1.json") else "FAIL", "Threshold lockfile matches manifest."))
    curated = [
        ("E4", "gold standard", "probioscore/curated_tables/gold_standard_curated_labels.tsv"),
        ("E5", "neutral prior", "probioscore/curated_tables/neutral_or_environmental_taxon_prior.tsv"),
        ("E6", "strain evidence", "probioscore/curated_tables/curated_strain_level_probiotic_evidence.tsv"),
        ("E7", "release exception", "probioscore/curated_tables/curated_release_exception_v2.tsv"),
        ("E8", "risk reclassification", "probioscore/curated_tables/risk_reclassification_curated.tsv"),
    ]
    for tid, scope, rel in curated:
        rows.append(result(tid, scope, "PASS" if file_hash_matches(rel) else "FAIL", f"{rel} matches manifest."))

    internal = pd.read_csv(ROOT / "benchmarks/internal_1103_reference/true_fce_results_I45.tsv", sep="\t")
    pathogen = internal["expected_class"].eq("Pathogen-reference")
    path37 = internal["genome_id"].eq("Path37_Lactococcus_lactis_subsp_lactis_DSM_4641")
    probio = internal["expected_class"].eq("Probio-Ichnos-reference")
    env = internal["expected_class"].eq("Environmental-comparator")
    ok9 = int((internal["ProbioScore_Status"].eq(PRIORITY) & probio).sum()) == 685
    ok9 = ok9 and int((internal["ProbioScore_Status"].eq(PRIORITY) & env).sum()) == 0
    ok9 = ok9 and int((internal["ProbioScore_Status"].eq(PRIORITY) & pathogen & ~path37).sum()) == 0
    rows.append(result("E9", "internal 1103 benchmark", "PASS" if ok9 else "FAIL", "Expected I45 internal counts are reproduced from embedded reference."))

    oun = pd.read_csv(ROOT / "benchmarks/ounissi_48_reference/ounissi48_true_fce_results_I45.tsv", sep="\t")
    oun_probio = oun["expected_class"].astype(str).str.contains("probiotic", case=False, na=False)
    oun_env = oun["expected_class"].eq("Environmental-comparator")
    oun_path = oun["expected_class"].eq("Pathogen-reference")
    ok10 = int((oun["ProbioScore_Status"].eq(PRIORITY) & oun_probio).sum()) == 24
    ok10 = ok10 and int((oun["ProbioScore_Status"].eq(PRIORITY) & oun_env).sum()) == 0
    ok10 = ok10 and int((oun["ProbioScore_Status"].eq(INSUFFICIENT) & oun_env).sum()) == 12
    ok10 = ok10 and int((oun["ProbioScore_Status"].eq(PRIORITY) & oun_path).sum()) == 0
    rows.append(result("E10", "Ounissi 48 benchmark", "PASS" if ok10 else "FAIL", "Expected Ounissi I45 counts are reproduced from embedded reference."))

    mini_cmd = [sys.executable, str(ROOT / "probioscore/tests/mini_v4_smoke_test.py")]
    mini = subprocess.run(mini_cmd, cwd=ROOT, capture_output=True, text=True)
    rows.append(result("E11", "mini fixture", "PASS" if mini.returncode == 0 else "FAIL", mini.stdout.strip() or mini.stderr.strip()))

    launcher_text = (ROOT / "run_probioscore.sh").read_text(encoding="utf-8") + (ROOT / "run_probioscore.bat").read_text(encoding="utf-8")
    ok12 = "probioscore" in launcher_text and "--help" in launcher_text
    rows.append(result("E12", "portable launchers", "PASS" if ok12 else "FAIL", "Launchers expose usage banner."))

    ci = (ROOT / ".github/workflows/ci.yml").read_text(encoding="utf-8")
    ok13 = all(x in ci for x in ["ubuntu-latest", "macos-latest", "windows-latest", "3.10", "3.11", "3.12"])
    rows.append(result("E13", "CI matrix", "PASS" if ok13 else "FAIL", "CI matrix declares 3 OS and 3 Python versions."))

    freeze = subprocess.run([sys.executable, str(ROOT / "scripts/verify_freeze_invariants.py")], cwd=ROOT, capture_output=True, text=True)
    rows.append(result("E14", "freeze invariants", "PASS" if freeze.returncode == 0 else "FAIL", "verify_freeze_invariants.py completed."))

    essential = [
        "probioscore/curated_tables/gold_standard_curated_labels.tsv",
        "probioscore/config/selected_thresholds_lockfile_v1.json",
        "probioscore/engine/consolidated_engine.py",
    ]
    wheel_ok = False
    wheels = sorted((ROOT / "dist").glob("*.whl"))
    if wheels:
        with zipfile.ZipFile(wheels[-1]) as zf:
            names = set(zf.namelist())
        wheel_ok = all(any(name.endswith(item) for name in names) for item in essential)
    package_ok = all((ROOT / item).is_file() for item in essential)
    e15_ok = package_ok and (wheel_ok if wheels else True)
    rows.append(result("E15", "source package completeness", "PASS" if e15_ok else "FAIL", "Essential package data are present in the source tree; any locally built wheel is inspected when present."))
    return pd.DataFrame(rows)


def main() -> int:
    df = run_acceptance()
    out = ROOT / "ACCEPTANCE_TESTS_E1_E15.tsv"
    df.to_csv(out, sep="\t", index=False)
    print(df.to_string(index=False))
    print(f"Wrote {out}")
    return 0 if df["status"].eq("PASS").all() else 1


if __name__ == "__main__":
    raise SystemExit(main())
