# Technical installation and user guide for ProbioScore I45 portable

## 1. Purpose

This guide describes the installation, execution, validation and interpretation of the portable ProbioScore I45 calculator. The calculator is a scoring and governance layer for outputs produced by the upstream V4 probiogenomic pipeline. It does not assemble genomes, annotate genes or rerun antiSMASH, BAGEL5, AMR screening, virulence screening, toxin detection or ProbML.

The scientific objective is conservative triage. The tool identifies genomes that are prioritized, rejected, sent to safety review or marked as insufficiently supported by probiogenomic evidence.

## 2. Scientific status of this release

The release corresponds to I45. The final decision is `benchmark_curated_locked` for the embedded benchmarks: the curated internal panel and the Ounissi 48 pilot are solved under the expert curation rules. The prospective frozen mode has since been applied to an independent external panel of 634 genomes, reported in Gomri et al. (2026); those external results are presented in the manuscript and its supplementary tables rather than embedded as a benchmark in this package. Score any new external panel with `--mode prospective_frozen`.

The exact sentence that must accompany benchmark results is:

This benchmark performance was achieved after expert curation of reference labels and is not a measure of prospective generalization.

## 3. Package contents

The package contains:

- `probioscore/engine`: frozen scoring and I45 governance code.
- `probioscore/config`: frozen configs, AHP weights, schema and threshold lockfile.
- `probioscore/curated_tables`: five expert curated TSV files and a SHA256 manifest.
- `benchmarks`: I45 internal 1103 and Ounissi 48 reference outputs.
- `scripts`: manifest, freeze invariant, release ZIP and scoring helpers.
- `probioscore/tests`: pytest tests, mini fixture and non regression scripts.
- `examples`: ready commands for the future 600 genome panel and for single input scoring.

## 4. Installation with pip

From the portable calculator directory:

```bash
python -m pip install .
probioscore --help
```

For development or audit:

```bash
python -m pip install -e .
python scripts/run_acceptance_tests_E1_E15.py
```

## 5. Installation with conda

```bash
conda env create -f environment.yml
conda activate probioscore-i45
python -m pip install .
probioscore --help
```

## 6. Installation with Docker

```bash
docker build -t probioscore-i45 .
docker run --rm probioscore-i45 probioscore --help
```

For a real panel, mount the V4 output directory and an output directory:

```bash
docker run --rm \
  -v /path/to/v4_outputs:/data/v4_outputs:ro \
  -v /path/to/results:/data/results \
  probioscore-i45 \
  probioscore --v4-tsv-dir /data/v4_outputs --out-dir /data/results --mode prospective_frozen --skip-visuals
```

## 7. Windows portable launcher

```bat
run_probioscore.bat --help
```

Example:

```bat
run_probioscore.bat --input-tsv benchmarks\ounissi_48_reference\ounissi48_true_fce_results_I45.tsv --out-dir out\ounissi_check --mode benchmark_curated --skip-visuals
```

## 8. POSIX portable launcher

```bash
./run_probioscore.sh --help
```

Example:

```bash
./run_probioscore.sh --input-tsv benchmarks/ounissi_48_reference/ounissi48_true_fce_results_I45.tsv --out-dir out/ounissi_check --mode benchmark_curated --skip-visuals
```

## 9. Supported inputs

Two input modes are supported.

### 9.1 Raw V4 directory

Use `--v4-tsv-dir` when the directory contains the expected V4 subfolders:

- `final/final_feature_matrix.tsv`
- `final/final_module_summary.tsv`
- `probiosml/by_category/*.tsv`
- `safety/*.tsv`
- `toxins/*.tsv`
- `mobile_safety/*.tsv`
- `bgc/*.tsv`
- `crispr/*.tsv`
- `ta_systems/*.tsv`

The raw V4 path runs the embedded AHP plus FCE base engine first, then applies I45 governance.

### 9.2 Prebuilt TSV

Use `--input-tsv` when the table already contains at least:

- `genome_id`
- `ProbioScore_Status`

Additional fields such as risk scores, ProbML columns, taxonomy and memberships are preserved when present. This mode is useful for non regression, benchmark replay and manually inspected V4 derived tables.

## 10. Modes

### 10.1 prospective_frozen

This is the default mode and the intended mode for future external panels. It uses only frozen curated tables and fixed thresholds. It does not read external labels. It writes `Ounissi_used_for_calibration=False` in the output.

### 10.2 benchmark_curated

This mode is for reproducing curated benchmark metrics. It adds curated benchmark target fields and should not be interpreted as prospective generalization.

## 11. Output structure

The output directory contains:

```text
1.raw_v4_direct_feature_build/
  raw_v4_input_audit.tsv
2.probioscore_results/
  probioscore_calculator_results.tsv
  status_distribution.tsv
  probml_comparison.tsv
  run_summary.json
  RUN_MANIFEST.tsv
  probioscore_human_report.xlsx
  status_distribution.png
```

If `--skip-xlsx` is used, the workbook is not created. If `--skip-visuals` is used, the PNG chart is not created.

## 12. Interpretation of statuses

`Genomically-prioritized-probiotic-candidate` means that the genome passes the current in silico triage logic and has enough positive evidence after safety governance.

`Safety-review` means that the genome is not rejected but carries enough risk signal to require expert safety review.

`Risk-profile-rejected` means that the risk profile is incompatible with prioritization under the current rules.

`Insufficient-probiogenomic-evidence` means that the genome does not have enough host specific or probiotic evidence, or that an independent neutral environmental prior blocks priority.

## 13. ProbML governance

ProbML is never a decision source in this release. If ProbML columns are present, the calculator writes comparator summaries. ProbML cannot promote, rescue, reject or reclassify a genome.

## 14. Curated tables

The five curated tables are:

- `gold_standard_curated_labels.tsv`
- `neutral_or_environmental_taxon_prior.tsv`
- `curated_strain_level_probiotic_evidence.tsv`
- `curated_release_exception_v2.tsv`
- `risk_reclassification_curated.tsv`

The tables are frozen for this release and tracked by SHA256 in `FINAL_RELEASE_MANIFEST.tsv` and `probioscore/curated_tables/CURATED_TABLES_MANIFEST.tsv`.

## 15. Freeze invariant checks

Run:

```bash
python scripts/verify_freeze_invariants.py
```

Expected result: all listed lines must be PASS.

## 16. Acceptance suite

Run:

```bash
python scripts/run_acceptance_tests_E1_E15.py
```

Expected result: E1 to E15 must be PASS.

The suite checks engine immutability, AHP weights, lockfile, curated tables, internal 1103 benchmark invariants, Ounissi 48 benchmark invariants, mini fixture behavior, launcher banners, CI matrix declaration, freeze invariants and package completeness.

## 17. Reproducibility record

Every run writes a `RUN_MANIFEST.tsv` with SHA256 hashes for outputs. The release itself contains `FINAL_RELEASE_MANIFEST.tsv`.

## 18. Known limitation

This calculator is not a regulatory approval tool and not a clinical safety certificate. It is an in silico triage model. Any candidate selected by this tool still requires expert biological review, strain authentication, assembly quality review, mobile AMR interpretation and laboratory validation.
