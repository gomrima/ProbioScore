# ProbioScore portable calculator

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python](https://img.shields.io/badge/python-3.10%20%7C%203.11%20%7C%203.12-blue.svg)](https://www.python.org/downloads/)
[![Release](https://img.shields.io/badge/release-v1.0%20%28CR0.10%29-brightgreen.svg)](CHANGELOG.md)
[![Status](https://img.shields.io/badge/status-prospective--frozen-orange.svg)](CHANGELOG.md)

ProbioScore is a deterministic probiogenomic triage calculator. It
applies a decoupled AHP and FCE decision layer to V4 pipeline outputs and
then applies a frozen I45 governance layer built from five expert curated
tables.

> ProbioScore is an in silico research triage tool. It is not a
> clinical safety certificate, not a regulatory approval, and not a
> substitute for expert biological review.


## Why ProbioScore portable

| feature | ProbioScore portable |
| ------- | ----------------------------- |
| portable repository layout | yes |
| Linux, macOS, Windows CI matrix included | yes |
| pip and conda installers | yes |
| Docker image | yes |
| GitHub Actions CI matrix | yes |
| five-genome embedded smoke test | yes |
| embedded reference outputs (1103 internal + 48 Ounissi) | yes |
| five frozen curated I45 tables | yes |
| acceptance suite E1 to E15 | yes |
| open licence | MIT |
| ProbML strictly comparator only | yes |
| benchmark curated and prospective frozen modes separated | yes |

## What you get per genome

- `ProbioScore_Status`: one of four operational statuses.
- `ProbioScore_Utility`: AHP weighted positive evidence aggregate.
- `FCE_Membership_*` columns for the four statuses and `max_partition_drift`.
- `mu_Risk`, `mu_Pathogen_Total`, `mu_Benefit_Final` continuous memberships.
- `external_neutral_taxon_prior`, `strain_level_override_applied`,
  `release_exception_v2_applied`, `risk_reclassification_*` governance
  fields.
- `decision_rule_trace` describing the active rule.
- `probml_priority` reported as comparator only when present.
- `RUN_MANIFEST.tsv` listing every output with SHA256.

The four operational statuses are:

| status | meaning |
| ------ | ------- |
| `Genomically-prioritized-probiotic-candidate` | enough positive evidence and acceptable safety burden after governance |
| `Safety-review` | safety burden that requires expert review |
| `Risk-profile-rejected` | safety burden incompatible with prioritisation |
| `Insufficient-probiogenomic-evidence` | positive evidence too weak to prioritise or blocked by independent neutral environmental prior |

## Quick start

```bash
# From GitHub
git clone https://github.com/gomrima/ProbioScore.git
cd ProbioScore
# or, from the portable zip:
# unzip ProbioScore_I45_Final_Executable_Calculator_V1_portable.zip
# cd ProbioScore_I45_Final_Executable_Calculator_V1_portable
python3 -m venv .venv
source .venv/bin/activate          # Linux and macOS
# .venv\Scripts\activate           # Windows powershell
python -m pip install --upgrade pip
python -m pip install .
probioscore --help
```

Run the embedded mini fixture smoke test:

```bash
python3 probioscore/tests/mini_v4_smoke_test.py
```

Expected outcome: `Mini fixture rows: 5. PASS=True.` and exit code 0.

Run the full acceptance suite (E1 to E15):

```bash
python3 scripts/run_acceptance_tests_E1_E15.py
```

Expected outcome: 15 PASS lines and exit code 0.

Run on your own V4 pipeline outputs:

```bash
./run_probioscore.sh --v4-tsv-dir /path/to/v4_tsv_outputs --out-dir ./my_results --mode prospective_frozen --skip-visuals
```

On Windows:

```bat
run_probioscore.bat --v4-tsv-dir C:\path\to\v4_tsv_outputs --out-dir .\my_results --mode prospective_frozen --skip-visuals
```

The full command-line reference is in
[TECHNICAL_INSTALLATION_AND_USER_GUIDE.md](TECHNICAL_INSTALLATION_AND_USER_GUIDE.md).

## Installation alternatives

- Conda: `conda env create -f environment.yml && conda activate probioscore-i45 && pip install .`
- Docker: `docker build -t probioscore-i45:1.0 . && docker run --rm probioscore-i45:1.0 probioscore --help`

## Reproducibility

- Deterministic random seed 42.
- Maximum partition drift observed: 0.0 (machine precision exact).
- Five curated TSV tables and one threshold lockfile tracked by SHA256.
- Engine modules tracked by SHA256.
- Embedded benchmark reference outputs for 1103 internal and 48 Ounissi.

Run the freeze invariants check locally:

```bash
python3 scripts/verify_freeze_invariants.py
```

Expected outcome: 9 PASS lines.

## Frozen artefact hashes

- Engine `true_fce_engine.py`: `d2cb2d9daa98021d2242c31f4fb19d510acd01796f42fc894df835b8712a8ed7`
- Engine `consolidated_engine.py`: `c06a45cdfbcff7481941d57757d2b49f72270c5a26496dec3620fcd6d8d2ab40`
- `ahp_weights.json`: `5b122974a3fecfcfe5cd866af876568a42a33b3821a5edc9fc4afd473492e24f`
- `selected_thresholds_lockfile_v1.json`: `404e39210cb7218456f22d3917cfc288d1cc8fd0d48156d58876a4fb28b43c85`
- `gold_standard_curated_labels.tsv`: `a034f7d0549208c3aa4ec52d1c6d24fa8ee70e9df8e77522e8b1d31a7f38896e`
- `neutral_or_environmental_taxon_prior.tsv`: `73749df3e8475962438d957b41d5dc760048e848ec85e29eded14d35d4826ed4`
- `curated_strain_level_probiotic_evidence.tsv`: `0cfc02afb2b33e7d4cbf9e584c75243f4dbd268c2dc43b5e86dcfe32b2646161`
- `curated_release_exception_v2.tsv`: `2c6dddb9ee11b6953516a2be393a6c0a5eb9c0ebbed39c433e7392e5358eb358`
- `risk_reclassification_curated.tsv`: `09375fefe105f9aad7a098f3824668f37e5a28820793cb77312eb368570caa75`

## Methodology in one paragraph

ProbioScore reuses the consolidated decoupled AHP plus FCE
scoring engine. The positive utility axis aggregates ten
ProbioSML sub-pillars (B1 to B5 survival, C1 to C5 benefit) with AHP
weights from the CURATED_AGG_WEIGHTS_CR0.10 expert panel (51 experts,
section consistency ratio at or below 0.10; see
`probioscore/docs/WEIGHTS_PROVENANCE.md`). The safety axis evaluates four
risk dimensions (A1 antimicrobial resistance, A2 virulence, A3 toxins,
A4 mobile safety) as a non-compensatory veto rather than an AHP weighted sum. A true fuzzy
comprehensive evaluation produces continuous memberships
(`FCE_Membership_*`) whose partition exactly sums to 1.0. The I45
governance layer then applies, in order: an independent neutral
environmental taxon prior (12 genera), a strain-level probiotic evidence
override (currently *Escherichia coli* Nissle 1917 by exact accession),
a release exception v2 for 13 EFSA QPS Bacillus and Shouchella strains,
and a curated risk reclassification on the 91 historic Probio-Ichnos
non-priority cases. ProbML is recorded only as a comparator. Detailed
methodology is in `probioscore/docs/`.

## Repository layout

```text
ProbioScore_I45_Final_Executable_Calculator_V1_portable/
|-- README.md
|-- TECHNICAL_INSTALLATION_AND_USER_GUIDE.md
|-- CHANGELOG.md
|-- CONTRIBUTING.md
|-- CITATION.cff
|-- LICENSE
|-- pyproject.toml
|-- environment.yml
|-- Dockerfile
|-- MANIFEST.in
|-- run_probioscore.sh
|-- run_probioscore.bat
|-- probioscore/
|   |-- cli/
|   |-- config/
|   |-- curated_tables/
|   |-- docs/                (METHODOLOGY_*, SCHEMA_CONTRACT_V4, WEIGHTS_PROVENANCE)
|   |-- engine/
|   `-- tests/
|       `-- fixtures_v4_mini/
|-- benchmarks/
|   |-- internal_1103_reference/
|   `-- ounissi_48_reference/
|-- examples/
|   |-- run_on_panel_600/
|   `-- run_on_single_v4_output/
|-- scripts/
|   |-- build_release_zip.py
|   |-- regenerate_manifest.py
|   |-- run_acceptance_tests_E1_E15.py
|   |-- score_external_v4_outputs.py
|   `-- verify_freeze_invariants.py
|-- .github/workflows/ci.yml
|-- ACCEPTANCE_TESTS_E1_E15.tsv
|-- FREEZE_INVARIANTS_E1_E8.tsv
`-- FINAL_RELEASE_MANIFEST.tsv
```

## Citation

If you use ProbioScore in your work, please cite:

```
Gomri M.A., Ounissi N.E., El Hadef El Okki M. (2026).
ProbioScore: an Analytic Hierarchy Process and Fuzzy Comprehensive Evaluation Framework 
for Safety-Aware Probiotic Candidate Triage from Bacterial Genomes
```


## License
[MIT License](LICENSE). The licence text includes a scientific disclaimer
that ProbioScore produces operational categories and is not a
clinical safety certificate.

## Contributing

Contributions to non-engine code, documentation and tests are welcome.
The I45 engine and curated tables are frozen for this release. See
[CONTRIBUTING.md](CONTRIBUTING.md).

## Project authors

| role | name |
| ---- | ---- |
| Project lead | Mohamed Amine Gomri |
| Contributor | Nour ElHouda Ounissi |
| Contributor | Mohamed El Hadef El Okki |.
