# Changelog

All notable changes to ProbioScore portable are documented in
this file.

The project adheres to a freeze-driven release process. Once a version is
frozen, the engine SHA256, the AHP weights SHA256 and the five curated
table SHA256s are locked. Later releases re-use the frozen engine and
only add packaging, documentation, tests, fixtures or CI infrastructure.

---

## [1.0] - 2026-06-27

First public GitHub release (CR0.10 expert weighting). Supersedes the internal
1.0.0 first portable build (2026-05-27) recorded below; the I45 engine and the
five curated tables are unchanged and remain frozen by SHA256.

### Changed
- AHP weights recalibrated to the `CURATED_AGG_WEIGHTS_CR0.10` panel
  (AHP generator v8, 51 experts, only section matrices with consistency
  ratio at or below 0.10 retained, exact principal eigenvector
  aggregation). Survival sub-pillars (B1 to B5) and benefit sub-pillars
  (C1 to C5) are taken directly; the two positive pillars are renormalized
  to sum to 1 (Survival 0.580462, Benefits 0.419538) because safety is a
  non-compensatory veto outside the positive utility.
- Configuration now records the weight provenance and the full
  three-pillar panel weights including safety, for traceability.

### Validated
- The recalibration changes the continuous positive utility only
  marginally (mean absolute change about 0.003 to 0.004, maximum about
  0.022) and produces zero status reassignments on the internal 1,103 and
  external 634 panels. All four-class statuses, and the ProbioScore versus
  ProbML comparison, are unchanged. See `probioscore/docs/WEIGHTS_PROVENANCE.md`.

### Cleaned
- Pruned the package to runtime and reproducibility essentials. Removed
  the historical build report, release sidecar readme, zip checksum
  sidecar, the superseded freeze decision sidecar and the I31 to I45
  documentation tree.

## [1.0.0] - 2026-05-27

### Added
- First portable release of ProbioScore calculator.
- Public CLI entry point `probioscore`.
- pip installable Python package `probioscore-i45`.
- Conda environment file `environment.yml`.
- Docker image with non-root user, build-time smoke test and
  OpenContainers labels (`python:3.12-slim` base).
- GitHub Actions CI matrix (`ubuntu-latest`, `macos-latest`,
  `windows-latest`) on Python 3.10, 3.11, 3.12.
- POSIX launcher `run_probioscore.sh` and Windows launcher
  `run_probioscore.bat`.
- Five frozen curated I45 tables under
  `probioscore/curated_tables/`:
  - `gold_standard_curated_labels.tsv` (1151 entries, Path37 reclassified)
  - `neutral_or_environmental_taxon_prior.tsv` (12 genera, independent)
  - `curated_strain_level_probiotic_evidence.tsv` (Nissle 1917 by accession)
  - `curated_release_exception_v2.tsv` (13 EFSA QPS strains)
  - `risk_reclassification_curated.tsv` (91 Probio-Ichnos cases)
- Threshold lockfile `selected_thresholds_lockfile_v1.json`.
- AHP weights `ahp_weights.json` and base FCE config
  `ahp_fce_base_config.json`.
- Embedded reference outputs for 1103 internal and 48 Ounissi.
- Mini V4 fixture under `probioscore/tests/fixtures_v4_mini/`.
- Acceptance suite E1 to E15 in
  `scripts/run_acceptance_tests_E1_E15.py`.
- Freeze invariants verifier in `scripts/verify_freeze_invariants.py`.
- Manifest regenerator in `scripts/regenerate_manifest.py`.
- Release ZIP builder in `scripts/build_release_zip.py`.
- External V4 scoring helper in `scripts/score_external_v4_outputs.py`.
- pytest suite in `probioscore/tests/` (engine + curated tables).
- Documentation: `README.md`, `TECHNICAL_INSTALLATION_AND_USER_GUIDE.md`,
  `CONTRIBUTING.md`, `CITATION.cff`, `LICENSE`,
  `RELEASE_SIDECARS_README.md`,
  `PORTABLE_CALCULATOR_BUILD_AND_VALIDATION_REPORT.md`,
  `probioscore/docs/SCHEMA_CONTRACT_V4.md`,
  `probioscore/docs/METHODOLOGY_PROSPECTIVE_FROZEN.md`,
  `probioscore/docs/METHODOLOGY_BENCHMARK_CURATED.md`.

### Reproducibility
- Random seed 42.
- Maximum partition drift observed: 0.0 (machine precision exact).
- Status distribution invariants captured by the embedded reference TSV
  files in `benchmarks/`.

### Scientific status
- Decision: `benchmark_curated_locked`.
- `prospective_release_candidate`: False.
- `external validation completed on 634 genomes in the 1.0 release`.
- ProbML is reported only as a comparator and never as a decision input.

### Validation
- Acceptance suite E1 to E15: 15 PASS.
- pytest: 4 PASS.
- Freeze invariants E1 to E8: 9 file checks PASS.

---

## I31 to I45 scientific lineage

The portable release re-uses the I31 to I45 consolidated lineage built
on top of the I00 to I30 baseline:

- **I31**: gold standard curated labels (reclassification of Path37 only,
  with biological justification documented).
- **I32**: robust taxonomy parsing for any V4 input.
- **I33**: independent neutral environmental taxon prior (12 genera,
  curated from BacDive and primary literature, no use of Ounissi labels).
- **I34**: strain-level probiotic evidence (Nissle 1917 by exact
  accession, `non_overrideable_risks` enforced).
- **I35**: forensic recovery of safe Probio-Ichnos (13 EFSA QPS Bacillus
  and Shouchella cases).
- **I36**: intrinsic versus mobile risk reclassification (91 cases).
- **I37**: host-specific quality gate hardening.
- **I38**: benchmark curated mode (label-aware evaluation only).
- **I39**: prospective frozen mode (label-blind production mode).
- **I40**: external 600 V4 panel registered as `pending_v4_outputs`.
- **I41**: strict taxon-aware robustness (LOGO + LOFO + bootstrap +
  Monte Carlo AHP + ablations).
- **I42**: ProbML discordance governance (comparator only, 88.4 percent
  concordance, kappa 0.760, 0 discordant priority safety cases).
- **I43**: schema and anti-leak hardening (permutation, anonymisation,
  ProbML removal tests).
- **I44**: scientific documentation release candidate.
- **I45**: final decision `benchmark_curated_locked`.

---

## Versioning policy

- Major.Minor.Patch (semantic versioning).
- The portable release carries its own version (`1.0.0`) and is linked
  to the scientific lineage (`I45`).
- Any change to the engine, AHP weights, lockfile or curated tables
  triggers a new scientific release line.
- Packaging-only changes increment the patch number.

---

## Authors

| role | name |
| ---- | ---- |
| Project lead | Nour El Houda Ounissi |
| Supervisor | Mohamed Amine Gomri |
| Co-supervisor | Mohamed El Hadef El Okki |
