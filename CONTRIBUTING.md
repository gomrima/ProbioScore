# Contributing

This portable calculator is a frozen scientific release. Contributions should not change the I45 decision logic unless a new version is explicitly created.

Acceptable changes for this release line:

- Documentation corrections.
- Packaging fixes.
- CI fixes that do not change results.
- Additional tests that confirm current behavior.

Changes that require a new scientific version:

- Any threshold change.
- Any curated table change.
- Any change to status assignment.
- Any change that uses ProbML as a decision source.
- Any change to benchmark interpretation.

Before submitting a change, run:

```bash
python scripts/regenerate_manifest.py
python scripts/verify_freeze_invariants.py
python scripts/run_acceptance_tests_E1_E15.py
```
