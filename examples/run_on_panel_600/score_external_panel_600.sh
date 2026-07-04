#!/usr/bin/env bash
set -euo pipefail
probioscore --v4-tsv-dir "${1:?V4 TSV directory required}" --out-dir "${2:-panel_600_results}" --mode prospective_frozen --external-panel-name panel_600_v4 --skip-visuals
