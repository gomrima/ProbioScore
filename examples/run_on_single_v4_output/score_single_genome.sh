#!/usr/bin/env bash
set -euo pipefail
probioscore --input-tsv "${1:?Input TSV required}" --out-dir "${2:-single_result}" --mode prospective_frozen --skip-visuals
