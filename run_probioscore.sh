#!/usr/bin/env bash
set -euo pipefail
if [ "${1:-}" = "" ]; then
  probioscore --help
else
  probioscore "$@"
fi
