from __future__ import annotations

import hashlib
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SKIP_DIRS = {"__pycache__", ".pytest_cache", ".mypy_cache", ".ruff_cache", "build", "dist", "htmlcov", "out", "probioscore_i45.egg-info", ".venv", "venv", ".git", ".idea", ".vscode"}
SKIP_FILES = {"FINAL_RELEASE_MANIFEST.tsv", "PROBIOSCORE_I45_PORTABLE_ZIP_SHA256.tsv"}


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def iter_files() -> list[Path]:
    files = []
    for path in ROOT.rglob("*"):
        if not path.is_file():
            continue
        rel_parts = set(path.relative_to(ROOT).parts)
        if rel_parts & SKIP_DIRS:
            continue
        if path.name in SKIP_FILES or path.name.startswith("~$"):
            continue
        if path.suffix.lower() in {".pyc", ".pyo", ".tmp"}:
            continue
        files.append(path)
    return sorted(files)


def main() -> int:
    rows = []
    for path in iter_files():
        rows.append({"relative_path": path.relative_to(ROOT).as_posix(), "bytes": path.stat().st_size, "sha256": sha256_file(path)})
    out = ROOT / "FINAL_RELEASE_MANIFEST.tsv"
    pd.DataFrame(rows).to_csv(out, sep="\t", index=False)
    print(f"Wrote {out}")
    print(f"Files indexed: {len(rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
