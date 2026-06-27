from __future__ import annotations

import hashlib
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def iter_files() -> list[Path]:
    skip_dirs = {"__pycache__", ".pytest_cache", ".mypy_cache", ".ruff_cache", "build", "dist", "htmlcov", "out", "probioscore_i45.egg-info", ".venv", "venv", ".git", ".idea", ".vscode"}
    skip_files = {"PROBIOSCORE_I45_PORTABLE_ZIP_SHA256.tsv"}
    files = []
    for path in ROOT.rglob("*"):
        if not path.is_file():
            continue
        rel = path.relative_to(ROOT)
        if set(rel.parts) & skip_dirs:
            continue
        if path.name in skip_files or path.name.startswith("~$") or path.suffix.lower() in {".pyc", ".pyo", ".tmp"}:
            continue
        files.append(path)
    return sorted(files)


def main() -> int:
    zip_path = ROOT.parent / f"{ROOT.name}.zip"
    if zip_path.exists():
        zip_path.unlink()
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for path in iter_files():
            zf.write(path, (Path(ROOT.name) / path.relative_to(ROOT)).as_posix())
    sidecar = ROOT / "PROBIOSCORE_I45_PORTABLE_ZIP_SHA256.tsv"
    sidecar.write_text("zip_file\tbytes\tsha256\n" + f"{zip_path.name}\t{zip_path.stat().st_size}\t{sha256_file(zip_path)}\n", encoding="utf-8")
    print(f"Wrote {zip_path}")
    print(f"SHA256 {sha256_file(zip_path)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
