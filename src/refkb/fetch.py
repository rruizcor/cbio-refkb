"""Download selected cBioPortal tables from Hugging Face at a pinned revision.

Every fetch is pinned to a commit SHA of the HF dataset repo and records a
checksum for each file, so a build can be tied to an exact data snapshot
(needed for CAP/CLIA document control and reproducibility).
"""

from __future__ import annotations

import hashlib
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path

from huggingface_hub import HfApi, hf_hub_download

from . import HF_REPO, PHASE1_TABLES

METADATA_FILES = ["manifest.json", "schema.sql", "publication.json"]


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def table_files(all_files: list[str], table: str) -> list[str]:
    """Parquet files belonging to one table.

    The repo exposes each table as its own dataset configuration; we match any
    parquet file with a path segment equal to the table name, so the code does
    not depend on the exact folder layout.
    """
    out = []
    for f in all_files:
        if not f.endswith(".parquet"):
            continue
        parts = f.split("/")
        if table in parts[:-1] or parts[-1].startswith(f"{table}-") or parts[-1] == f"{table}.parquet":
            out.append(f)
    return sorted(out)


def fetch(
    out_root: Path,
    tables: list[str] | None = None,
    revision: str | None = None,
    token: str | None = None,
) -> Path:
    """Fetch tables into ``out_root/<sha>/<table>/`` and write ``fetch_manifest.json``.

    If ``revision`` is None the current ``main`` is resolved to its commit SHA
    first, so the download itself is always pinned.
    """
    tables = tables or PHASE1_TABLES
    api = HfApi(token=token)
    info = api.dataset_info(HF_REPO, revision=revision or "main")
    sha = info.sha
    all_files = api.list_repo_files(HF_REPO, repo_type="dataset", revision=sha)

    snap_dir = Path(out_root) / sha
    snap_dir.mkdir(parents=True, exist_ok=True)
    record: dict = {
        "hf_repo": HF_REPO,
        "revision_requested": revision or "main",
        "revision_sha": sha,
        "fetched_at_utc": datetime.now(timezone.utc).isoformat(),
        "tables": {},
        "metadata_files": {},
    }

    for meta in METADATA_FILES:
        if meta in all_files:
            p = Path(hf_hub_download(HF_REPO, meta, repo_type="dataset", revision=sha, token=token))
            dest = snap_dir / meta
            shutil.copyfile(p, dest)
            record["metadata_files"][meta] = _sha256(dest)

    missing = []
    for table in tables:
        files = table_files(all_files, table)
        if not files:
            missing.append(table)
            continue
        tdir = snap_dir / table
        tdir.mkdir(exist_ok=True)
        record["tables"][table] = []
        for i, f in enumerate(files):
            p = Path(hf_hub_download(HF_REPO, f, repo_type="dataset", revision=sha, token=token))
            dest = tdir / f"part-{i:05d}.parquet"
            shutil.copyfile(p, dest)
            record["tables"][table].append({"source": f, "file": dest.name, "sha256": _sha256(dest)})

    if missing:
        raise RuntimeError(
            f"No parquet files found for tables {missing} at revision {sha}. "
            "Check the repo layout with HfApi().list_repo_files()."
        )

    (snap_dir / "fetch_manifest.json").write_text(json.dumps(record, indent=2))
    return snap_dir
