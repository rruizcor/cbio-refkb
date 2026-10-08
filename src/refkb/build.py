"""Build Phase 1 reference tables from a fetched snapshot."""

from __future__ import annotations

import json
from datetime import datetime, timezone
from importlib import resources
from pathlib import Path

import duckdb

from . import PHASE1_TABLES, __version__


def _sql(name: str, **params) -> str:
    text = resources.files("refkb").joinpath(f"sql/{name}").read_text()
    return text.format(**params) if params else text


def connect(snapshot_dir: Path, db_path: str = ":memory:", tables: list[str] | None = None):
    """Open DuckDB with one view per fetched table."""
    snapshot_dir = Path(snapshot_dir)
    con = duckdb.connect(db_path)
    for table in tables or PHASE1_TABLES:
        tdir = snapshot_dir / table
        if not any(tdir.glob("*.parquet")):
            raise FileNotFoundError(f"No parquet files for table '{table}' in {tdir}")
        glob = str(tdir / "*.parquet").replace("'", "''")
        con.execute(f"CREATE OR REPLACE VIEW {table} AS SELECT * FROM read_parquet('{glob}', union_by_name=true)")
    con.execute(_sql("samples.sql"))
    con.execute(_sql("mutations.sql"))
    return con


def build(
    snapshot_dir: Path,
    out_root: Path,
    min_normal_depth: int = 10,
    threads: int | None = None,
    memory_limit: str | None = None,
) -> Path:
    snapshot_dir = Path(snapshot_dir)
    min_normal_depth = int(min_normal_depth)
    fetch_manifest_path = snapshot_dir / "fetch_manifest.json"
    fetch_manifest = json.loads(fetch_manifest_path.read_text()) if fetch_manifest_path.exists() else {}
    snapshot_id = fetch_manifest.get("revision_sha", snapshot_dir.name)

    out_dir = Path(out_root) / snapshot_id
    out_dir.mkdir(parents=True, exist_ok=True)

    con = connect(snapshot_dir)
    if threads:
        con.execute(f"SET threads = {int(threads)}")
    if memory_limit:
        # Cap RAM so large joins spill to disk instead of being killed.
        con.execute("SET memory_limit = ?", [memory_limit])
        spill = (out_dir / ".duckdb_tmp").as_posix()
        con.execute("SET temp_directory = ?", [spill])

    outputs = {}

    def export(name: str, query_or_table: str):
        dest = out_dir / f"{name}.parquet"
        con.execute(f"COPY ({query_or_table}) TO '{dest.as_posix()}' (FORMAT parquet, COMPRESSION zstd)")
        n = con.execute(f"SELECT count(*) FROM read_parquet('{dest.as_posix()}')").fetchone()[0]
        outputs[name] = {"file": dest.name, "rows": n}

    export("germline_profile_by_study", _sql("profile_germline.sql", min_normal_depth=min_normal_depth).rstrip().rstrip(";"))
    export("panel_profile", _sql("panel_profile.sql").rstrip().rstrip(";"))

    con.execute(_sql("recurrence.sql", min_normal_depth=min_normal_depth))
    export("variant_recurrence", "SELECT * FROM variant_recurrence ORDER BY tumor_type, gene, start_pos")

    manifest = {
        "refkb_version": __version__,
        "built_at_utc": datetime.now(timezone.utc).isoformat(),
        "snapshot_id": snapshot_id,
        "source": fetch_manifest.get("hf_repo"),
        "params": {"min_normal_depth": min_normal_depth},
        "outputs": outputs,
        "notes": [
            "Coordinates are as deposited (mixed GRCh37/GRCh38, see 'build' column); liftover is not yet applied.",
            "Annotation resource only: not a filter, not a validation truth set, not an actionability knowledgebase.",
        ],
    }
    (out_dir / "build_manifest.json").write_text(json.dumps(manifest, indent=2))
    return out_dir
