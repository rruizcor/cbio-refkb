"""Command line: refkb fetch | build | query."""

from __future__ import annotations

import argparse
import os
from pathlib import Path


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(prog="refkb", description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)

    f = sub.add_parser("fetch", help="Download Phase 1 tables from Hugging Face at a pinned revision")
    f.add_argument("--out", default="data/raw", type=Path)
    f.add_argument("--revision", default=None, help="HF commit SHA to pin (default: resolve current main)")
    f.add_argument("--tables", nargs="*", default=None)

    b = sub.add_parser("build", help="Build reference tables from a fetched snapshot")
    b.add_argument("snapshot_dir", type=Path)
    b.add_argument("--out", default="builds", type=Path)
    b.add_argument("--min-normal-depth", default=10, type=int)
    b.add_argument("--threads", default=None, type=int)
    b.add_argument("--memory-limit", default=None, help="DuckDB RAM cap, e.g. 12GB (spills to disk beyond it)")

    q = sub.add_parser("query", help="Run ad hoc SQL against a fetched snapshot (views per table)")
    q.add_argument("snapshot_dir", type=Path)
    q.add_argument("sql")

    args = ap.parse_args(argv)

    if args.cmd == "fetch":
        from .fetch import fetch

        snap = fetch(args.out, tables=args.tables, revision=args.revision, token=os.environ.get("HF_TOKEN"))
        print(f"Snapshot written to {snap}")
    elif args.cmd == "build":
        from .build import build

        out = build(args.snapshot_dir, args.out, min_normal_depth=args.min_normal_depth,
                    threads=args.threads, memory_limit=args.memory_limit)
        print(f"Build written to {out}")
    elif args.cmd == "query":
        from .build import connect

        print(connect(args.snapshot_dir).sql(args.sql))


if __name__ == "__main__":
    main()
