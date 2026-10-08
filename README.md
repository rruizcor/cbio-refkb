# cbio-refkb

A versioned cancer-cohort reference layer built from the public cBioPortal database
([Hugging Face mirror](https://huggingface.co/datasets/cBioPortal/publicDatabase)), used to
annotate somatic calls from both the Illumina WES/WTS (Agilent capture) pipeline and the ONT
somatic pipeline.

It is an **annotation resource**: population priors and context. It is not a call filter,
not a validation truth set, and not an actionability knowledgebase.

## Phase 1 (this release)

| Output | What it is |
|---|---|
| `germline_profile_by_study.parquet` | Per study: germline/somatic/LOH-flagged calls, % of calls with matched-normal read counts, normal-confirmed somatic calls, genome builds. Run first on every snapshot. |
| `panel_profile.parquet` | Denominator sanity checks: panels in use, mutated samples with no panel row, how tumor type was assigned. |
| `variant_recurrence.parquet` | Patient-level, panel-aware SNV/indel recurrence by tumor type (+ `ALL`), with normal-confirmed, normal-alt-read and germline-flag counts. |
| `build_manifest.json` | Snapshot SHA, code version, parameters, row counts. |

Key definitions:
- **Tumor type**: sample-level `ONCOTREE_CODE` when available, else the study's `type_of_cancer_id` (upper-cased).
- **Denominator**: distinct patients whose sequencing covered the gene. `WES` samples cover all genes; panel samples only their panel genes; patients with both are counted once.
- **Normal-confirmed somatic**: not germline-flagged, `normalAltCount = 0`, `normalRefCount >= --min-normal-depth` (default 10).

## Usage

```bash
pip install -e ".[dev]"

# 1. Pull the Phase 1 tables at a pinned revision (resolves main -> commit SHA)
refkb fetch --out data/raw                    # or --revision <sha>

# 2. Build
refkb build data/raw/<sha> --out builds --min-normal-depth 10

# 3. Ad hoc questions against the raw snapshot
refkb query data/raw/<sha> "SELECT * FROM samples_base LIMIT 5"
```

Phase 1 pulls 7 of the 60 tables, a small fraction of the full 87 GB dump. Set `HF_TOKEN` if you hit rate limits.

## Data and licensing

- No data is committed to this repo (`data/` and `builds/` are git-ignored).
- The source database is distributed under **AGPL-3.0**, and individual studies carry their own terms.
  Review with your institution before redistributing builds or embedding them in a clinical or
  commercial service. Choose this repo's own code license separately.

## Roadmap

See [docs/DESIGN.md](docs/DESIGN.md).
