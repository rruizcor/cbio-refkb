# Phase 1 cloud run (instructions for Claude Code)

Run the Phase 1 build on the real cBioPortal snapshot, check the results make sense,
and open a pull request with the reports. Work without asking unless a step says to stop.

## 0. Preflight
- Check network: `python -c "from huggingface_hub import HfApi; print(HfApi().dataset_info('cBioPortal/publicDatabase').sha)"`.
  If it fails with a proxy/connection error, **stop** and report the exact host(s) refused,
  so they can be added to the environment's Custom allowlist.
- Record resources with `nproc`, `free -g`, `df -h .`.

## 1. Install and test
- `pip install -e ".[dev]"` (add `--break-system-packages` if pip requires it), then `pytest -q`.
- If tests fail, fix the cause before going on.

## 2. Fetch
- `refkb fetch --out data/raw` and note the snapshot SHA it prints.
- If a table is reported missing or the parquet layout differs from what `fetch.table_files`
  expects, inspect `HfApi().list_repo_files(...)`, fix the matcher, add a test case, continue.

## 3. Build
- `refkb build data/raw/<sha> --out builds --memory-limit 12GB --threads 4`
- If a column referenced in `src/refkb/sql/` does not exist in the real data, compare with
  `schema.sql` in the snapshot, fix the SQL, update the synthetic fixtures in
  `tests/test_build.py` to match the real schema, and re-run tests.

## 4. Reports (commit these, never the raw data or full builds)
Write to `reports/<sha>/`:
- `germline_profile_by_study.csv` (full)
- `panel_profile.csv` (full)
- `top_recurrent_pan_cancer.csv`: tumor_type = 'ALL', n_patients >= 10, top 500 by patient_frequency
- `top_recurrent_by_tumor_type.csv`: top 25 per tumor type with n_patients_profiled >= 200
- `normal_alt_suspects.csv`: variants with n_patients_with_normal_counts >= 5 where
  n_patients_normal_alt / n_patients_with_normal_counts >= 0.3 (germline leak / CHIP candidates)
- `build_manifest.json` (copy)
Keep every CSV under 5 MB; trim row counts if needed and say so.

## 5. Plausibility checks
Look up patient_frequency for these and flag any that differ from the expected
range by more than about 2x; if one does, investigate (tumor-type mapping, denominators,
duplicate studies) before writing the summary:
- KRAS G12C in LUAD: roughly 10–15%
- BRAF V600E in SKCM: roughly 35–50%
- IDH1 R132H in low-grade glioma codes (e.g. ODG, ASTR, or study-level LGG): roughly 50–80%
- EGFR L858R in LUAD: roughly 5–15%
Also check whether the same patients appear under several overlapping studies
(e.g. TCGA PanCanAtlas vs legacy TCGA vs per-paper cohorts). If they do, report how much
this inflates counts and propose a de-duplication rule; do not silently change the method.

## 6. SUMMARY.md
Write `reports/<sha>/SUMMARY.md` (concise, numbers first) answering:
1. Which studies carry germline-flagged calls (top 10 by count) and the total germline-flagged records.
2. Share of all mutation records with matched-normal read counts; studies with the most.
3. Size of the normal-confirmed somatic set (records and distinct variants).
4. Does the "WES = all genes" denominator assumption hold? Mutated samples with no panel row?
5. Share of samples with sample-level ONCOTREE_CODE vs study-level tumor type.
6. GRCh37 vs GRCh38 mix.
7. Plausibility check results and any study-overlap inflation found.
8. Code changes made during the run and why.
9. Recommended next steps for Phase 1b (liftover/normalization) and anything that should change first.

## 7. Deliver
- Branch `phase1-first-build`; commit code fixes and `reports/` (data/ and builds/ stay ignored).
- Push and open a pull request whose description is the SUMMARY.md content.
