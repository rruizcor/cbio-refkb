# Phase 1 first build: cBioPortal snapshot `3e5b3b5d`

- **Snapshot:** `cBioPortal/publicDatabase@3e5b3b5de1ee243a3e19f93aff62601f88e1e22b`. The 7 Phase 1 tables are 658 MB; 542 studies have mutations.
- **Build:** `refkb build … --memory-limit 12GB --threads 4` on 4 vCPU / 16 GB RAM: **99 s, peak RSS 12.2 GB**. Outputs: 20,456,158 recurrence rows (9,076,083 pan-cancer variants, 709 tumor types plus `ALL`), 542 study profiles, 53 panel-profile rows.
- **Data:** 18,838,434 called mutation records. The 18,839,393 raw records include 959 from an uncalled profile, which is dropped.
- **Numbers:** every number here comes from `scripts/phase1_reports.py` and is in `summary_stats.json`.

> **Headline.** The hotspot frequencies are plausible: all 4 checks are within range or within 2×. The **counts are inflated about 2.3× by patients who appear in several studies**. The denominators also include ~50k clonal-hematopoiesis (blood) profiles and ~7k cell-line/PDX models. The normal-alt flag is dominated by 1–2-read noise in deep normals. Those three things should change before Phase 1b (see §9).

## 1. Germline-flagged calls
**737 germline-flagged records (0.004%) in 16 of 542 studies.** They are spelled `Germline`, `GERMLINE` or `germline`, all matched case-insensitively. Germline calls are effectively absent: cBioPortal filters them on import.

| study | type | germline records | all records |
|---|---|---:|---:|
| coadread_tcga_pub | coadread | 101 | 62,520 |
| lung_nci_2022 | luad | 100 | 8,398 |
| pog570_bcgsc_2020 | mixed | 74 | 60,112 |
| coadread_tcga | coadread | 70 | 62,797 |
| hnsc_a5consortium_2025 | hnsc | 67 | 3,674 |
| brca_tcga_pub | brca | 64 | 23,068 |
| kirc_tcga_pub | ccrcc | 57 | 22,273 |
| all_stjude_2013 | myeloid | 52 | 527 |
| ov_tcga_pub | hgsoc | 51 | 14,771 |
| prad_su2c_2019 | prad | 40 | 40,055 |

`mutationStatus` is mostly uninformative:
- `NA`: 9.75M
- `Somatic` / `SOMATIC` / `somatic`: 6.59M
- `.`: 2.21M
- `Unknown` / `__UNKNOWN__`: 0.27M

## 2. Matched-normal read counts
- **6,726,741 of 18,838,434 records (35.7%) have usable matched-normal counts**, in 259 studies (234 of them at ≥90%).
- As deposited, 6,924,355 records (36.8%) carried *some* normal value. 197,614 of those are now treated as missing because they were negative, 0/0 depth, or only one of ref/alt (§8).
- Studies with the most records carrying normal counts:

| study | records with normal counts | % |
|---|---:|---:|
| ucec_tcga_pan_can_atlas_2018 | 538,948 | 100 |
| msk_impact_50k_2026 | 478,311 | 99.9 |
| skcm_tcga_pan_can_atlas_2018 | 325,843 | 100 |
| nsclc_tcga_broad_2016 | 302,103 | 100 |
| crc_sysucc_2022 | 263,996 | 96.9 |
| msk_met_2021 | 227,931 | 100 |
| coadread_tcga_pan_can_atlas_2018 | 208,739 | 100 |
| msk_chord_2024 | 208,024 | 99.9 |

**The normal-alt flag is mostly noise.** 735,884 non-germline records have `normalAltCount > 0`, and 516,193 of them (70%) have normal VAF < 1%. Typically that's 1–2 reads in a ~600× MSK-IMPACT normal.

TP53 R175H shows "alt in normal" for 2,886 of 5,143 records with normal counts, a median of 2 reads at 621× (0.3%). That is tumor-in-normal or cross-sample noise, not germline or CH. The same applies to MSI homopolymer indels (RNF43 G659fs, MSH3 K383fs, TCF7L2 K468fs: 86–95% "normal-alt").

So `normal_alt_suspects.csv` follows the runbook's definition and has **30,347 variants**. Under a ≥3-read **and** ≥2% normal-VAF rule it would have **804**, 516 of them in CH or germline-predisposition genes (DNMT3A, TET2, ASXL1, PPM1D, JAK2, SF3B1, SRSF2, CHEK2, TP53). For example, DNMT3A R882H keeps 310 of 420 normal-alt calls.

## 3. Normal-confirmed somatic set
(Normal-confirmed = not germline-flagged, 0 alt reads, ≥10× normal.)

- **5,969,343 records (31.7%)**
- **3,538,431 distinct variants**
- 186,672 study-patients in 254 studies

"0 alt reads" is strict for deep normals: any 1-read noise call drops out (see §2).

## 4. Does "WES = all genes" hold?
- **Mutated samples with no `MUTATION_EXTENDED` panel row: 0.** The `WES` panel is explicit: 44,824 genes.
- **But 45 of the 394 WES-labelled studies are not exome-wide in practice.** They have ≥50 samples but <1,000 distinct mutated genes, and cover **27,963 of the 112,325 WES samples (25%)**. They are either targeted panels without a panel matrix or driver-only deposits:

| study | WES-labelled samples | genes with mutations | notes |
|---|---:|---:|---|
| pan_origimed_2020 | 10,194 | 450 | a 450-gene panel |
| msk_ch_2020 | 5,280 | 428 | CH, IMPACT |
| hdcn_msk_2025 | 1,270 | 215 | |
| nbl_target_2018_pub | 1,089 | 151 | |
| dlbcl_duke_2017 | 1,001 | 150 | |
| ccrcc_sjuh_2023 | 943 | 12 | |

- For genes outside their real targets, these studies inflate denominators by up to ~28k patients, so frequencies for non-panel genes are underestimated. The full list with a `looks_like_panel` flag is in `wes_labelled_studies.csv`.
- This is a heuristic: a few small, low-TMB pediatric WES studies are also caught.

## 5. Tumor-type source
- **393,118 of 414,659 samples (94.8%) have a sample-level `ONCOTREE_CODE`.** Among mutated samples it's 310,256 of 329,207 (94.2%).
- The other 21,541 fall back to the study-level type. That fallback is dominated by **`MIXED` (12,090 samples)**, which is not a tumor type. The rest (MDS 3,323, SOFT_TISSUE 2,138, SKCM 1,068, OVARY 567, PROSTATE 256) are partly non-ONCOTREE codes, so they don't merge with sample-level codes (e.g. PROSTATE vs PRAD).

## 6. GRCh37 vs GRCh38
**Records:**
- GRCh37: 15,097,597 (80.1%)
- GRCh38: 3,738,567 (19.8%)
- `NA`: 2,270 (these have GRCh37 coordinates)

**Studies:**
- GRCh37 only: 268
- GRCh38 only: 5
- Mixed: 269

Within-study mixing is a sliver: 13,772 GRCh37 records sit in 65 GRCh38-majority studies, and 2,141 GRCh38 records sit in 175 GRCh37-majority studies. **The per-record labels are trustworthy.** All 13,772, and 2,104 of the 2,141, are sites seen under the labelled build elsewhere. So liftover can be per record.

Until liftover, the same variant splits into two rows. For example, KRAS G12C in LUAD has 3,572 patients on GRCh37 and 106 on GRCh38. Pan-cancer variant rows: 5,986,284 GRCh37, 3,088,877 GRCh38, 922 NA.

## 7. Plausibility and study overlap
Counts are distinct patients over the shared (tumor type, gene) denominator, summed over all genomic encodings of the amino-acid change. "Per person" applies the de-duplication rule below.

| check | expected | build method | per person (dedup) | within 2× |
|---|---|---:|---:|---|
| KRAS G12C, LUAD | 10–15% | **11.1%** (3,727/33,643) | 12.4% (1,262/10,191) | yes |
| BRAF V600E, SKCM | 35–50% | **28.0%** (1,658/5,917) | 34.3% (818/2,386) | yes (1.25× low) |
| IDH1 R132H, ODG | 50–80% | **73.0%** (571/782) | 81.5% | yes |
| IDH1 R132H, ASTR | 50–80% | **58.1%** (389/670) | 65.0% | yes |
| IDH1 R132H, ODG+ASTR+OAST+DASTR+LGG | 50–80% | **66.4%** (1,286/1,937) | 73.1% | yes |
| EGFR L858R, LUAD | 5–15% | **7.9%** (2,656/33,643) | 9.1% (930/10,191) | yes |

**Why BRAF/SKCM is low.**
- The TCGA cohorts give 36–44%; the MSK-IMPACT cohorts give ~25%.
- `msk_ch_2023` adds 795 SKCM "patients" with 0 V600E, because its calls come from blood. Dropping them alone gives ~32%.
- Model systems add more noise (`ccle_broad_2019`: 93 SKCM cell lines).

**Study overlap is large.** I keyed persons conservatively, merging only ID families that are unambiguous across studies; everything else stays study-scoped. On that basis, **355,690 profiled study-patients are 152,543 persons (2.33×)**:

| ID family | study-patients | persons | inflation | studies |
|---|---:|---:|---:|---:|
| MSK DMP (`P-#######`) | 237,701 | 63,349 | **3.75×** | 127 |
| TCGA (12-char barcode) | 38,321 | 10,960 | **3.50×** | 128 |
| CPTAC (`C3L/C3N-`) | 1,822 | 1,276 | 1.43× | 12 |
| TARGET | 4,244 | 3,356 | 1.26× | 12 |
| other (not merged) | 73,602 | 73,602 | 1.00× | 311 |

- **Of the 78,941 persons with a cross-study ID, 62,388 appear in ≥2 studies, and up to 15.**
- By tumor type: LUAD 33,643 → 10,191 (3.3×), IDC 3.2×, PRAD 3.8×. Per-type figures are in `study_overlap_by_tumor_type.csv`.
- Further real overlaps are not merged here: ICGC, SU2C/P1000, GLASS, and CPTAC/TCGA legacy-vs-GDC re-releases.
- Frequencies move 10–25% relative (table above), but counts and confidence are overstated about 2–4×. Hotspot weights also skew toward cohorts that are re-deposited many times.

**Proposed de-duplication rule** (not applied; method unchanged in this build):
1. Add `person_key` to `samples_base`:
   - TCGA: the 12-character barcode.
   - MSK: the DMP `P-#######` ID.
   - TARGET: the barcode.
   - CPTAC: the case ID.
   - Otherwise the study-scoped `patient_unique_id`. Extend the list only with ID families curated as unambiguous.
2. Count distinct `person_key` instead of `patient_unique_id`:
   - **Numerator:** the person carries the variant in any study.
   - **Denominator:** the person was profiled for the gene in any study, using the union of their panels. The existing panel-set denominator already handles this.
3. Keep sample-level tumor type. A person with discordant types counts once in each type, as multi-type patients do today.
4. Ship it as a versioned method change (v0.2), with this report as the before/after baseline.

## 8. Code changes in this run
No column-level schema mismatches. The fetch matcher handled the `data/<table>.parquet` layout and all tests passed at the start. The changes were needed for the build to run at all and for the values in the real data:

1. **`recurrence.sql` restructured.**
   - **Before:** the original was killed by the kernel OOM killer at `--memory-limit 12GB`. In test runs it filled ~30 GB of disk while spilling. DuckDB's cap was not the problem; the memory and spill went to:
     - a `COUNT(DISTINCT)` over **225M** patient×panel-gene rows for panel denominators;
     - five `COUNT(DISTINCT patient) FILTER` aggregates;
     - `mode()`/`median()` heap state across 20M groups. That memory is outside DuckDB's cap and is never spilled.
   - **Changes:**
     - Panel denominators are now computed per *panel set* (≈0.5M rows). This is identical row-for-row on the full snapshot, and runs in 0.8 s instead of 132 s.
     - Recurrence collapses to one row per variant × tumor type × patient, then uses plain counts.
     - `mode`/`median` run only on the ~1%/10% of groups whose answer isn't fixed by min/max. The rest use fixed-size states.
   - **Verified against the original** on a 35-study subset (1.26M records): all patient-count columns are identical. The only differences are arbitrary picks that were already arbitrary before: `any_value(gene)` at 158 two-gene loci (which shifts `n_patients_profiled` in 23 rows) and `mode` ties.
2. **Two documented semantic changes** in that rewrite:
   - `median_tumor_vaf` is now the median over patients of each patient's mean VAF; before, it was over records, so multi-sample patients counted more than once.
   - The annotation mode is counted over patients.
3. **New `mutations.sql` view**, which all builds read:
   - Drops `*_mutations_uncalled` (959 records in `glioma_msk_2018`; these are read counts at other samples' sites, not calls).
   - Maps chr `23`/`24` → `X`/`Y` (722,187 records; the snapshot has no `X` strings) and strips `chr`.
   - Normalizes build aliases (`37` → GRCh37).
   - Sets negative read counts to NULL (1,237 records; 1,147 have negative normal counts, 975 of them in `uvm_tcga`, plus `angs_painter_2025`'s −2/−2 placeholders).
   - Counts a normal as present only if ref and alt are both ≥0 and depth > 0. Before, alt-only normals (e.g. `coad_silu_2022`, `mel_ucla_2016`) counted as normal-alt but not as "with normal counts", so the suspect ratio could exceed 1. Ref-only and 0/0 normals (`blca_tcga_pan_can_atlas_2018` 99,709; `dlbcl_dfci_2018` all 0/0) diluted the denominator.
4. **Tests:** 7 → 12. A new `quirks_snapshot` fixture reproduces the snapshot's encodings: 23/chrX, the uncalled profile, negative, 0/0 and alt-only normals, build `37`, a multi-sample patient, and a patient on two panels. 4 of the new tests fail on the original code. The fifth (two-panel patient) guards the denominator rewrite.
5. **`scripts/phase1_reports.py`** regenerates every file in this folder.

The extra report files are `plausibility_checks.csv`, `study_overlap_by_*.csv`, `wes_labelled_studies.csv` and `summary_stats.json`. No CSV needed trimming; the largest is `normal_alt_suspects.csv` at 4.6 MB. That file has one row per variant, using pan-cancer (`ALL`) evidence.

## 9. Next steps for Phase 1b
**Change first.** These affect every frequency, so fix them before spending effort on liftover:
1. **Study selection.**
   - Exclude the CH studies (`msk_ch_2020`, `msk_ch_2023`, `msk_ch_ped_2021`: 50,587 profiled patients, 14% of the denominator). They could become a separate CH prior, which would help tumor-only triage.
   - Exclude the 18 cell-line/PDX/organoid studies (5.43M records, 29% of all records; `pancan_pdmr_2025` alone is 3.3M).
2. **Person-level de-duplication** using the rule in §7.
3. **Normal-read thresholds.**
   - `normal_alt`: require ≥3 reads and ≥2% normal VAF (or a binomial test against a tumor-in-normal rate).
   - `normal_confirmed`: use normal VAF <1% at depth ≥10, not 0 reads.
4. **WES-labelled panels.** For the 45 flagged studies, either attach the real gene list (OrigiMed 450, CH-IMPACT, …), restrict their denominators to genes observed in that study, or drop them from off-target denominators.
5. **Tumor types.** Map study-level fallbacks to ONCOTREE codes, and exclude `MIXED` from per-type rows (it stays in `ALL`).

**Phase 1b itself.**
- Lift over per record: the labels are verified, and `NA` records are GRCh37.
- Normalize alleles (left-align, trim) against the target FASTA, and decide how to handle MNVs. Today KRAS G12C appears as an SNV and three DNV encodings in separate rows.
- Re-key recurrence on normalized GRCh38 alleles so the build rows merge.
- Re-run `scripts/phase1_reports.py` on every snapshot and keep the plausibility table as a regression gate.
