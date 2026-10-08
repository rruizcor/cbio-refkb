# Design

## Why a separate repo

One reference layer serves three consumers: the Illumina WES/WTS clinical pipeline, the ONT
somatic multi-agent pipeline, and the MolDX CGP dossier. Keeping it separate gives:

- **Locked, versioned releases** a CAP/CLIA pipeline can pin (`refkb` version + snapshot SHA).
- **Platform neutrality**: neither pipeline's conventions leak into the reference.
- **License separation**: AGPL-licensed source data stays out of pipeline repos.

## Principles

1. **Pin everything.** Every build is tied to an HF commit SHA with per-file checksums.
2. **Annotate, don't filter.** Outputs are scores and flags for triage and interpretation.
   Anything that changes a call becomes part of analytical validation.
3. **Patient-level, panel-aware counts.** Avoid multi-sample inflation and MSK-IMPACT gene bias.
4. **Platform adapters at the edge.** Each caller maps into one schema; matching tolerances are
   set per platform, not baked into the reference.

## Phases

| Phase | Scope | Status |
|---|---|---|
| 1 | SNV/indel recurrence; germline flags; normal-confirmed somatic set | **done (v0.1)** |
| 1b | Liftover GRCh37 -> GRCh38 (and T2T if needed); left-align/normalize alleles; canonical transcripts | next |
| 2 | SV/fusion priors at gene-pair + exon level, DNA vs RNA support (`structural_variant`) | planned |
| 3 | CNV expectations: arm-level and focal gain/loss frequencies, GISTIC peaks (`copy_number_seg`, `gistic`) | planned |
| 4 | Annotation CLI + platform adapters (Illumina VCF, Manta/Arriba, ONT Sniffles2/Spectre); read-only query tool for pipeline agents | planned |

## Known limitations

- Mixed genome builds until Phase 1b; the `build` column must be respected when joining.
- Germline-flagged calls are sparse and study-dependent; use `germline_profile_by_study` to see where they exist.
- Calls come from short-read panels/WES: under-represents repeat-region SVs and STRs.
- Not a substitute for an in-house Agilent panel of normals, gnomAD, or orthogonal reference materials.
