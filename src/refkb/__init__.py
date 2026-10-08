"""cbio-refkb: versioned cancer-cohort reference layer built from public cBioPortal."""

__version__ = "0.1.0"

HF_REPO = "cBioPortal/publicDatabase"

# Tables needed for Phase 1 (SNV/indel recurrence + germline/normal-count profiling).
PHASE1_TABLES = [
    "mutation_derived",
    "sample_derived",
    "cancer_study",
    "type_of_cancer",
    "clinical_data_derived",
    "sample_to_gene_panel_derived",
    "gene_panel_to_gene_derived",
]
