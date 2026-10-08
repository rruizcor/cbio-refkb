"""Synthetic-fixture tests for the Phase 1 build logic.

Cohort (`snapshot`):
  study s_wes   (type_of_cancer luad, WES):    patients P1, P2 (P2 has two samples)
  study s_panel (mixed, panel IMPACT_MINI = KRAS, TP53):
        P3 (ONCOTREE LUAD), P4 (PAAD), P5 (LUAD)

Cohort (`quirks_snapshot`): encodings seen in the real snapshot.
  study s_quirk (skcm, WES): patients Q1, Q2, Q3 (two samples), Q4
  study s_multi (coad): M1 sequenced on PANEL_A and PANEL_B, M2 on PANEL_A only
"""

import json

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from refkb.build import build

KRAS = dict(chr="12", startPosition=25398285, endPosition=25398285, referenceAllele="C",
            tumorSeqAllele="A", gene="KRAS", entrez=3845, proteinChange="G12C")
TP53 = dict(chr="17", startPosition=7577120, endPosition=7577120, referenceAllele="C",
            tumorSeqAllele="T", gene="TP53", entrez=7157, proteinChange="R273H")
EGFR = dict(chr="7", startPosition=55259515, endPosition=55259515, referenceAllele="T",
            tumorSeqAllele="G", gene="EGFR", entrez=1956, proteinChange="L858R")
BRAF = dict(chr="7", startPosition=140453136, endPosition=140453136, referenceAllele="A",
            tumorSeqAllele="T", gene="BRAF", entrez=673, proteinChange="V600E")
KDM6A = dict(chr="23", startPosition=44941819, endPosition=44941819, referenceAllele="C",
             tumorSeqAllele="T", gene="KDM6A", entrez=7403, proteinChange="R1111C")
PIK3CA = dict(chr="3", startPosition=178936091, endPosition=178936091, referenceAllele="G",
              tumorSeqAllele="A", gene="PIK3CA", entrez=5290, proteinChange="E545K")
APC = dict(chr="5", startPosition=112175639, endPosition=112175639, referenceAllele="C",
           tumorSeqAllele="T", gene="APC", entrez=324, proteinChange="R1450*")


def _write(snap, table, rows):
    d = snap / table
    d.mkdir(parents=True)
    pq.write_table(pa.Table.from_pylist(rows), d / "part-00000.parquet")


def _mut(study, sample, patient, v, status="Somatic", t_alt=20, t_ref=80, n_alt=None, n_ref=None, **overrides):
    return {
        "molecularProfileId": f"{study}_mutations", "sampleId": sample, "sampleInternalId": 0,
        "patientId": patient, "entrezGeneId": v["entrez"], "studyId": study,
        "mutationStatus": status, "tumorAltCount": t_alt, "tumorRefCount": t_ref,
        "normalAltCount": n_alt, "normalRefCount": n_ref, "chr": v["chr"],
        "startPosition": v["startPosition"], "endPosition": v["endPosition"],
        "referenceAllele": v["referenceAllele"], "tumorSeqAllele": v["tumorSeqAllele"],
        "proteinChange": v["proteinChange"], "mutationType": "Missense_Mutation",
        "ncbiBuild": "GRCh37", "variantType": "SNP", "GENE.hugoGeneSymbol": v["gene"],
    } | overrides


@pytest.fixture
def snapshot(tmp_path):
    snap = tmp_path / "raw" / "abc123"
    samples = [
        ("s_wes", "S1", "P1"), ("s_wes", "S2a", "P2"), ("s_wes", "S2b", "P2"),
        ("s_panel", "S3", "P3"), ("s_panel", "S4", "P4"), ("s_panel", "S5", "P5"),
    ]
    _write(snap, "sample_derived", [
        {"sample_unique_id": f"{st}_{s}", "sample_stable_id": s, "patient_unique_id": f"{st}_{p}",
         "cancer_study_identifier": st, "sample_type": "Primary Solid Tumor"} for st, s, p in samples
    ])
    _write(snap, "cancer_study", [
        {"cancer_study_identifier": "s_wes", "type_of_cancer_id": "luad", "name": "WES"},
        {"cancer_study_identifier": "s_panel", "type_of_cancer_id": "mixed", "name": "Panel"},
    ])
    _write(snap, "type_of_cancer", [{"type_of_cancer_id": "luad", "name": "Lung adenocarcinoma"}])
    _write(snap, "clinical_data_derived", [
        {"sample_unique_id": "s_panel_S3", "attribute_name": "ONCOTREE_CODE", "attribute_value": "LUAD", "type": "sample"},
        {"sample_unique_id": "s_panel_S4", "attribute_name": "ONCOTREE_CODE", "attribute_value": "PAAD", "type": "sample"},
        {"sample_unique_id": "s_panel_S5", "attribute_name": "ONCOTREE_CODE", "attribute_value": "LUAD", "type": "sample"},
    ])
    _write(snap, "sample_to_gene_panel_derived", [
        {"sample_unique_id": f"{st}_{s}", "alteration_type": "MUTATION_EXTENDED",
         "gene_panel_id": "WES" if st == "s_wes" else "IMPACT_MINI", "cancer_study_identifier": st}
        for st, s, _ in samples
    ])
    _write(snap, "gene_panel_to_gene_derived", [
        {"gene_panel_id": "IMPACT_MINI", "gene": "KRAS"}, {"gene_panel_id": "IMPACT_MINI", "gene": "TP53"},
    ])
    _write(snap, "mutation_derived", [
        _mut("s_wes", "S1", "P1", KRAS, n_alt=0, n_ref=30),
        _mut("s_wes", "S2a", "P2", KRAS),
        _mut("s_wes", "S2b", "P2", KRAS),          # same patient, second sample: counted once
        _mut("s_panel", "S3", "P3", KRAS),
        _mut("s_panel", "S4", "P4", KRAS, n_alt=2, n_ref=40),   # alt reads in normal
        _mut("s_panel", "S5", "P5", TP53, status="Germline", t_alt=50, t_ref=50),
        _mut("s_wes", "S1", "P1", EGFR),
    ])
    (snap / "fetch_manifest.json").write_text(json.dumps({"revision_sha": "abc123", "hf_repo": "test"}))
    return snap


@pytest.fixture
def quirks_snapshot(tmp_path):
    snap = tmp_path / "raw" / "quirks"
    samples = [
        ("s_quirk", "Q1-T", "Q1", "WES"), ("s_quirk", "Q2-T", "Q2", "WES"),
        ("s_quirk", "Q3-T1", "Q3", "WES"), ("s_quirk", "Q3-T2", "Q3", "WES"), ("s_quirk", "Q4-T", "Q4", "WES"),
        ("s_multi", "M1-A", "M1", "PANEL_A"), ("s_multi", "M1-B", "M1", "PANEL_B"), ("s_multi", "M2-A", "M2", "PANEL_A"),
    ]
    _write(snap, "sample_derived", [
        {"sample_unique_id": f"{st}_{s}", "sample_stable_id": s, "patient_unique_id": f"{st}_{p}",
         "patient_stable_id": p, "cancer_study_identifier": st, "sample_type": "Primary Solid Tumor"}
        for st, s, p, _ in samples
    ])
    _write(snap, "cancer_study", [
        {"cancer_study_identifier": "s_quirk", "type_of_cancer_id": "skcm", "name": "Quirks"},
        {"cancer_study_identifier": "s_multi", "type_of_cancer_id": "coad", "name": "Two panels"},
    ])
    _write(snap, "type_of_cancer", [{"type_of_cancer_id": "skcm", "name": "Melanoma"}])
    _write(snap, "clinical_data_derived", [
        {"sample_unique_id": "", "patient_unique_id": "s_quirk_Q1", "attribute_name": "SEX",
         "attribute_value": "F", "type": "patient"},
    ])
    _write(snap, "sample_to_gene_panel_derived", [
        {"sample_unique_id": f"{st}_{s}", "alteration_type": "MUTATION_EXTENDED", "gene_panel_id": panel,
         "cancer_study_identifier": st} for st, s, _, panel in samples
    ])
    _write(snap, "gene_panel_to_gene_derived", [
        {"gene_panel_id": "PANEL_A", "gene": "KRAS"}, {"gene_panel_id": "PANEL_A", "gene": "APC"},
        {"gene_panel_id": "PANEL_B", "gene": "KRAS"}, {"gene_panel_id": "PANEL_B", "gene": "TP53"},
    ])
    q = "s_quirk"
    _write(snap, "mutation_derived", [
        # BRAF: VAFs 0.1, 0.3, Q3 = mean(0.5, 0.7) = 0.6, Q4 unusable (negative tumor count)
        _mut(q, "Q1-T", "Q1", BRAF, t_alt=10, t_ref=90, n_alt=0, n_ref=30),
        _mut(q, "Q2-T", "Q2", BRAF, t_alt=30, t_ref=70, n_alt=-2, n_ref=-2),         # deposit error
        _mut(q, "Q3-T1", "Q3", BRAF, t_alt=50, t_ref=50, n_alt=0, n_ref=0,           # zero-depth normal
             proteinChange="V640E"),
        _mut(q, "Q3-T2", "Q3", BRAF, t_alt=70, t_ref=30, proteinChange="V640E"),
        _mut(q, "Q4-T", "Q4", BRAF, t_alt=-777, t_ref=80, n_alt=3, n_ref=None),     # alt-only normal
        # KDM6A on X, deposited as 23, chrX and with build alias '37'
        _mut(q, "Q1-T", "Q1", KDM6A, t_alt=20, t_ref=80),
        _mut(q, "Q2-T", "Q2", KDM6A, t_alt=40, t_ref=60, chr="chrX", ncbiBuild="37"),
        # PIK3CA: one call, plus a record from the uncalled profile
        _mut(q, "Q1-T", "Q1", PIK3CA),
        _mut(q, "Q2-T", "Q2", PIK3CA, status="UNCALLED", molecularProfileId=f"{q}_mutations_uncalled"),
        _mut("s_multi", "M1-B", "M1", TP53),
        _mut("s_multi", "M2-A", "M2", KRAS),
    ])
    (snap / "fetch_manifest.json").write_text(json.dumps({"revision_sha": "quirks", "hf_repo": "test"}))
    return snap


def _row(out, tumor_type, gene):
    con = duckdb.connect()
    return con.execute(
        f"SELECT * FROM read_parquet('{out / 'variant_recurrence.parquet'}') WHERE tumor_type=? AND gene=?",
        [tumor_type, gene],
    ).fetchdf().iloc[0]


def test_panel_aware_denominator_and_patient_dedup(snapshot, tmp_path):
    out = build(snapshot, tmp_path / "builds")
    r = _row(out, "LUAD", "KRAS")
    assert r.n_patients == 3            # P1, P2 (deduped), P3
    assert r.n_patients_profiled == 4   # WES P1,P2 + panel P3,P5
    assert r.patient_frequency == pytest.approx(0.75)


def test_off_panel_gene_uses_wes_only(snapshot, tmp_path):
    out = build(snapshot, tmp_path / "builds")
    r = _row(out, "LUAD", "EGFR")
    assert r.n_patients_profiled == 2
    assert r.patient_frequency == pytest.approx(0.5)


def test_pan_cancer_row(snapshot, tmp_path):
    out = build(snapshot, tmp_path / "builds")
    r = _row(out, "ALL", "KRAS")
    assert r.n_patients == 4
    assert r.n_patients_profiled == 5


def test_normal_count_and_germline_flags(snapshot, tmp_path):
    out = build(snapshot, tmp_path / "builds")
    luad = _row(out, "LUAD", "KRAS")
    assert luad.n_patients_normal_confirmed == 1
    assert luad.n_patients_with_normal_counts == 1
    paad = _row(out, "PAAD", "KRAS")
    assert paad.n_patients_normal_alt == 1
    assert paad.n_patients_normal_confirmed == 0
    tp53 = _row(out, "LUAD", "TP53")
    assert tp53.n_patients == 0
    assert tp53.n_patients_germline_flag == 1


def test_profiles_and_manifest(snapshot, tmp_path):
    out = build(snapshot, tmp_path / "builds")
    con = duckdb.connect()
    prof = con.execute(f"SELECT * FROM read_parquet('{out / 'germline_profile_by_study.parquet'}')").fetchdf()
    panel = prof.set_index("study_id").loc["s_panel"]
    assert panel.n_germline == 1 and panel.n_somatic == 2
    m = json.loads((out / "build_manifest.json").read_text())
    assert m["snapshot_id"] == "abc123"
    assert m["outputs"]["variant_recurrence"]["rows"] > 0


def test_memory_limit_option(snapshot, tmp_path):
    out = build(snapshot, tmp_path / "builds", memory_limit="1GB", threads=2)
    assert (out / "variant_recurrence.parquet").exists()


def _rows(out, tumor_type, gene):
    return duckdb.connect().execute(
        f"SELECT * FROM read_parquet('{out / 'variant_recurrence.parquet'}') WHERE tumor_type=? AND gene=?",
        [tumor_type, gene],
    ).fetchdf()


def test_chrom_and_build_aliases_merge(quirks_snapshot, tmp_path):
    out = build(quirks_snapshot, tmp_path / "builds")
    kdm6a = _rows(out, "SKCM", "KDM6A")
    assert len(kdm6a) == 1                       # '23' + 'chrX', 'GRCh37' + '37' -> one variant
    r = kdm6a.iloc[0]
    assert (r.chrom, r.build, r.n_patients) == ("X", "GRCh37", 2)
    assert r.median_tumor_vaf == pytest.approx(0.3)


def test_uncalled_profile_excluded(quirks_snapshot, tmp_path):
    out = build(quirks_snapshot, tmp_path / "builds")
    assert _row(out, "SKCM", "PIK3CA").n_patients == 1
    prof = duckdb.connect().execute(
        f"SELECT * FROM read_parquet('{out / 'germline_profile_by_study.parquet'}')").fetchdf().set_index("study_id")
    assert prof.loc["s_quirk"].n_mutation_records == 8
    assert prof.loc["s_quirk"].ncbi_builds == "GRCh37"


def test_unusable_read_counts_are_missing(quirks_snapshot, tmp_path):
    out = build(quirks_snapshot, tmp_path / "builds")
    r = _row(out, "SKCM", "BRAF")
    assert r.n_patients == 4
    assert r.n_patients_with_normal_counts == 1  # only Q1: negative, 0/0 and alt-only normals carry no evidence
    assert r.n_patients_normal_confirmed == 1
    assert r.n_patients_normal_alt == 0
    prof = duckdb.connect().execute(
        f"SELECT * FROM read_parquet('{out / 'germline_profile_by_study.parquet'}')").fetchdf().set_index("study_id")
    assert prof.loc["s_quirk"].pct_with_normal_counts == pytest.approx(12.5)


def test_patient_level_vaf_and_annotation(quirks_snapshot, tmp_path):
    out = build(quirks_snapshot, tmp_path / "builds")
    r = _row(out, "SKCM", "BRAF")
    assert r.median_tumor_vaf == pytest.approx(0.3)   # median of per-patient VAFs 0.1, 0.3, 0.6
    assert r.protein_change == "V600E"                # 3 patients vs 1 (Q3 has two V640E records)


def test_patient_on_two_panels_counted_once(quirks_snapshot, tmp_path):
    out = build(quirks_snapshot, tmp_path / "builds")
    assert _row(out, "COAD", "TP53").n_patients_profiled == 1   # PANEL_B only: M1
    kras = _row(out, "COAD", "KRAS")                             # both panels: M1 once + M2
    assert kras.n_patients_profiled == 2
    assert kras.patient_frequency == pytest.approx(0.5)
