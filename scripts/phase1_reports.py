"""Write the Phase 1 review reports (docs/CLOUD_RUN.md step 4-5) for one build.

    python scripts/phase1_reports.py data/raw/<sha> builds/<sha> reports/<sha>

Outputs small CSVs that can be committed (never the raw data or full builds),
plus summary_stats.json with every number quoted in SUMMARY.md.
"""

from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path

from refkb.build import connect

MAX_CSV_BYTES = 5 * 1024 * 1024

# Cross-study person identity, used only to measure study-overlap inflation; the
# build itself still counts study-scoped patient_unique_id. Only ID families that
# are unambiguous across studies are merged; everything else stays study-scoped.
PERSON_KEY_SQL = """
CREATE OR REPLACE TEMP VIEW person AS
SELECT
    sample_unique_id,
    CASE
        WHEN regexp_matches(patient_stable_id, '^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}') THEN 'TCGA'
        WHEN regexp_matches(patient_stable_id, '^P-[0-9]{7}$')                THEN 'MSK-DMP'
        WHEN regexp_matches(patient_stable_id, '^TARGET-[0-9]{2}-[A-Z0-9]+$') THEN 'TARGET'
        WHEN regexp_matches(patient_stable_id, '^C3[LN]-[0-9]{5}$')           THEN 'CPTAC'
        ELSE 'study-scoped'
    END AS id_family,
    CASE
        WHEN regexp_matches(patient_stable_id, '^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}') THEN left(patient_stable_id, 12)
        WHEN regexp_matches(patient_stable_id, '^(P-[0-9]{7}|TARGET-[0-9]{2}-[A-Z0-9]+|C3[LN]-[0-9]{5})$')
             THEN patient_stable_id
        ELSE patient_unique_id
    END AS person_key
FROM sample_derived
"""

# (label, gene, protein_change, tumor types, expected low, expected high)
PLAUSIBILITY = [
    ("KRAS G12C in LUAD", "KRAS", "G12C", ["LUAD"], 0.10, 0.15),
    ("BRAF V600E in SKCM", "BRAF", "V600E", ["SKCM"], 0.35, 0.50),
    ("IDH1 R132H in ODG", "IDH1", "R132H", ["ODG"], 0.50, 0.80),
    ("IDH1 R132H in ASTR", "IDH1", "R132H", ["ASTR"], 0.50, 0.80),
    ("IDH1 R132H in grade 2 diffuse glioma (ODG, ASTR, OAST, DASTR, LGG)", "IDH1", "R132H",
     ["ODG", "ASTR", "OAST", "DASTR", "LGG"], 0.50, 0.80),
    ("EGFR L858R in LUAD", "EGFR", "L858R", ["LUAD"], 0.05, 0.15),
]


def copy_csv(con, query: str, dest: Path) -> int:
    con.execute(f"COPY ({query}) TO '{dest.as_posix()}' (HEADER, DELIMITER ',')")
    size = dest.stat().st_size
    if size > MAX_CSV_BYTES:
        raise SystemExit(f"{dest.name} is {size / 1e6:.1f} MB, over the 5 MB limit: trim it")
    return con.execute(f"SELECT count(*) FROM read_csv_auto('{dest.as_posix()}')").fetchone()[0]


def frequency(con, gene: str, protein_change: str, tumor_types: list[str], key: str) -> tuple[int, int]:
    """(n carriers, n profiled for the gene) counting distinct `key` (patient or person)."""
    num = con.execute(f"""
        SELECT count(DISTINCT {key})
        FROM mutations m
        JOIN samples_x sx ON sx.study_id = m.studyId AND sx.sample_stable_id = m.sampleId
        JOIN person USING (sample_unique_id)
        WHERE sx.tumor_type IN (SELECT unnest(?)) AND m."GENE.hugoGeneSymbol" = ? AND m.proteinChange = ?
          AND lower(coalesce(m.mutationStatus, '')) <> 'germline'
    """, [tumor_types, gene, protein_change]).fetchone()[0]
    den = con.execute(f"""
        SELECT count(DISTINCT {key})
        FROM sample_to_gene_panel_derived p
        JOIN samples_x sx USING (sample_unique_id)
        JOIN person USING (sample_unique_id)
        WHERE p.alteration_type = 'MUTATION_EXTENDED' AND sx.tumor_type IN (SELECT unnest(?))
          AND (p.gene_panel_id = 'WES'
               OR p.gene_panel_id IN (SELECT gene_panel_id FROM gene_panel_to_gene_derived WHERE gene = ?))
    """, [tumor_types, gene]).fetchone()[0]
    return num, den


def main(snapshot_dir: str, build_dir: str, out_dir: str) -> None:
    build_dir, out = Path(build_dir), Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    con = connect(Path(snapshot_dir))
    con.execute("SET threads = 4")
    con.execute("SET memory_limit = '8GB'")
    con.execute(PERSON_KEY_SQL)
    vr = f"read_parquet('{(build_dir / 'variant_recurrence.parquet').as_posix()}')"
    order = "patient_frequency DESC, n_patients DESC, gene, build, chrom, start_pos, ref, alt"
    rows = {}

    # Step 4: reports
    rows["germline_profile_by_study.csv"] = copy_csv(
        con, f"SELECT * FROM read_parquet('{(build_dir / 'germline_profile_by_study.parquet').as_posix()}')",
        out / "germline_profile_by_study.csv")
    rows["panel_profile.csv"] = copy_csv(
        con, f"SELECT * FROM read_parquet('{(build_dir / 'panel_profile.parquet').as_posix()}')",
        out / "panel_profile.csv")
    rows["top_recurrent_pan_cancer.csv"] = copy_csv(
        con, f"SELECT * FROM {vr} WHERE tumor_type = 'ALL' AND n_patients >= 10 ORDER BY {order} LIMIT 500",
        out / "top_recurrent_pan_cancer.csv")
    rows["top_recurrent_by_tumor_type.csv"] = copy_csv(con, f"""
        SELECT * EXCLUDE (rk) FROM (
            SELECT *, row_number() OVER (PARTITION BY tumor_type ORDER BY {order}) AS rk
            FROM {vr} WHERE tumor_type <> 'ALL' AND n_patients_profiled >= 200 AND n_patients > 0)
        WHERE rk <= 25 ORDER BY tumor_type, rk""", out / "top_recurrent_by_tumor_type.csv")
    # One row per variant: pan-cancer ('ALL') evidence.
    rows["normal_alt_suspects.csv"] = copy_csv(con, f"""
        SELECT *, n_patients_normal_alt / n_patients_with_normal_counts AS frac_normal_alt
        FROM {vr}
        WHERE tumor_type = 'ALL' AND n_patients_with_normal_counts >= 5
          AND n_patients_normal_alt / n_patients_with_normal_counts >= 0.3
        ORDER BY n_patients_normal_alt DESC, frac_normal_alt DESC, gene, build, chrom, start_pos, ref, alt""",
        out / "normal_alt_suspects.csv")
    shutil.copy(build_dir / "build_manifest.json", out / "build_manifest.json")

    # Step 5: plausibility, per study-scoped patient (the build's method) and per person.
    con.execute("CREATE OR REPLACE TEMP TABLE plaus (check_name VARCHAR, tumor_types VARCHAR, expected VARCHAR, "
                "n_patients INT, n_profiled INT, frequency DOUBLE, n_persons INT, n_persons_profiled INT, "
                "frequency_person_dedup DOUBLE, within_2x BOOLEAN)")
    for label, gene, pc, tts, lo, hi in PLAUSIBILITY:
        n, d = frequency(con, gene, pc, tts, "sx.patient_unique_id")
        pn, pd = frequency(con, gene, pc, tts, "person.person_key")
        f = n / d if d else None
        con.execute("INSERT INTO plaus VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", [
            label, "+".join(tts), f"{lo:.0%}-{hi:.0%}", n, d, f, pn, pd, pn / pd if pd else None,
            f is not None and lo / 2 <= f <= hi * 2])
    rows["plausibility_checks.csv"] = copy_csv(con, "SELECT * FROM plaus", out / "plausibility_checks.csv")

    # Study overlap: profiled study-patients vs distinct persons, by ID family and by tumor type.
    con.execute("""
        CREATE OR REPLACE TEMP TABLE prof AS
        SELECT DISTINCT sx.tumor_type, sx.patient_unique_id, sx.study_id, person.id_family, person.person_key
        FROM sample_to_gene_panel_derived p
        JOIN samples_x sx USING (sample_unique_id)
        JOIN person USING (sample_unique_id)
        WHERE p.alteration_type = 'MUTATION_EXTENDED'""")
    rows["study_overlap_by_id_family.csv"] = copy_csv(con, """
        SELECT id_family, count(DISTINCT patient_unique_id) AS n_study_patients, count(DISTINCT person_key) AS n_persons,
               round(count(DISTINCT patient_unique_id) / count(DISTINCT person_key), 2) AS inflation,
               count(DISTINCT study_id) AS n_studies
        FROM prof WHERE tumor_type = 'ALL' GROUP BY id_family
        UNION ALL
        SELECT 'total', count(DISTINCT patient_unique_id), count(DISTINCT person_key),
               round(count(DISTINCT patient_unique_id) / count(DISTINCT person_key), 2), count(DISTINCT study_id)
        FROM prof WHERE tumor_type = 'ALL'
        ORDER BY n_study_patients DESC""", out / "study_overlap_by_id_family.csv")
    rows["study_overlap_by_tumor_type.csv"] = copy_csv(con, """
        SELECT tumor_type, count(DISTINCT patient_unique_id) AS n_study_patients, count(DISTINCT person_key) AS n_persons,
               round(count(DISTINCT patient_unique_id) / count(DISTINCT person_key), 2) AS inflation,
               count(DISTINCT study_id) AS n_studies
        FROM prof GROUP BY tumor_type HAVING count(DISTINCT patient_unique_id) >= 200
        ORDER BY n_study_patients DESC""", out / "study_overlap_by_tumor_type.csv")

    # Numbers for SUMMARY.md
    s = {"rows_written": rows}
    one = lambda q, *a: con.execute(q, list(a)).fetchone()
    gp = f"read_parquet('{(build_dir / 'germline_profile_by_study.parquet').as_posix()}')"
    s["n_mutation_records_raw"] = one("SELECT count(*) FROM mutation_derived")[0]
    s["n_mutation_records_called"] = one("SELECT count(*) FROM mutations")[0]
    s["n_studies_with_mutations"] = one("SELECT count(DISTINCT studyId) FROM mutations")[0]
    s["germline_total_records"] = one(f"SELECT sum(n_germline) FROM {gp}")[0]
    s["germline_n_studies"] = one(f"SELECT count(*) FROM {gp} WHERE n_germline > 0")[0]
    s["germline_top10"] = con.execute(f"""SELECT study_id, study_cancer_type, n_germline, n_mutation_records
        FROM {gp} WHERE n_germline > 0 ORDER BY n_germline DESC, study_id LIMIT 10""").fetchall()
    s["status_values"] = con.execute("""SELECT coalesce(mutationStatus, 'NULL'), count(*) FROM mutations
        GROUP BY 1 ORDER BY 2 DESC LIMIT 12""").fetchall()
    # Usable normal counts (mutations view) vs any normal value present as deposited.
    usable = one("SELECT count(*) FROM mutations WHERE normalRefCount IS NOT NULL")[0]
    raw_ref, raw_any = one("""
        SELECT count(*) FILTER (WHERE normalRefCount IS NOT NULL),
               count(*) FILTER (WHERE normalRefCount IS NOT NULL OR normalAltCount IS NOT NULL)
        FROM mutation_derived WHERE molecularProfileId NOT LIKE '%\\_mutations\\_uncalled' ESCAPE '\\'""")
    s["normal_counts"] = {"records_with_usable_normal_counts": usable,
                          "records_with_normal_ref_as_deposited": raw_ref,
                          "records_with_any_normal_value_as_deposited": raw_any}
    s["normal_counts"]["n_studies_with_any"] = one(f"SELECT count(*) FROM {gp} WHERE pct_with_normal_counts > 0")[0]
    s["normal_counts"]["n_studies_ge_90pct"] = one(f"SELECT count(*) FROM {gp} WHERE pct_with_normal_counts >= 90")[0]
    s["normal_counts_top10"] = con.execute(f"""SELECT study_id, round(n_mutation_records * pct_with_normal_counts / 100)::BIGINT,
        pct_with_normal_counts, n_mutation_records FROM {gp} ORDER BY 2 DESC LIMIT 10""").fetchall()
    s["normal_confirmed"] = dict(zip(["records", "distinct_variants", "patients", "studies"], one(f"""
        SELECT count(*), count(DISTINCT (ncbiBuild, chr, startPosition, endPosition, referenceAllele, tumorSeqAllele)),
               count(DISTINCT (studyId, patientId)), count(DISTINCT studyId)
        FROM mutations
        WHERE lower(coalesce(mutationStatus, '')) <> 'germline' AND normalAltCount = 0
          AND normalRefCount >= {json.loads((build_dir / 'build_manifest.json').read_text())['params']['min_normal_depth']}""")))
    s["normal_alt_any"] = one("""SELECT count(*) FROM mutations
        WHERE lower(coalesce(mutationStatus, '')) <> 'germline' AND normalAltCount > 0""")[0]
    s["panel_profile"] = con.execute(
        f"SELECT * FROM read_parquet('{(build_dir / 'panel_profile.parquet').as_posix()}')").fetchall()
    s["wes_gene_list_size"] = one("SELECT count(DISTINCT gene) FROM gene_panel_to_gene_derived WHERE gene_panel_id = 'WES'")[0]
    # WES-labelled studies whose deposited calls span < 1000 genes look like panels / driver-only deposits.
    con.execute("""
        CREATE OR REPLACE TEMP TABLE wes_studies AS
        WITH w AS (SELECT cancer_study_identifier AS study_id, count(DISTINCT sample_unique_id) AS n_wes_samples
                   FROM sample_to_gene_panel_derived WHERE alteration_type = 'MUTATION_EXTENDED' AND gene_panel_id = 'WES'
                   GROUP BY 1),
             g AS (SELECT m.studyId AS study_id, count(DISTINCT m."GENE.hugoGeneSymbol") AS n_genes_mutated,
                          count(*) / count(DISTINCT m.sampleId) AS muts_per_mutated_sample
                   FROM mutations m JOIN samples_base sb ON sb.study_id = m.studyId AND sb.sample_stable_id = m.sampleId
                   JOIN sample_to_gene_panel_derived p ON p.sample_unique_id = sb.sample_unique_id
                   WHERE p.alteration_type = 'MUTATION_EXTENDED' AND p.gene_panel_id = 'WES' GROUP BY 1)
        SELECT w.*, coalesce(g.n_genes_mutated, 0) AS n_genes_mutated, g.muts_per_mutated_sample
        FROM w LEFT JOIN g USING (study_id)""")
    s["wes_studies"] = one("SELECT count(*), sum(n_wes_samples) FROM wes_studies")
    s["wes_suspect_panels"] = one("""SELECT count(*), sum(n_wes_samples) FROM wes_studies
        WHERE n_genes_mutated < 1000 AND n_wes_samples >= 50""")
    s["wes_suspect_panels_top"] = con.execute("""SELECT study_id, n_wes_samples, n_genes_mutated,
        round(muts_per_mutated_sample, 1) FROM wes_studies WHERE n_genes_mutated < 1000 AND n_wes_samples >= 50
        ORDER BY n_wes_samples DESC LIMIT 12""").fetchall()
    rows["wes_labelled_studies.csv"] = copy_csv(con, """SELECT *, n_genes_mutated < 1000 AND n_wes_samples >= 50
        AS looks_like_panel FROM wes_studies ORDER BY looks_like_panel DESC, n_wes_samples DESC""",
        out / "wes_labelled_studies.csv")
    s["tumor_type_source"] = con.execute("""SELECT tumor_type_is_sample_level, count(*) FROM samples_base
        GROUP BY 1 ORDER BY 1""").fetchall()
    s["tumor_type_source_mutated"] = con.execute("""SELECT sb.tumor_type_is_sample_level, count(DISTINCT sb.sample_unique_id)
        FROM samples_base sb SEMI JOIN mutations m ON m.studyId = sb.study_id AND m.sampleId = sb.sample_stable_id
        GROUP BY 1 ORDER BY 1""").fetchall()
    s["study_level_types_top"] = con.execute("""SELECT tumor_type, count(*) FROM samples_base
        WHERE NOT tumor_type_is_sample_level GROUP BY 1 ORDER BY 2 DESC LIMIT 8""").fetchall()
    s["builds_records"] = con.execute("SELECT ncbiBuild, count(*) FROM mutations GROUP BY 1 ORDER BY 2 DESC").fetchall()
    s["builds_studies"] = con.execute("""SELECT builds, count(*) FROM (
        SELECT studyId, string_agg(DISTINCT ncbiBuild, '+' ORDER BY ncbiBuild) AS builds FROM mutations GROUP BY 1)
        GROUP BY 1 ORDER BY 2 DESC""").fetchall()
    s["builds_recurrence_rows_all"] = con.execute(f"""SELECT build, count(*) FROM {vr} WHERE tumor_type = 'ALL'
        GROUP BY 1 ORDER BY 2 DESC""").fetchall()
    s["chrom_values"] = con.execute("SELECT chr, count(*) FROM mutations GROUP BY 1 ORDER BY 2 DESC").fetchall()
    s["recurrence_rows"] = one(f"SELECT count(*), count(*) FILTER (WHERE tumor_type = 'ALL'), "
                               f"count(DISTINCT tumor_type) FROM {vr}")
    s["plausibility"] = con.execute("SELECT * FROM plaus").fetchall()
    s["skcm_braf_by_study"] = con.execute("""
        WITH d AS (SELECT sx.study_id, count(DISTINCT sx.patient_unique_id) AS n_profiled
                   FROM sample_to_gene_panel_derived p JOIN samples_x sx USING (sample_unique_id)
                   WHERE p.alteration_type = 'MUTATION_EXTENDED' AND sx.tumor_type = 'SKCM' GROUP BY 1),
             n AS (SELECT sx.study_id, count(DISTINCT sx.patient_unique_id) AS n_v600e
                   FROM mutations m JOIN samples_x sx ON sx.study_id = m.studyId AND sx.sample_stable_id = m.sampleId
                   WHERE sx.tumor_type = 'SKCM' AND m."GENE.hugoGeneSymbol" = 'BRAF' AND m.proteinChange = 'V600E'
                   GROUP BY 1)
        SELECT d.study_id, d.n_profiled, coalesce(n.n_v600e, 0), round(coalesce(n.n_v600e, 0) / d.n_profiled, 3)
        FROM d LEFT JOIN n USING (study_id) ORDER BY d.n_profiled DESC LIMIT 15""").fetchall()
    s["overlap_studies_msk"] = con.execute("""SELECT study_id, count(DISTINCT person_key) FROM prof
        WHERE tumor_type = 'ALL' AND id_family = 'MSK-DMP' GROUP BY 1 ORDER BY 2 DESC LIMIT 10""").fetchall()
    s["overlap_studies_tcga"] = con.execute("""SELECT study_id, count(DISTINCT person_key) FROM prof
        WHERE tumor_type = 'ALL' AND id_family = 'TCGA' GROUP BY 1 ORDER BY 2 DESC LIMIT 10""").fetchall()
    s["persons_by_n_studies"] = con.execute("""SELECT n_studies, count(*) FROM (
        SELECT person_key, count(DISTINCT study_id) AS n_studies FROM prof
        WHERE tumor_type = 'ALL' AND id_family <> 'study-scoped' GROUP BY 1) GROUP BY 1 ORDER BY 1""").fetchall()
    (out / "summary_stats.json").write_text(json.dumps(s, indent=1, default=str))
    print(json.dumps(rows, indent=1))


if __name__ == "__main__":
    main(*sys.argv[1:4])
