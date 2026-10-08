-- Patient-level, panel-aware SNV/indel recurrence by tumor type (+ 'ALL').
--
-- Numerator  : distinct patients with the variant (germline-flagged calls excluded).
-- Denominator: distinct patients in that tumor type whose sequencing covered the gene
--              (WES samples cover every gene; panel samples only their panel's genes).
-- Extra columns support tumor-only germline/somatic triage:
--   n_patients_normal_confirmed : somatic call, 0 alt reads in matched normal at >= {min_normal_depth}x
--   n_patients_normal_alt       : called somatic but alt reads present in the normal
--                                 (germline leak, CHIP, tumor-in-normal)
--   n_patients_germline_flag    : deposited with mutation_status = Germline
--   median_tumor_vaf            : median over somatic patients of each patient's mean VAF
CREATE OR REPLACE TABLE variant_recurrence AS
WITH
profiled AS (
    SELECT DISTINCT sx.tumor_type, sx.patient_unique_id, p.gene_panel_id
    FROM sample_to_gene_panel_derived p
    JOIN samples_x sx USING (sample_unique_id)
    WHERE p.alteration_type = 'MUTATION_EXTENDED'
),
wes_patients AS (
    SELECT DISTINCT tumor_type, patient_unique_id FROM profiled WHERE gene_panel_id = 'WES'
),
wes_n AS (
    SELECT tumor_type, count(*) AS n FROM wes_patients GROUP BY tumor_type
),
-- Panel-profiled patients that have no WES sample in the same tumor type
-- (so they are not double counted). Patients are grouped by the set of panels
-- they were sequenced with, and only those few distinct sets are expanded to
-- genes: expanding every patient to every panel gene is ~225M rows.
panel_sets AS (
    SELECT pr.tumor_type, pr.patient_unique_id, list_sort(list(DISTINCT pr.gene_panel_id)) AS panels
    FROM profiled pr
    ANTI JOIN wes_patients w ON w.tumor_type = pr.tumor_type AND w.patient_unique_id = pr.patient_unique_id
    WHERE pr.gene_panel_id <> 'WES'
    GROUP BY pr.tumor_type, pr.patient_unique_id
),
panel_set_n AS (
    SELECT tumor_type, panels, count(*) AS n FROM panel_sets GROUP BY tumor_type, panels
),
panel_set_genes AS (
    SELECT DISTINCT s.panels, g.gene
    FROM (SELECT DISTINCT panels, unnest(panels) AS gene_panel_id FROM panel_set_n) s
    JOIN gene_panel_to_gene_derived g ON g.gene_panel_id = s.gene_panel_id
),
panel_only_n AS (
    SELECT sn.tumor_type, sg.gene, sum(sn.n) AS n
    FROM panel_set_n sn
    JOIN panel_set_genes sg ON sg.panels = sn.panels
    GROUP BY sn.tumor_type, sg.gene
),
muts AS (
    SELECT
        sx.tumor_type,
        sx.patient_unique_id,
        m.ncbiBuild                           AS build,
        m.chr                                 AS chrom,
        m.startPosition                       AS start_pos,
        m.endPosition                         AS end_pos,
        m.referenceAllele                     AS ref,
        m.tumorSeqAllele                      AS alt,
        m."GENE.hugoGeneSymbol"               AS gene,
        m.proteinChange                       AS protein_change,
        m.mutationType                        AS mutation_type,
        m.variantType                         AS variant_type,
        lower(coalesce(m.mutationStatus, '')) = 'germline' AS germline,
        m.tumorAltCount::DOUBLE / nullif(m.tumorAltCount + m.tumorRefCount, 0) AS vaf,
        m.normalAltCount, m.normalRefCount
    FROM mutations m
    JOIN samples_x sx ON sx.study_id = m.studyId AND sx.sample_stable_id = m.sampleId
),
-- One row per variant x tumor type x patient (collapses multi-sample patients and
-- repeated records), so the counts below are plain filtered counts. A
-- COUNT(DISTINCT patient) per filter needs one distinct hash table each and does
-- not fit a 16 GB machine on the full snapshot.
pv AS MATERIALIZED (
    SELECT
        build, chrom, start_pos, end_pos, ref, alt, tumor_type, patient_unique_id,
        any_value(gene)           AS gene,
        any_value(protein_change) AS protein_change,
        any_value(mutation_type)  AS mutation_type,
        any_value(variant_type)   AS variant_type,
        bool_or(NOT germline)     AS somatic,
        bool_or(germline)         AS germline,
        bool_or(NOT germline AND normalAltCount = 0 AND normalRefCount >= {min_normal_depth}) AS normal_confirmed,
        bool_or(NOT germline AND normalAltCount > 0)                                          AS normal_alt,
        bool_or(NOT germline AND normalRefCount IS NOT NULL)                                  AS with_normal,
        avg(vaf) FILTER (WHERE NOT germline) AS vaf
    FROM muts
    GROUP BY build, chrom, start_pos, end_pos, ref, alt, tumor_type, patient_unique_id
),
-- mode(), median() and string min()/max() keep per-group heap state that DuckDB
-- neither caps nor spills; over all ~20M groups that exhausts a 16 GB machine.
-- So this pass uses fixed-size states only: when every patient row carries the
-- same annotation (equal min/max hash) any_value() is the mode, and with <= 2
-- VAFs (min + max) / 2 is the median. mode()/median() then run only on the
-- groups that need them (~1% and ~10% of groups).
counts AS MATERIALIZED (
    SELECT
        build, chrom, start_pos, end_pos, ref, alt, tumor_type,
        any_value(gene)                          AS gene,
        any_value(protein_change)                AS protein_change,
        any_value(mutation_type)                 AS mutation_type,
        any_value(variant_type)                  AS variant_type,
        count(*) FILTER (WHERE somatic)          AS n_patients,
        count(*) FILTER (WHERE germline)         AS n_patients_germline_flag,
        count(*) FILTER (WHERE normal_confirmed) AS n_patients_normal_confirmed,
        count(*) FILTER (WHERE normal_alt)       AS n_patients_normal_alt,
        count(*) FILTER (WHERE with_normal)      AS n_patients_with_normal_counts,
        (min(vaf) + max(vaf)) / 2                AS median_tumor_vaf,
        min(hash(protein_change)) <> max(hash(protein_change))
          OR min(hash(mutation_type)) <> max(hash(mutation_type))
          OR min(hash(variant_type)) <> max(hash(variant_type)) AS needs_mode,
        count(vaf) >= 3                          AS needs_median
    FROM pv
    GROUP BY build, chrom, start_pos, end_pos, ref, alt, tumor_type
),
holistic AS (
    SELECT
        pv.build, pv.chrom, pv.start_pos, pv.end_pos, pv.ref, pv.alt, pv.tumor_type,
        mode(pv.protein_change) FILTER (WHERE c.needs_mode) AS protein_change,
        mode(pv.mutation_type)  FILTER (WHERE c.needs_mode) AS mutation_type,
        mode(pv.variant_type)   FILTER (WHERE c.needs_mode) AS variant_type,
        median(pv.vaf)          FILTER (WHERE c.needs_median) AS median_tumor_vaf
    FROM pv
    JOIN (SELECT * FROM counts WHERE needs_mode OR needs_median) c
        ON c.build IS NOT DISTINCT FROM pv.build
        AND c.chrom IS NOT DISTINCT FROM pv.chrom
        AND c.start_pos IS NOT DISTINCT FROM pv.start_pos
        AND c.end_pos IS NOT DISTINCT FROM pv.end_pos
        AND c.ref IS NOT DISTINCT FROM pv.ref
        AND c.alt IS NOT DISTINCT FROM pv.alt
        AND c.tumor_type IS NOT DISTINCT FROM pv.tumor_type
    GROUP BY pv.build, pv.chrom, pv.start_pos, pv.end_pos, pv.ref, pv.alt, pv.tumor_type
),
agg AS (
    SELECT c.* EXCLUDE (needs_mode, needs_median) REPLACE (
        coalesce(h.protein_change, c.protein_change)     AS protein_change,
        coalesce(h.mutation_type, c.mutation_type)       AS mutation_type,
        coalesce(h.variant_type, c.variant_type)         AS variant_type,
        coalesce(h.median_tumor_vaf, c.median_tumor_vaf) AS median_tumor_vaf
    )
    FROM counts c
    LEFT JOIN holistic h
        ON h.build IS NOT DISTINCT FROM c.build
        AND h.chrom IS NOT DISTINCT FROM c.chrom
        AND h.start_pos IS NOT DISTINCT FROM c.start_pos
        AND h.end_pos IS NOT DISTINCT FROM c.end_pos
        AND h.ref IS NOT DISTINCT FROM c.ref
        AND h.alt IS NOT DISTINCT FROM c.alt
        AND h.tumor_type IS NOT DISTINCT FROM c.tumor_type
)
SELECT
    a.*,
    coalesce(w.n, 0) + coalesce(p.n, 0) AS n_patients_profiled,
    CASE WHEN coalesce(w.n, 0) + coalesce(p.n, 0) > 0
         THEN a.n_patients::DOUBLE / (coalesce(w.n, 0) + coalesce(p.n, 0)) END AS patient_frequency,
    CASE WHEN a.n_patients_with_normal_counts > 0
         THEN a.n_patients_normal_confirmed::DOUBLE / a.n_patients_with_normal_counts END AS frac_normal_confirmed
FROM agg a
LEFT JOIN wes_n w        ON w.tumor_type = a.tumor_type
LEFT JOIN panel_only_n p ON p.tumor_type = a.tumor_type AND p.gene = a.gene;
