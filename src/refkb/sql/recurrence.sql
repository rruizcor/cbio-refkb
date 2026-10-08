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
-- (so they are not double counted).
panel_only_n AS (
    SELECT pr.tumor_type, g.gene, count(DISTINCT pr.patient_unique_id) AS n
    FROM profiled pr
    JOIN gene_panel_to_gene_derived g ON g.gene_panel_id = pr.gene_panel_id
    ANTI JOIN wes_patients w ON w.tumor_type = pr.tumor_type AND w.patient_unique_id = pr.patient_unique_id
    WHERE pr.gene_panel_id <> 'WES'
    GROUP BY pr.tumor_type, g.gene
),
muts AS (
    SELECT
        sx.tumor_type,
        sx.patient_unique_id,
        CASE
            WHEN m.ncbiBuild IN ('GRCh37', '37', 'hg19') THEN 'GRCh37'
            WHEN m.ncbiBuild IN ('GRCh38', '38', 'hg38') THEN 'GRCh38'
            ELSE coalesce(m.ncbiBuild, 'NA')
        END                                   AS build,
        regexp_replace(m.chr, '^chr', '')     AS chrom,
        m.startPosition                       AS start_pos,
        m.endPosition                         AS end_pos,
        m.referenceAllele                     AS ref,
        m.tumorSeqAllele                      AS alt,
        m."GENE.hugoGeneSymbol"               AS gene,
        m.proteinChange                       AS protein_change,
        m.mutationType                        AS mutation_type,
        m.variantType                         AS variant_type,
        lower(coalesce(m.mutationStatus, '')) AS status,
        m.tumorAltCount, m.tumorRefCount, m.normalAltCount, m.normalRefCount
    FROM mutation_derived m
    JOIN samples_x sx ON sx.study_id = m.studyId AND sx.sample_stable_id = m.sampleId
),
agg AS (
    SELECT
        build, chrom, start_pos, end_pos, ref, alt, tumor_type,
        any_value(gene)           AS gene,
        mode(protein_change)      AS protein_change,
        mode(mutation_type)       AS mutation_type,
        mode(variant_type)        AS variant_type,
        count(DISTINCT patient_unique_id) FILTER (WHERE status <> 'germline') AS n_patients,
        count(DISTINCT patient_unique_id) FILTER (WHERE status = 'germline')  AS n_patients_germline_flag,
        count(DISTINCT patient_unique_id) FILTER (
            WHERE status <> 'germline' AND normalAltCount = 0 AND normalRefCount >= {min_normal_depth}
        ) AS n_patients_normal_confirmed,
        count(DISTINCT patient_unique_id) FILTER (
            WHERE status <> 'germline' AND normalAltCount > 0
        ) AS n_patients_normal_alt,
        count(DISTINCT patient_unique_id) FILTER (
            WHERE status <> 'germline' AND normalRefCount IS NOT NULL
        ) AS n_patients_with_normal_counts,
        median(tumorAltCount::DOUBLE / nullif(tumorAltCount + tumorRefCount, 0)) FILTER (
            WHERE status <> 'germline'
        ) AS median_tumor_vaf
    FROM muts
    GROUP BY build, chrom, start_pos, end_pos, ref, alt, tumor_type
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
