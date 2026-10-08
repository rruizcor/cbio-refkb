-- How much germline-flagged and matched-normal information exists, per study.
-- Run this first on a new snapshot: it tells you which studies can feed the
-- "normal-confirmed somatic" set and how sparse the germline flags really are.
SELECT
    m.studyId                                                        AS study_id,
    any_value(cs.type_of_cancer_id)                                  AS study_cancer_type,
    count(*)                                                         AS n_mutation_records,
    count(DISTINCT m.sampleId)                                       AS n_samples_with_mutations,
    count(*) FILTER (WHERE lower(m.mutationStatus) = 'somatic')      AS n_somatic,
    count(*) FILTER (WHERE lower(m.mutationStatus) = 'germline')     AS n_germline,
    count(*) FILTER (WHERE lower(m.mutationStatus) = 'loh')          AS n_loh,
    count(*) FILTER (WHERE m.mutationStatus IS NULL
                       OR lower(m.mutationStatus) NOT IN ('somatic','germline','loh')) AS n_status_other_or_null,
    round(100.0 * count(*) FILTER (WHERE m.normalRefCount IS NOT NULL) / count(*), 1)
                                                                     AS pct_with_normal_counts,
    count(*) FILTER (WHERE lower(coalesce(m.mutationStatus,'')) <> 'germline'
                       AND m.normalAltCount = 0
                       AND m.normalRefCount >= {min_normal_depth})   AS n_normal_confirmed_somatic,
    count(*) FILTER (WHERE lower(coalesce(m.mutationStatus,'')) <> 'germline'
                       AND m.normalAltCount > 0)                     AS n_somatic_with_normal_alt_reads,
    string_agg(DISTINCT coalesce(m.ncbiBuild, 'NA'), ',')            AS ncbi_builds
FROM mutations m
LEFT JOIN cancer_study cs ON cs.cancer_study_identifier = m.studyId
GROUP BY m.studyId
ORDER BY n_germline DESC, pct_with_normal_counts DESC, n_mutation_records DESC;
