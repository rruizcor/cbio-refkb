-- Sanity check for the denominator logic: which gene panels samples were
-- profiled with for mutations, and how many mutated samples have no panel row.
-- The recurrence build assumes gene_panel_id = 'WES' means "all genes covered".
SELECT 'panel' AS check_type, gene_panel_id AS key, count(DISTINCT sample_unique_id) AS n_samples
FROM sample_to_gene_panel_derived
WHERE alteration_type = 'MUTATION_EXTENDED'
GROUP BY gene_panel_id

UNION ALL

SELECT 'mutated_sample_without_panel_row', NULL, count(DISTINCT sb.sample_unique_id)
FROM mutation_derived m
JOIN samples_base sb ON sb.study_id = m.studyId AND sb.sample_stable_id = m.sampleId
WHERE sb.sample_unique_id NOT IN (
    SELECT sample_unique_id FROM sample_to_gene_panel_derived WHERE alteration_type = 'MUTATION_EXTENDED'
)

UNION ALL

SELECT 'tumor_type_source', CASE WHEN tumor_type_is_sample_level THEN 'sample_oncotree' ELSE 'study_level' END,
       count(*)
FROM samples_base
GROUP BY 2

ORDER BY check_type, n_samples DESC;
