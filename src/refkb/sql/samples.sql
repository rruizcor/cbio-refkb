-- One row per sample with its tumor type, plus a duplicate 'ALL' row for pan-cancer stats.
-- Tumor type = sample-level ONCOTREE_CODE when present, else the study's type_of_cancer_id.
-- Upper-cased so sample-level and study-level codes share one namespace.
CREATE OR REPLACE TEMP VIEW samples_base AS
WITH oc AS (
    SELECT sample_unique_id, any_value(attribute_value) AS oncotree_code
    FROM clinical_data_derived
    WHERE attribute_name = 'ONCOTREE_CODE' AND type = 'sample'
    GROUP BY sample_unique_id
)
SELECT
    s.sample_unique_id,
    s.sample_stable_id,
    s.patient_unique_id,
    s.cancer_study_identifier AS study_id,
    cs.type_of_cancer_id      AS study_cancer_type,
    upper(coalesce(nullif(trim(oc.oncotree_code), ''), cs.type_of_cancer_id)) AS tumor_type,
    (oc.oncotree_code IS NOT NULL AND trim(oc.oncotree_code) <> '') AS tumor_type_is_sample_level
FROM sample_derived s
JOIN cancer_study cs ON cs.cancer_study_identifier = s.cancer_study_identifier
LEFT JOIN oc USING (sample_unique_id);

CREATE OR REPLACE TEMP VIEW samples_x AS
SELECT * FROM samples_base
UNION ALL
SELECT * REPLACE ('ALL' AS tumor_type) FROM samples_base;
