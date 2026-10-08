-- Called mutation records with source encodings normalized. Every build query
-- reads this view, not mutation_derived directly.
--   * <study>_mutations_uncalled profiles hold read counts at sites called in
--     other samples, not calls: dropped.
--   * chr: 'chr' prefix stripped; cBioPortal's numeric 23/24 -> X/Y.
--   * ncbiBuild aliases (37, hg19, ...) -> GRCh37 / GRCh38; missing -> 'NA'.
--   * Read counts: a negative value is a deposit error, so the pair is set to NULL.
--     A normal with ref + alt = 0 (or only one of the two counts) carries no
--     matched-normal evidence, so the normal pair is NULL too. "Has normal
--     counts" is then simply normalRefCount IS NOT NULL.
CREATE OR REPLACE TEMP VIEW mutations AS
WITH m AS (
    SELECT * REPLACE (regexp_replace(chr, '^chr', '') AS chr)
    FROM mutation_derived
    WHERE molecularProfileId NOT LIKE '%\_mutations\_uncalled' ESCAPE '\'
)
SELECT * REPLACE (
    CASE chr WHEN '23' THEN 'X' WHEN '24' THEN 'Y' ELSE chr END AS chr,
    CASE
        WHEN ncbiBuild IN ('GRCh37', '37', 'hg19') THEN 'GRCh37'
        WHEN ncbiBuild IN ('GRCh38', '38', 'hg38') THEN 'GRCh38'
        ELSE coalesce(nullif(ncbiBuild, ''), 'NA')
    END AS ncbiBuild,
    CASE WHEN tumorAltCount >= 0 AND tumorRefCount >= 0 THEN tumorAltCount END AS tumorAltCount,
    CASE WHEN tumorAltCount >= 0 AND tumorRefCount >= 0 THEN tumorRefCount END AS tumorRefCount,
    CASE WHEN normalAltCount >= 0 AND normalRefCount >= 0 AND normalAltCount + normalRefCount > 0
         THEN normalAltCount END AS normalAltCount,
    CASE WHEN normalAltCount >= 0 AND normalRefCount >= 0 AND normalAltCount + normalRefCount > 0
         THEN normalRefCount END AS normalRefCount
)
FROM m;
