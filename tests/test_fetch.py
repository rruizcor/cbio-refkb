from refkb.fetch import table_files

FILES = [
    "mutation_derived/latest-00000-of-00002.parquet",
    "mutation_derived/latest-00001-of-00002.parquet",
    "mutation/latest-00000-of-00001.parquet",
    "data/sample_derived/part-0.parquet",
    "manifest.json",
]


def test_table_files_matches_exact_table_only():
    assert table_files(FILES, "mutation") == ["mutation/latest-00000-of-00001.parquet"]
    assert len(table_files(FILES, "mutation_derived")) == 2
    assert table_files(FILES, "sample_derived") == ["data/sample_derived/part-0.parquet"]
    assert table_files(FILES, "gene") == []
