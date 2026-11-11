from panx_data.release_check import metadata_differences

S3_TABLE = b"accession\tstrain\tcountry\thost\nA\ts1\tunknown\tcow\nB\ts2\tCH\tunknown\n"


def test_metadata_check_takes_unknown_as_an_empty_field() -> None:
    ours = b"accession\tstrain\tcountry\thost\nB\ts2\tCH\t\nA\ts1\t\tcow\n"

    assert list(metadata_differences(S3_TABLE, ours)) == []


def test_metadata_check_reports_other_values_columns_and_accessions() -> None:
    ours = b"accession\tstrain\tcountry\tdate\nA\ts1\tFR\t2013\nC\ts3\t\t\n"

    assert list(metadata_differences(S3_TABLE, ours)) == [
        "column 'host' of the S3 metainfo.tsv is missing",
        "column 'date' is not in the S3 metainfo.tsv",
        "accession 'B' is in one table only",
        "accession 'C' is in one table only",
        "A 'country': S3 '', release 'FR'",
    ]


def test_metadata_check_reads_quoted_fields() -> None:
    theirs = b'accession\tstrain\tecotype\nA\ts1\t"LLII\tIII"\n'

    assert list(metadata_differences(theirs, theirs)) == []
