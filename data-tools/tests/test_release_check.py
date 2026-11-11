from panx_data.release_check import Kind, metadata_differences

S3_TABLE = b"accession\tstrain\tcountry\thost\nA\ts1\tunknown\tcow\nB\ts2\tCH\tunknown\n"


def test_metadata_check_takes_unknown_as_an_empty_field() -> None:
    ours = b"accession\tstrain\tcountry\thost\nB\ts2\tCH\t\nA\ts1\t\tcow\n"

    assert list(metadata_differences(S3_TABLE, ours)) == []


def test_metadata_check_lists_left_out_columns_and_fails_on_other_differences() -> None:
    ours = b"accession\tstrain\tcountry\tdate\nA\ts1\tFR\t2013\nC\ts3\t\t\n"

    assert list(metadata_differences(S3_TABLE, ours)) == [
        (Kind.LEFT_OUT_COLUMN, "column 'host' of the S3 metainfo.tsv is not in meta_display_order"),
        (Kind.METADATA_DIFFERENCE, "column 'date' is not in the S3 metainfo.tsv"),
        (Kind.METADATA_DIFFERENCE, "accession 'B' is in one table only"),
        (Kind.METADATA_DIFFERENCE, "accession 'C' is in one table only"),
        (Kind.METADATA_DIFFERENCE, "A 'country': S3 '', release 'FR'"),
    ]
    assert not Kind.LEFT_OUT_COLUMN.fails


def test_metadata_check_reads_quoted_fields() -> None:
    theirs = b'accession\tstrain\tecotype\nA\ts1\t"LLII\tIII"\n'

    assert list(metadata_differences(theirs, theirs)) == []
