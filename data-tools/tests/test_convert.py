import gzip
import io
import tarfile

import pytest

from panx_data.archive import GZIP_HEADER, TAR_END, GzipWriter, gzip_member, tar_entry
from panx_data.convert import (
    MetadataConfig,
    cluster_columns,
    cluster_entry,
    date_number,
    fasta_ids,
    gene_tree_json,
    gene_tree_newick,
    metadata_and_strains,
    metadata_tsv,
    newick_tips,
    older_fasta,
    strain_tree,
)

STRAINS = frozenset({"A", "B", "C"})


@pytest.mark.parametrize(
    ("text", "number"),
    [
        ("1922", 1922.5),
        ("2013-01-01", 2013.0),
        ("2012-07-01", 2012 + 182 / 366),
        ("2013-02", 2013 + (31 + 14) / 365),
        ("Feb-2013", 2013 + (31 + 14) / 365),
        ("22-Sep-2015", 2015 + 264 / 365),
        ("2018-05-11T12:00:00Z", 2018 + 130.5 / 365),
        ("1800/2014", (1800.5 + 2014.5) / 2),
    ],
)
def test_date_number_gives_the_decimal_year_or_the_middle_of_the_interval(text: str, number: float) -> None:
    assert date_number(text) == pytest.approx(number, abs=1e-9)


@pytest.mark.parametrize("text", ["unknown", "26", "2.5", "Sept-2013", "2013/unknown", "USA"])
def test_date_number_is_none_for_other_text(text: str) -> None:
    assert date_number(text) is None


@pytest.mark.parametrize("text", ["2013-02-30", "2013-13", "2018-05-11T25:00:00Z", "31-Apr-2015"])
def test_date_number_fails_on_a_date_form_without_calendar_date(text: str) -> None:
    with pytest.raises(ValueError, match="no calendar date|no time of day"):
        date_number(text)


def test_metadata_columns_follow_the_display_order_and_type_rules() -> None:
    rows = [
        {"accession": "A", "strain": "unknown", "date": "2013-01-01", "mic": ">2.00", "size": "1.5", "kind": "x"},
        {"accession": "B", "strain": "s2", "date": "1922", "mic": "unknown", "size": "big", "kind": "y", "geo": "1"},
    ]
    config = MetadataConfig(
        {"mic": [[">2.00", 1.0], ["<0.03", -5.06], ["unknown", "unknown"]]},
        {
            "meta_display_order": ["date", "mic", "size", "kind", "hidden"],
            "color_options": {
                "date": {"type": "continuous"},
                "mic": {"type": "mixed_continuous"},
                "size": {"type": "continuous"},
                "kind": {"type": "discrete"},
                "hidden": {"type": "discrete", "display": "no"},
            },
        },
    )

    metadata, strains = metadata_and_strains(rows, config)

    assert metadata == [
        {"key": "date", "type": "continuous", "scale": "linear", "numbers": {"2013-01-01": 2013.0, "1922": 1922.5}},
        {"key": "mic", "type": "continuous", "scale": "log2", "numbers": {">2.00": 1.0}},
        {"key": "size", "type": "continuous", "scale": "linear", "numbers": {"1.5": 1.5}},
        {"key": "kind", "type": "discrete"},
        {"key": "hidden", "type": "text"},
    ]
    assert strains[0] == {
        "accession": "A",
        "strain": None,
        "values": {"date": "2013-01-01", "mic": ">2.00", "size": "1.5", "kind": "x", "hidden": None},
    }
    assert strains[1]["strain"] == "s2" and strains[1]["values"]["mic"] is None


def test_metadata_tsv_quotes_fields_with_tabs_and_quotes_and_leaves_unknown_empty() -> None:
    metadata = [{"key": "ecotype", "type": "discrete"}, {"key": "note", "type": "text"}]
    strains = [
        {"accession": "A", "strain": None, "values": {"ecotype": "LLII\tIII", "note": 'say "hi"'}},
        {"accession": "B", "strain": "s2", "values": {"ecotype": None, "note": "plain"}},
    ]

    assert metadata_tsv(metadata, strains) == (
        b'accession\tstrain\tecotype\tnote\nA\t\t"LLII\tIII"\t"say ""hi"""\nB\ts2\t\tplain\n'
    )


def test_cluster_entry_maps_the_fields_of_a_cluster_table_row() -> None:
    row = {
        "geneId": 23,
        "geneLen": "486",
        "count": 11,
        "dupli": "yes",
        "dup_detail": "NZ_1#2@NZ_2#3",
        "ann": "dihydrodipicolinate_reductase",
        "msa": "dd40_refined",
        "divers": "0.1",
        "event": "0",
        "allAnn": "dihydrodipicolinate_reductase#6@hypothetical#protein#5@subunit B @ subunit C#2",
        "GName": "None",
        "allGName": "None#9@dapB#5",
        "AZM BA": "0.5",
    }

    assert cluster_entry(row, "protein", ["AZM BA", "CRO BA"]) == {
        "id": "dd40_refined",
        "kind": "protein",
        "strains": 11,
        "length": 486,
        "annotation": "dihydrodipicolinate_reductase",
        "gene": None,
        "annotations": [
            ["dihydrodipicolinate_reductase", 6],
            ["hypothetical#protein", 5],
            ["subunit B @ subunit C", 2],
        ],
        "genes": [[None, 9], ["dapB", 5]],
        "duplications": [["NZ_1", 2], ["NZ_2", 3]],
        "events": 0,
        "diversity": 0.1,
        "extra": {"AZM BA": "0.5", "CRO BA": None},
    }


def test_cluster_entry_fails_when_dupli_and_dup_detail_disagree() -> None:
    row = {"msa": "x", "count": 1, "geneLen": 3, "event": 0, "dupli": "no", "dup_detail": "A#2"}

    with pytest.raises(ValueError, match="disagree"):
        cluster_entry(row, "rna", [])


def test_cluster_columns_skip_the_group_headers_of_the_old_table() -> None:
    config = [
        {"insertion_pos": "Length", "new_col": {"data": "", "group name": "branch association", "tooltip": ""}},
        {"new_col": {"data": "AZM BA", "group": "branch association", "name": "AZM BA", "tooltip": "AZM ba"}},
    ]

    assert cluster_columns(config) == [
        {"key": "AZM BA", "label": "AZM BA", "group": "branch association", "description": "AZM ba"}
    ]


def test_older_fasta_headers_become_accession_and_locus_with_free_text_kept() -> None:
    text = (
        ">consensus <unknown description>\nAC\n"
        ">ECDC_BE13_487-GKMF_01361-9-gph_2-like <unknown description>\nA.\n"
        ">CF02-L_1-1-hyp\nAG\n"
    )

    converted = older_fasta(text)

    assert converted == (
        ">consensus <unknown description>\nAC\n"
        ">ECDC_BE13_487|GKMF_01361 9-gph_2-like <unknown description>\nA.\n"
        ">CF02|L_1 1-hyp\nAG\n"
    )
    assert fasta_ids(converted) == ["consensus", "ECDC_BE13_487|GKMF_01361", "CF02|L_1"]


def test_older_fasta_fails_on_a_header_of_another_form() -> None:
    with pytest.raises(ValueError, match="is not"):
        older_fasta(">A|locus\nAC\n")


def test_strain_tree_keeps_names_lengths_and_children_with_a_zero_root_length() -> None:
    core = {
        "name": "N0",
        "clade": 0,
        "children": [
            {"name": "A", "clade": 1, "branch_length": 0.1, "attr": {"strain": "x"}},
            {"name": "B", "clade": 2, "branch_length": "0.2", "attr": {}},
        ],
    }

    assert strain_tree(core) == {
        "name": "N0",
        "branch_length": 0,
        "children": [{"name": "A", "branch_length": 0.1}, {"name": "B", "branch_length": 0.2}],
    }


def test_gene_tree_json_removes_tips_without_strain_and_joins_the_branches() -> None:
    tree = {
        "name": "N0",
        "branch_length": 0.001,
        "children": [
            {"name": "A|1", "accession": "A", "branch_length": 0.1},
            {
                "name": "N1",
                "branch_length": "0.2",
                "muts": "A1T",
                "children": [
                    {"name": "X|1", "accession": "X", "branch_length": 0.3},
                    {"name": "B|1", "accession": "B", "branch_length": 0.1, "muts": "C2G", "annotation": "hyp"},
                ],
            },
        ],
    }

    pruned, removed = gene_tree_json(tree, STRAINS.__contains__)

    assert removed == 1
    assert pruned == {
        "name": "N0",
        "branch_length": 0.001,
        "children": [
            {"name": "A|1", "accession": "A", "branch_length": 0.1},
            {"name": "B|1", "accession": "B", "branch_length": 0.3, "muts": "A1T,C2G", "annotation": "hyp"},
        ],
    }


def test_gene_tree_newick_removes_tips_without_strain_from_quoted_and_supported_labels() -> None:
    text = "(('A|1-1-3\\'(2\\'),5\\'-x':0.1,'X|2-1-y':0.2)N1:0.99:0.05,B|3-1-z:0.3)N0:0.0001;\n"

    pruned, removed = gene_tree_newick(text, STRAINS.__contains__)

    assert removed == 1
    assert pruned == "('A|1-1-3\\'(2\\'),5\\'-x':0.15,B|3-1-z:0.3)N0:0.0001;\n"
    assert newick_tips(pruned) == ["'A|1-1-3\\'(2\\'),5\\'-x'", "B|3-1-z"]


def test_gene_tree_newick_keeps_a_tree_whose_tips_are_strains() -> None:
    text = "((A|1:0.1,B|1:0.2)NODE_1:0.3,C|1:0.4):0.0;\n"

    assert gene_tree_newick(text, STRAINS.__contains__) == (text, 0)


def test_gzip_member_has_the_fixed_header_and_decompresses_alone() -> None:
    data = b">A|1\nACGT\n" * 100

    member = gzip_member(data)

    assert member[:10] == GZIP_HEADER
    assert gzip.decompress(member) == data


def test_gzip_writer_gives_one_member_of_the_whole_data() -> None:
    stream = io.BytesIO()
    writer = GzipWriter(stream)
    for part in (b"first ", b"second"):
        writer.write(part)
    writer.close()

    assert stream.getvalue() == gzip_member(b"first second")


def test_members_of_tar_entries_read_as_one_archive_and_each_alone() -> None:
    first = gzip_member(tar_entry("D/geneCluster/a_patterns.json", b'{"patterns": "33"}', 1_657_000_000))
    second = gzip_member(tar_entry("D/geneCluster/" + "b" * 80 + "_patterns.json", b"{}", 7))
    pack = first + second + gzip_member(TAR_END)

    with tarfile.open(fileobj=io.BytesIO(pack), mode="r:gz") as tar:
        entries = [(info.name, info.mode, info.uid, info.gid, info.uname, info.gname, info.mtime) for info in tar]
    with tarfile.open(fileobj=io.BytesIO(second), mode="r:gz") as tar:
        alone = [(info.name, tar.extractfile(info).read()) for info in tar]  # type: ignore[union-attr]

    assert entries == [
        ("D/geneCluster/a_patterns.json", 0o644, 0, 0, "", "", 1_657_000_000),
        ("D/geneCluster/" + "b" * 80 + "_patterns.json", 0o644, 0, 0, "", "", 7),
    ]
    assert alone == [("D/geneCluster/" + "b" * 80 + "_patterns.json", b"{}")]
