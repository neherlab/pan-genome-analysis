import pytest

from panx_data.catalog import NCBI_NUCCORE, Catalog, CatalogDataset, Category, js_values, read_catalog

SPECIES_LIST = """var species_dt = {
 "optgroup1": {"label": "Case studies:", "members": ["Study_one"]},
 "optgroup2": {"label": "Pan-genomes from RefSeq records:", "members": ["Testus_alpha", "Testus_beta"]},
 "Study_one": "A study (Someone et al)",
 "Testus_alpha": "Testus alpha (3)",
 "Testus_beta": "Testus beta (2)",
 "Alias": "A study (Someone et al)"
};
var species_search_dt = [{"species": "Testus alpha", "linkName": "Testus_alpha"}]"""
CONFIG = {
    "categories": {"optgroup1": "case_study", "optgroup2": "refseq"},
    "datasets": {"Study_one": {"name": "A study (Someone et al)", "source": "Someone et al"}},
}


def _downloads(*names: str) -> dict:
    return {"data": [{"panX link": name, "species name": name.replace("_", " "), "source": "RefSeq"} for name in names]}


def test_js_values_reads_statements_separated_by_comma_or_semicolon() -> None:
    text = 'var meta_details={"a": [1]},\nmeta_display={"b": "x;y"};\nvar association_columns=[];\n'

    assert js_values(text) == {"meta_details": {"a": [1]}, "meta_display": {"b": "x;y"}, "association_columns": []}


@pytest.mark.parametrize(
    "text",
    [
        "var a = {'single': 1};",
        "var a = {b: 1};",
        "// comment\nvar a = 1;",
        "var a = 1 var b = 2",
        "var a = 1; var a = 2;",
        "",
    ],
)
def test_js_values_rejects_other_javascript(text: str) -> None:
    with pytest.raises(ValueError):
        js_values(text)


def test_catalog_joins_groups_rows_and_configuration_in_display_order() -> None:
    catalog = read_catalog(SPECIES_LIST, _downloads("Testus_alpha", "Testus_beta"), CONFIG)

    assert catalog == Catalog(
        (Category("case_study", "Case studies"), Category("refseq", "Pan-genomes from RefSeq records")),
        (
            CatalogDataset(
                "Study_one", "A study (Someone et al)", "A study (Someone et al)", "case_study", "Someone et al", None
            ),
            CatalogDataset("Testus_alpha", "Testus alpha", "Testus alpha (3)", "refseq", "RefSeq", NCBI_NUCCORE),
            CatalogDataset("Testus_beta", "Testus beta", "Testus beta (2)", "refseq", "RefSeq", NCBI_NUCCORE),
        ),
    )


@pytest.mark.parametrize(
    ("rows", "message"),
    [
        (("Testus_alpha",), "'Testus_beta' has neither a downloads row nor a configured name and source"),
        (("Testus_alpha", "Testus_beta", "Testus_gamma"), "downloads row 'Testus_gamma' belongs to no dropdown group"),
        (("Testus_alpha", "Testus_beta", "Study_one"), "'Study_one' has a downloads row and a configured name"),
    ],
)
def test_catalog_fails_on_members_without_name_or_rows_without_member(rows: tuple[str, ...], message: str) -> None:
    with pytest.raises(ValueError, match=message):
        read_catalog(SPECIES_LIST, _downloads(*rows), CONFIG)


def test_catalog_fails_on_a_group_without_category_id() -> None:
    config = {**CONFIG, "categories": {"optgroup1": "case_study"}}

    with pytest.raises(ValueError, match="configured categories"):
        read_catalog(SPECIES_LIST, _downloads("Testus_alpha", "Testus_beta"), config)
