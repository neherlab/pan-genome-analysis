"""The dataset catalog of a release: the dropdown groups, their labels, and the datasets in display order.

The update zip of the November 2022 run holds the catalog of the old site: `species-list-info.js` with the dropdown
groups and texts, and `all_downloads_table.json` with the species name and source of each 2022 dataset.
`catalog/datasets.toml` adds the category ids and the name and source of the older datasets, which no archive holds.
"""

import json
import re
import tomllib
import zipfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from panx_data.layout import UPDATE_ZIP_PREFIX

CATALOG_CONFIG = Path(__file__).resolve().parent.parent / "catalog" / "datasets.toml"
SPECIES_LIST = "species-list-info.js"
DOWNLOADS_TABLE = "all_downloads_table.json"
NCBI_NUCCORE = "https://www.ncbi.nlm.nih.gov/nuccore/{accession}"

_ASSIGNMENT = re.compile(r"\s*(?:var\s+)?([A-Za-z_$][\w$]*)\s*=\s*")
_SEPARATOR = re.compile(r"\s*([,;]?)\s*")


def js_values(text: str) -> dict[str, Any]:
    """Return the values of JavaScript statements `[var] <name> = <JSON value>`, separated by `,` or `;`.

    Any other syntax fails, so a JavaScript input that is no list of such statements is never misread.
    """
    values: dict[str, Any] = {}
    decoder = json.JSONDecoder()
    position = 0
    while True:
        assignment = _ASSIGNMENT.match(text, position)
        if assignment is None:
            raise ValueError(f"expected '[var] <name> = <JSON value>' at offset {position}: {text[position:][:40]!r}")
        name = assignment.group(1)
        if name in values:
            raise ValueError(f"'{name}' is assigned twice")
        try:
            values[name], position = decoder.raw_decode(text, assignment.end())
        except json.JSONDecodeError as error:
            raise ValueError(f"the value of '{name}' is no JSON value: {error}") from error
        separator = _SEPARATOR.match(text, position)
        assert separator is not None
        position = separator.end()
        if position == len(text):
            return values
        if not separator.group(1):
            raise ValueError(f"expected ',' or ';' at offset {position}: {text[position:][:40]!r}")


@dataclass(frozen=True)
class Category:
    """One dropdown group."""

    id: str
    label: str


@dataclass(frozen=True)
class CatalogDataset:
    """One dataset of the catalog, with the fields of its `index.json` entry and the link of its accessions."""

    id: str
    name: str
    label: str
    category: str
    source: str
    accession_url: str | None


@dataclass(frozen=True)
class Catalog:
    """The categories and the datasets, each in display order."""

    categories: tuple[Category, ...]
    datasets: tuple[CatalogDataset, ...]

    def dataset(self, dataset_id: str) -> CatalogDataset:
        for dataset in self.datasets:
            if dataset.id == dataset_id:
                return dataset
        raise KeyError(f"'{dataset_id}' is no dataset of the catalog")


def load_catalog(update_zip: Path, config: Path = CATALOG_CONFIG) -> Catalog:
    """Read the catalog from the update zip and the catalog configuration."""
    with zipfile.ZipFile(update_zip) as update:
        species_list = update.read(f"{UPDATE_ZIP_PREFIX}{SPECIES_LIST}").decode()
        downloads = json.loads(update.read(f"{UPDATE_ZIP_PREFIX}{DOWNLOADS_TABLE}"))
    return read_catalog(species_list, downloads, tomllib.loads(config.read_text()))


def read_catalog(species_list: str, downloads: dict[str, Any], config: dict[str, Any]) -> Catalog:
    """Join the dropdown groups of `species_list`, the rows of the downloads table, and the configuration.

    Every group member needs a name and a source, from its downloads row or from the configuration, never both, and
    every downloads row must belong to a group member.
    """
    species = js_values(species_list)["species_dt"]
    groups = {key: value for key, value in species.items() if isinstance(value, dict) and "members" in value}
    category_ids: dict[str, str] = config["categories"]
    if set(groups) != set(category_ids):
        raise ValueError(
            f"dropdown groups {sorted(groups)} differ from the configured categories {sorted(category_ids)}"
        )
    rows = {row["panX link"]: row for row in downloads["data"]}
    if len(rows) != len(downloads["data"]):
        raise ValueError(f"'{DOWNLOADS_TABLE}' has two rows for one dataset")
    configured: dict[str, dict[str, str]] = config.get("datasets", {})
    categories = tuple(Category(category_ids[key], group["label"].removesuffix(":")) for key, group in groups.items())
    datasets: list[CatalogDataset] = []
    errors: list[str] = []
    for key, group in groups.items():
        for member in group["members"]:
            row, extra = rows.get(member), configured.get(member)
            label = species.get(member)
            if not isinstance(label, str):
                errors.append(f"'{member}' has no dropdown text")
            elif row is not None and extra is not None:
                errors.append(f"'{member}' has a downloads row and a configured name")
            elif row is not None:
                name, source, link = row["species name"], row["source"], NCBI_NUCCORE
                datasets.append(CatalogDataset(member, name, label, category_ids[key], source, link))
            elif extra is not None:
                name, source, link = extra["name"], extra["source"], extra.get("accession_url")
                datasets.append(CatalogDataset(member, name, label, category_ids[key], source, link))
            else:
                errors.append(f"'{member}' has neither a downloads row nor a configured name and source")
    members = [dataset.id for dataset in datasets]
    if len(set(members)) != len(members):
        errors.append("a dataset is a member of two dropdown groups")
    errors += [f"downloads row '{name}' belongs to no dropdown group" for name in sorted(set(rows) - set(members))]
    errors += [
        f"configured dataset '{name}' belongs to no dropdown group" for name in sorted(set(configured) - set(members))
    ]
    if errors:
        raise ValueError("the catalog is inconsistent:\n" + "\n".join(errors))
    return Catalog(categories, tuple(datasets))
