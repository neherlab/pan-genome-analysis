"""Pure conversions of the pipeline output of one dataset to the files of the data format (`docs/data-format.md` of
the app repository `neherlab/pan-genome-visualization`).

The inputs are the parsed files of the old layout: `coreGenomeTree.json`, `strainMetainfo.json`, the metadata
configuration of `metaConfiguration.js`, and the rows of `geneCluster.json`, and the texts of the cluster files.
"""

import calendar
import json
import math
import re
from collections.abc import Callable, Iterable, Mapping
from dataclasses import dataclass
from datetime import date
from decimal import Decimal, InvalidOperation
from typing import Any

UNKNOWN = "unknown"
NONE = "None"
CONSENSUS = "consensus"
MONTHS = {name: number for number, name in enumerate(calendar.month_abbr) if name}
_YEAR_MONTH_DAY = re.compile(r"(\d{4})(?:-(\d{2})(?:-(\d{2})(?:T(\d{2}):(\d{2}):(\d{2})Z)?)?)?")
_DAY_MONTH_YEAR = re.compile(r"(?:(\d{1,2})-)?([A-Z][a-z]{2})-(\d{4})")
_NUMBER = re.compile(r"[+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?")
_INTEGER = re.compile(r"[+-]?\d+")
_PAIR = re.compile(r"(.*?)#(\d+)(?:@|$)", re.DOTALL)


@dataclass(frozen=True)
class MetadataConfig:
    """The metadata configuration of a dataset: the values of each key, and the display options."""

    details: Mapping[str, Any]
    display: Mapping[str, Any]


def strain_tree(node: Mapping[str, Any]) -> dict[str, Any]:
    """Keep `name`, `branch_length` (0 on the root), and `children` of the nodes of `coreGenomeTree.json`."""

    def convert(node: Mapping[str, Any], length: float) -> dict[str, Any]:
        result: dict[str, Any] = {"name": node["name"], "branch_length": length}
        if node.get("children"):
            result["children"] = [convert(child, _branch_length(child)) for child in node["children"]]
        return result

    return convert(node, 0)


def tree_tips(tree: Mapping[str, Any]) -> list[str]:
    """Return the tip names of a tree of nested `children`, in pre-order."""
    tips: list[str] = []
    stack = [tree]
    while stack:
        node = stack.pop()
        children = node.get("children")
        if children:
            stack.extend(reversed(children))
        else:
            tips.append(node["name"])
    return tips


def non_root_nodes(tree: Mapping[str, Any]) -> int:
    """Count the nodes of a tree without its root."""
    count, stack = 0, list(tree.get("children") or [])
    while stack:
        node = stack.pop()
        count += 1
        stack.extend(node.get("children") or [])
    return count


def metadata_and_strains(
    rows: list[Mapping[str, Any]], config: MetadataConfig
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Return the `metadata` columns and the `strains` of `dataset.json`.

    The columns are the keys of `meta_display_order`, in its order. A value `unknown` becomes `null`, also as a
    strain name.
    """
    keys = list(config.display["meta_display_order"])
    if len(set(keys)) != len(keys):
        raise ValueError("meta_display_order names a key twice")
    strains = [
        {
            "accession": row["accession"],
            "strain": _value(row.get("strain")),
            "values": {key: _value(row.get(key)) for key in keys},
        }
        for row in rows
    ]
    options = config.display.get("color_options", {})
    columns = []
    for key in keys:
        values = list(dict.fromkeys(v for strain in strains if (v := strain["values"][key]) is not None))
        columns.append(_column(key, options.get(key, {}), config.details.get(key, []), values))
    return columns, strains


def date_number(text: str) -> float | None:
    """Return the decimal year of a date, a year, a month, or a range `<a>/<b>` of two of them; None for other text.

    A year or a month gives the middle of its interval. Text of date form that is no calendar date fails.
    """
    first, separator, second = text.partition("/")
    if separator:
        a, b = _date_interval(first), _date_interval(second)
        return None if a is None or b is None else (sum(a) / 2 + sum(b) / 2) / 2
    interval = _date_interval(text)
    return None if interval is None else sum(interval) / 2


def plain_number(text: str) -> float | None:
    """Return the number that a text writes, or None."""
    if not _NUMBER.fullmatch(text):
        return None
    number = float(text)
    return number if math.isfinite(number) else None


def cluster_entry(row: Mapping[str, Any], kind: str, extra_keys: Iterable[str]) -> dict[str, Any]:
    """Return the `clusters.json` entry of a row of `geneCluster.json`, without its `range`."""
    duplications = _pairs(row.get("dup_detail"))
    if bool(duplications) != (row.get("dupli") == "yes"):
        raise ValueError(f"cluster '{row['msa']}': dupli '{row.get('dupli')}' and dup_detail disagree")
    return {
        "id": row["msa"],
        "kind": kind,
        "strains": _integer(row["count"], "count"),
        "length": _integer(row["geneLen"], "geneLen"),
        "annotation": _name(row.get("ann")),
        "gene": _name(row.get("GName")),
        "annotations": _pairs(row.get("allAnn")),
        "genes": _pairs(row.get("allGName")),
        "duplications": duplications,
        "events": _integer(row["event"], "event"),
        "diversity": _optional_number(row.get("divers")),
        "extra": {key: _extra_value(row.get(key), key) for key in extra_keys},
    }


def cluster_columns(new_columns: Iterable[Mapping[str, Any]]) -> list[dict[str, str]]:
    """Return the `cluster_columns` of `dataset.json` from the entries of a `newColumnConfig.js`.

    Entries without `data` are the group headers of the old table, which the group of each column replaces.
    """
    columns = [
        {"key": column["data"], "label": column["name"], "group": column["group"], "description": column["tooltip"]}
        for entry in new_columns
        if (column := entry["new_col"])["data"]
    ]
    keys = [column["key"] for column in columns]
    if len(set(keys)) != len(keys):
        raise ValueError("two cluster columns have one key")
    return columns


def fasta_ids(text: str) -> list[str]:
    """Return the record identifiers of a FASTA text: each header without `>` and without the free text after it."""
    return [line[1:].split(" ", 1)[0] for line in text.split("\n") if line.startswith(">")]


def older_fasta(text: str) -> str:
    """Write the FASTA headers `<accession>-<locus>-<number>-<annotation>` of an older dataset as
    `<accession>|<locus> <number>-<annotation>`, and keep any free text after them."""

    def header(line: str) -> str:
        if not line.startswith(">"):
            return line
        identifier, space, free = line[1:].partition(" ")
        if identifier == CONSENSUS:
            return line
        parts = identifier.split("-", 3)
        if len(parts) != 4:
            raise ValueError(f"FASTA header '{line}' is not '<accession>-<locus>-<number>-<annotation>'")
        accession, locus, number, annotation = parts
        return f">{accession}|{locus} {number}-{annotation}{space}{free}"

    return "\n".join(header(line) for line in text.split("\n"))


def gene_tree_json(tree: dict[str, Any], keep: Callable[[str], bool]) -> tuple[dict[str, Any] | None, int]:
    """Return a `_tree.json` tree with numbers as numbers and without the tips whose accession `keep` rejects, and
    the number of removed tips.

    A node left with one child is joined with it: the branch lengths add up, and the mutations of both branches
    are kept. A tree without tips gives None.
    """
    removed = 0

    def convert(node: dict[str, Any]) -> dict[str, Any] | None:
        nonlocal removed
        result = dict(node)
        if "branch_length" in result:
            result["branch_length"] = _number(result["branch_length"], "branch_length")
        children = node.get("children")
        if not children:
            if not keep(node["accession"]):
                removed += 1
                return None
            return result
        kept = [child for child in map(convert, children) if child is not None]
        if not kept:
            return None
        if len(kept) == 1 and len(children) > 1:
            return _join_json(result, kept[0])
        result["children"] = kept
        return result

    return convert(tree), removed


def gene_tree_newick(text: str, keep: Callable[[str], bool]) -> tuple[str | None, int]:
    """Return a Newick gene tree without the tips whose accession `keep` rejects, and the number of removed tips.

    A tip's accession is its label up to the first `|` or `:`. A node left with one child is joined with it, and the
    branch lengths add up. A tree without tips gives None.
    """
    root, tail = _parse_newick(text)
    removed = 0

    def prune(node: _NewickNode) -> _NewickNode | None:
        nonlocal removed
        if node.children is None:
            if keep(_tip_accession(node.label)):
                return node
            removed += 1
            return None
        kept = [child for child in map(prune, node.children) if child is not None]
        if not kept:
            return None
        if len(kept) == 1 and len(node.children) > 1:
            return _join_newick(node, kept[0])
        return _NewickNode(node.label, kept)

    pruned = prune(root)
    return (None if pruned is None else _write_newick(pruned) + tail), removed


def metadata_tsv(metadata: list[Mapping[str, Any]], strains: list[Mapping[str, Any]]) -> bytes:
    """Render `metadata.tsv`: a header row, then one row per strain; unknown values are empty fields."""
    keys = [column["key"] for column in metadata]
    rows = [["accession", "strain", *keys]] + [
        [strain["accession"], strain["strain"], *(strain["values"][key] for key in keys)] for strain in strains
    ]
    return "".join("\t".join(map(tsv_field, row)) + "\n" for row in rows).encode()


def tsv_field(value: str | None) -> str:
    """Write a TSV field: empty for None, in double quotes when it holds a tab, a line break, or `"`."""
    if value is None:
        return ""
    if any(character in value for character in '\t\n\r"'):
        return '"' + value.replace('"', '""') + '"'
    return value


def _column(key: str, option: Mapping[str, Any], details: list[Any], values: list[str]) -> dict[str, Any]:
    kind = option.get("type")
    if kind is None or option.get("display") == "no":
        return {"key": key, "type": "text"}
    if kind == "discrete":
        return {"key": key, "type": "discrete"}
    if kind == "continuous":
        return {"key": key, "type": "continuous", "scale": "linear", "numbers": _linear_numbers(values)}
    if kind == "mixed_continuous":
        return {"key": key, "type": "continuous", "scale": "log2", "numbers": _log2_numbers(key, details, values)}
    raise ValueError(f"metadata column '{key}' has the unknown type '{kind}'")


def _linear_numbers(values: list[str]) -> dict[str, float]:
    """Give a date column the decimal year of each value, every other column the number that a value writes."""
    dates = {value: date_number(value) for value in values}
    if values and all(number is not None for number in dates.values()):
        return {value: number for value, number in dates.items() if number is not None}
    return {value: number for value in values if (number := plain_number(value)) is not None}


def _log2_numbers(key: str, details: list[Any], values: list[str]) -> dict[str, float]:
    """Take the pairs `[text, log2 number]` of a column's values, without the pairs that hold no finite number."""
    numbers: dict[str, float] = {}
    present = set(values)
    for pair in details:
        text, number = pair
        if text not in present or isinstance(number, bool) or not isinstance(number, int | float):
            continue
        if not math.isfinite(number):
            continue
        if numbers.get(text, number) != number:
            raise ValueError(f"metadata column '{key}' gives '{text}' two numbers")
        numbers[text] = number
    return numbers


def _date_interval(text: str) -> tuple[float, float] | None:
    """Return the decimal years of the start and end of a date's interval (equal for a full date), or None."""
    if match := _YEAR_MONTH_DAY.fullmatch(text):
        year, month, day, hour, minute, second = (None if part is None else int(part) for part in match.groups())
    elif (match := _DAY_MONTH_YEAR.fullmatch(text)) and match.group(2) in MONTHS:
        day = None if match.group(1) is None else int(match.group(1))
        year, month, hour, minute, second = int(match.group(3)), MONTHS[match.group(2)], None, None, None
    else:
        return None
    assert year is not None
    if month is None:
        _calendar_date(text, year, 1, 1)
        return float(year), float(year + 1)
    days = 366 if calendar.isleap(year) else 365
    if day is None:
        start = _calendar_date(text, year, month, 1).timetuple().tm_yday - 1
        length = calendar.monthrange(year, month)[1]
        return year + start / days, year + (start + length) / days
    seconds = 0 if hour is None else _time_of_day(text, hour, minute or 0, second or 0)
    point = year + (_calendar_date(text, year, month, day).timetuple().tm_yday - 1 + seconds / 86400) / days
    return point, point


def _calendar_date(text: str, year: int, month: int, day: int) -> date:
    try:
        return date(year, month, day)
    except ValueError as error:
        raise ValueError(f"'{text}' has the form of a date but is no calendar date: {error}") from error


def _time_of_day(text: str, hour: int, minute: int, second: int) -> int:
    if hour > 23 or minute > 59 or second > 59:
        raise ValueError(f"'{text}' has the form of a date-time but is no time of day")
    return hour * 3600 + minute * 60 + second


def _value(value: Any) -> str | None:
    if value is None or value == UNKNOWN:
        return None
    if isinstance(value, str):
        return value
    if isinstance(value, int | float) and not isinstance(value, bool):
        return json.dumps(value)
    raise ValueError(f"metadata value {value!r} is no string or number")


def _name(value: Any) -> str | None:
    return None if value is None or value == NONE else str(value)


def _pairs(text: Any) -> list[list[Any]]:
    """Split `name#count@name#count` into `[name, count]` pairs; a name `None` becomes null.

    An entry ends at the first `#<digits>` followed by `@` or the end, because some names hold `@` or `#`, such as
    `Aspartyl-tRNA(Asn) amidotransferase subunit B @ Glutamyl-tRNA(Gln) amidotransferase subunit B#25`.
    """
    if text is None or text in ("", NONE):
        return []
    text = str(text)
    pairs, position = [], 0
    while position < len(text):
        entry = _PAIR.match(text, position)
        if entry is None:
            raise ValueError(f"'{text[position:]}' of '{text}' is not '<name>#<count>'")
        pairs.append([_name(entry.group(1)), int(entry.group(2))])
        position = entry.end()
    return pairs


def _integer(value: Any, field: str) -> int:
    if isinstance(value, int) and not isinstance(value, bool):
        return value
    if isinstance(value, float) and value.is_integer():
        return int(value)
    if isinstance(value, str) and _INTEGER.fullmatch(value.strip()):
        return int(value)
    raise ValueError(f"{field} {value!r} is no integer")


def _number(value: Any, field: str) -> float | int:
    if isinstance(value, int | float) and not isinstance(value, bool) and math.isfinite(value):
        return value
    if isinstance(value, str) and (number := plain_number(value.strip())) is not None:
        return number
    raise ValueError(f"{field} {value!r} is no finite number")


def _optional_number(value: Any) -> float | int | None:
    return None if value is None or value in ("", NONE) else _number(value, "divers")


def _branch_length(node: Mapping[str, Any]) -> float | int:
    length = _number(node["branch_length"], "branch_length")
    if length < 0:
        raise ValueError(f"node '{node['name']}' has the negative branch length {length}")
    return length


def _extra_value(value: Any, key: str) -> Any:
    if value is None or isinstance(value, str) or (isinstance(value, int | float) and not isinstance(value, bool)):
        return value
    raise ValueError(f"cluster column '{key}' has the value {value!r}, which is no number or string")


def _join_json(parent: dict[str, Any], child: dict[str, Any]) -> dict[str, Any]:
    joined = dict(child)
    if "branch_length" in parent or "branch_length" in child:
        joined["branch_length"] = _add(parent.get("branch_length", 0), child.get("branch_length", 0))
    for key in ("muts", "aa_muts"):
        if key in parent or key in child:
            joined[key] = ",".join(part for part in (parent.get(key), child.get(key)) if part)
    return joined


def _add(a: float, b: float) -> float:
    """Return the sum of the decimal readings of `a` and `b`, so that `0.1` and `0.2` give `0.3`."""
    return float(Decimal(repr(a)) + Decimal(repr(b)))


@dataclass(frozen=True)
class _NewickNode:
    label: str
    children: list["_NewickNode"] | None


def _parse_newick(text: str) -> tuple[_NewickNode, str]:
    """Parse a Newick text into nodes whose labels keep name, support, and length as written; return the text from
    the final `;` on. A label may be quoted with `'`, with `''` or `\\'` for a quote inside it."""
    end = text.rindex(";")
    stack: list[list[_NewickNode]] = []
    pending: _NewickNode | None = None
    position = 0
    while position < end:
        character = text[position]
        if character == "(":
            stack.append([])
            position += 1
            continue
        if character in ",)":
            if not stack:
                raise ValueError(f"unbalanced Newick text at offset {position}")
            stack[-1].append(pending or _NewickNode("", None))
            pending = None
            position += 1
            if character == ",":
                continue
            label_end = _label_end(text, position, end)
            pending = _NewickNode(text[position:label_end], stack.pop())
            position = label_end
            continue
        label_end = _label_end(text, position, end)
        if label_end == position:
            raise ValueError(f"unexpected '{character}' in Newick text at offset {position}")
        pending = _NewickNode(text[position:label_end], None)
        position = label_end
    if stack or pending is None:
        raise ValueError("unbalanced Newick text")
    return pending, text[end:]


def _label_end(text: str, position: int, end: int) -> int:
    """Return the offset after a Newick label that starts at `position`, with its quoted parts."""
    while position < end:
        character = text[position]
        if character in "(),;":
            return position
        position += 1
        if character != "'":
            continue
        while position < end:
            if text[position] == "\\" or text.startswith("''", position):
                position += 2
            elif text[position] == "'":
                position += 1
                break
            else:
                position += 1
    return position


def _write_newick(node: _NewickNode) -> str:
    if node.children is None:
        return node.label
    return "(" + ",".join(map(_write_newick, node.children)) + ")" + node.label


def _join_newick(parent: _NewickNode, child: _NewickNode) -> _NewickNode:
    _, parent_length = _newick_length(parent.label)
    name, child_length = _newick_length(child.label)
    if parent_length is None:
        return child
    total = parent_length + (child_length or Decimal(0))
    return _NewickNode(f"{name}:{total}", child.children)


def _tip_accession(label: str) -> str:
    """Return the accession of a gene-tree tip label: the label without quotes up to the first `|` or `:`."""
    return re.split(r"[|:]", label.removeprefix("'"), maxsplit=1)[0]


def _newick_length(label: str) -> tuple[str, Decimal | None]:
    name, separator, length = label.rpartition(":")
    if not separator:
        return label, None
    try:
        return name, Decimal(length)
    except InvalidOperation:
        return label, None


def newick_tips(text: str) -> list[str]:
    """Return the tip labels of a Newick text, without their branch lengths."""
    root, _ = _parse_newick(text)
    tips: list[str] = []
    stack = [root]
    while stack:
        node = stack.pop()
        if node.children is None:
            tips.append(_newick_length(node.label)[0])
        else:
            stack.extend(reversed(node.children))
    return tips
