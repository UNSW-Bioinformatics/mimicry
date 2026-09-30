#!/usr/bin/env python3
"""Join modified match tables to high-binding HLA prediction results.

Input match-table filenames
---------------------------

By default, the script processes only files ending exactly in::

    _match_tables0.csv
    _match_tables1.csv

Examples from a split dataset are::

    ligandmhc2_5001-6000_match_tables0.csv
    ligandmhc2_5001-6000_match_tables1.csv
    ligandmhc2_9001-10000_match_tables0.csv
    ligandmhc2_9001-10000_match_tables1.csv

Every input remains a separate output, so numbered chunks are never combined.
Table 0 and table 1 also remain separate. For example::

    ligandmhc2_9001-10000_match_tables0.csv
        -> ligandmhc2_9001-10000_match_tables0_with_hla.csv

    ligandmhc2_9001-10000_match_tables1.csv
        -> ligandmhc2_9001-10000_match_tables1_with_hla.csv

HLA-file selection
------------------

Normally, the script automatically pairs each dataset family with a file ending
in ``_predictions_long.csv`` under ``--hla-dir``. Project aliases such as
``ligandmhc2 -> taa_mhc2lig``, ``ise_* -> taa_ise``, ``mizukoshi -> taa_mizu``
and ``reparaz/reparatz -> taa_repa`` are supported.

You can instead state the exact HLA file with ``--hla-file``. This is useful for
one family of split inputs, for example::

    python3 combine_match_tables01_with_high_binding_hla.py --match-dir MATCH_DIR --match-prefix ligandmhc2 --hla-dir HLA_DIR --hla-file taa_mhc2lig_predictions_long.csv --output-dir OUTPUT_DIR

A relative ``--hla-file`` is resolved relative to ``--hla-dir``; an absolute
path is also accepted. ``--match-prefix`` can be repeated and limits which
modified match-table files are selected.

Join behaviour
--------------

The peptide column is detected automatically (or can be supplied with
``--match-peptide-column``). Matching is case-insensitive and ignores peptide
whitespace. The join is always an INNER JOIN: only peptides present in the
selected HLA-results file are written. Therefore, when ``hla_combine_results``
contains only high-binding peptides, every output peptide is high binding.

If a peptide has multiple distinct HLA prediction rows, the match-table row is
repeated once per HLA row. Exact duplicate HLA rows are removed. HLA columns
are prefixed with ``hla_`` and provenance columns retain source filenames and
CSV record numbers.

Large files are streamed and relevant HLA rows are indexed in a temporary
SQLite database. Existing outputs are skipped by default; use ``--overwrite``
to recreate them.
"""

from __future__ import print_function

import argparse
import csv
import hashlib
import json
import os
import re
import sqlite3
import sys
import tempfile
from collections import defaultdict
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Sequence, Set, Tuple


DEFAULT_MATCH_DIR = Path("match_tables")
DEFAULT_HLA_DIR = Path("hla_combine_results")
DEFAULT_OUTPUT_DIR = Path("match_tables01_with_high_binding_hla")

HLA_SUFFIX = "_predictions_long.csv"
MATCH_FILE_RE = re.compile(
    r"^(?P<dataset>.+)_match_tables(?P<table_number>[01]).*\.csv$",
    re.IGNORECASE,
)

# Project-specific filename aliases used only when pairing files. Names are
# normalised, a leading taa_ is ignored, and a terminal numeric chunk range is
# removed before these aliases are applied.
#
# Some match-table datasets use a study/author name while the corresponding HLA
# file uses a shortened dataset name. These aliases are deliberately applied
# only to filename pairing; they never change the output filename or CSV data.
DATASET_ALIASES = {
    "ligandmhc1": "mhc1lig",
    "ligandmhc2": "mhc2lig",
    "mizukoshi": "mizu",
    "reparaz": "repa",
    "reparatz": "repa",
}

# This pattern also permits ambiguity/non-standard amino-acid symbols sometimes
# present in prediction inputs. Peptide matching itself does not require the
# value to pass this pattern; it is used only for fallback column detection.
PEPTIDE_RE = re.compile(r"^[ACDEFGHIKLMNPQRSTVWYBXZJUO]+$", re.IGNORECASE)
QUERY_ID_RE = re.compile(
    r"^(?P<peptide>[A-Za-z]+)_.+_(?:ClassI|ClassII)_\d+$",
    re.IGNORECASE,
)


@dataclass(frozen=True)
class CsvSpec:
    fieldnames: Tuple[str, ...]
    dialect: object


@dataclass(frozen=True)
class MatchName:
    dataset: str
    table_number: int


@dataclass
class MatchJob:
    match_path: Path
    hla_path: Path
    output_path: Path
    match_name: MatchName
    match_spec: CsvSpec
    match_peptide_index: int


@dataclass
class HlaDatabaseInfo:
    connection: sqlite3.Connection
    database_path: Path
    spec: CsvSpec
    peptide_index: int
    allele_index: Optional[int]
    scanned_rows: int
    relevant_rows: int
    unique_rows: int
    duplicate_rows: int
    blank_peptide_rows: int


def positive_int(value: str) -> int:
    try:
        number = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError("must be an integer")
    if number <= 0:
        raise argparse.ArgumentTypeError("must be greater than zero")
    return number


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--match-dir",
        type=Path,
        default=DEFAULT_MATCH_DIR,
        help=(
            "Directory containing *_match_tables0.csv and "
            "*_match_tables1.csv files. Default: {}"
        ).format(DEFAULT_MATCH_DIR),
    )
    parser.add_argument(
        "--match-prefix",
        action="append",
        default=[],
        metavar="PREFIX",
        help=(
            "Process only modified match-table files whose dataset name starts "
            "with PREFIX. For example, --match-prefix ligandmhc2 selects all "
            "ligandmhc2 numbered chunks and both table 0 and table 1. May be "
            "supplied more than once. Default: process all *.csv match "
            "tables in --match-dir."
        ),
    )
    parser.add_argument(
        "--hla-dir",
        type=Path,
        default=DEFAULT_HLA_DIR,
        help=(
            "Directory containing *_predictions_long.csv files. This can be "
            "hla_combine_results or another HLA-output directory. Default: {}"
        ).format(DEFAULT_HLA_DIR),
    )
    parser.add_argument(
        "--hla-file",
        type=Path,
        help=(
            "Use this exact HLA-results CSV for every selected match-table "
            "file, instead of automatic filename pairing. A relative path is "
            "resolved under --hla-dir; an absolute path is accepted. Combine "
            "this with --match-prefix when --match-dir contains multiple "
            "dataset families."
        ),
    )
    parser.add_argument(
        "--allow-hla-file-for-all",
        action="store_true",
        help=(
            "Allow one explicit --hla-file to be applied to more than one "
            "canonical match-table dataset family. Without this option, the "
            "script stops to prevent an accidental cross-dataset join."
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=(
            "Output directory. Each modified input remains a separate output, "
            "with match_tables0_mod or match_tables1_mod retained in its filename. "
            "Default: {}"
        ).format(DEFAULT_OUTPUT_DIR),
    )
    parser.add_argument(
        "--recursive-hla",
        action="store_true",
        help=(
            "Search recursively below --hla-dir for *_predictions_long.csv. "
            "Useful if the HLA files are inside dataset subdirectories."
        ),
    )
    parser.add_argument(
        "--match-peptide-column",
        help=(
            "Peptide column in the match-table files. May be a header name or "
            "a 1-based column number. Default: detect automatically; the "
            "sample files resolve to 'Peptide Sequence'."
        ),
    )
    parser.add_argument(
        "--hla-peptide-column",
        help=(
            "Peptide column in the HLA prediction files. May be a header name "
            "or a 1-based column number. Default: detect automatically."
        ),
    )
    parser.add_argument(
        "--hla-allele-column",
        help=(
            "Optional HLA allele column, used only for deterministic output "
            "ordering. May be a header name or 1-based column number."
        ),
    )
    parser.add_argument(
        "--join",
        choices=("inner",),
        default="inner",
        help=(
            "Compatibility option. Only 'inner' is allowed: output rows are "
            "written only for peptides present in the paired HLA-results file."
        ),
    )
    parser.add_argument(
        "--file-map",
        action="append",
        default=[],
        metavar="MATCH_DATASET=HLA_DATASET",
        help=(
            "Manual filename-family mapping, excluding the standard suffixes. "
            "For example --file-map ligandmhc1=taa_mhc1lig. A base mapping "
            "automatically applies to every numbered chunk and to table 0/1. "
            "May be supplied more than once."
        ),
    )
    parser.add_argument(
        "--temp-dir",
        type=Path,
        help=(
            "Directory for temporary SQLite indexes. Default: --output-dir. "
            "A local scratch disk may be faster."
        ),
    )
    parser.add_argument(
        "--batch-size",
        type=positive_int,
        default=10000,
        help="SQLite insertion batch size. Default: 10000.",
    )
    parser.add_argument(
        "--cache-peptides",
        type=positive_int,
        default=5000,
        help="Number of peptide lookups cached per output file. Default: 5000.",
    )
    parser.add_argument(
        "--progress-every",
        type=positive_int,
        default=1000000,
        help="Print progress after every N rows. Default: 1000000.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help=(
            "Overwrite output files that already exist. By default an existing "
            "output is skipped before its CSV data are read."
        ),
    )
    parser.add_argument(
        "--strict-file-matching",
        action="store_true",
        help=(
            "Stop if any match-table file cannot be paired uniquely. Default: "
            "report and skip unpaired/ambiguous files."
        ),
    )
    parser.add_argument(
        "--strict-peptides",
        action="store_true",
        help=(
            "Stop at the first blank peptide in either a match-table or a "
            "relevant HLA row. Default: retain/report blank match peptides and "
            "skip blank HLA peptide rows."
        ),
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show filename pairings and output names, then exit.",
    )
    return parser.parse_args()


def configure_csv_field_limit() -> None:
    """Set the largest CSV field size accepted by this Python build."""
    limit = sys.maxsize
    while True:
        try:
            csv.field_size_limit(limit)
            return
        except OverflowError:
            limit //= 10


def normalise_header(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", value.strip().lower()).strip("_")


def normalise_dataset_name(value: str) -> str:
    value = Path(value).name.lower().strip()
    value = re.sub(r"[^a-z0-9]+", "_", value)
    return re.sub(r"_+", "_", value).strip("_")


def without_taa_prefix(value: str) -> str:
    value = normalise_dataset_name(value)
    return value[4:] if value.startswith("taa_") else value


def strip_terminal_chunk_range(value: str) -> str:
    """Remove a terminal numeric chunk label for filename pairing only."""
    value = normalise_dataset_name(value)
    # Handles 7001-8000, 7001_8000, 7001-to-8000, and 7001to8000 after
    # normalisation. The original filename/output name is never changed.
    value = re.sub(r"_\d+_(?:to_)?\d+$", "", value)
    value = re.sub(r"_\d+to\d+$", "", value)
    return value


def canonical_dataset_name(value: str) -> str:
    """Return the project-level dataset family used for filename pairing.

    Exact filename matches are scored before this family-level comparison, so
    a more specific HLA file still wins when one exists. The family rules below
    provide fallbacks for known naming differences in this project:

    * ise_aml / ise_ball / ise_tall / other ise_* chunks -> taa_ise
    * any ise_* dataset containing ``extra``             -> taa_ise_extra
    * mizukoshi                                          -> taa_mizu
    * reparaz or reparatz                                -> taa_repa
    * ligandmhc1 / ligandmhc2                            -> taa_mhc1lig / taa_mhc2lig

    The leading ``taa_`` is removed from HLA filenames before comparison.
    """
    value = without_taa_prefix(value)
    value = strip_terminal_chunk_range(value)
    value = re.sub(r"_all$", "", value)

    # ISE match tables are divided into disease-specific and numbered chunks,
    # whereas the HLA predictions are stored as taa_ise_predictions_long.csv.
    # Extra ISE records use taa_ise_extra_predictions_long.csv.
    if value == "ise" or value.startswith("ise_"):
        tokens = set(value.split("_"))
        return "ise_extra" if "extra" in tokens else "ise"


    # Author/study names may include an additional suffix; all such chunks use
    # the same abbreviated HLA prediction family.
    if value == "mizukoshi" or value.startswith("mizukoshi_"):
        return "mizu"
    if (
        value == "reparaz"
        or value.startswith("reparaz_")
        or value == "reparatz"
        or value.startswith("reparatz_")
    ):
        return "repa"

    return DATASET_ALIASES.get(value, value)


def parse_match_filename(path: Path) -> Optional[MatchName]:
    match = MATCH_FILE_RE.fullmatch(path.name)
    if match is None:
        return None
    return MatchName(
        dataset=match.group("dataset"),
        table_number=int(match.group("table_number")),
    )


def dataset_name_from_hla(path: Path) -> str:
    """Return the dataset portion of an HLA filename.

    Automatically discovered files normally end in ``_predictions_long.csv``.
    An explicitly supplied ``--hla-file`` may have another CSV filename, so in
    that case its stem is used rather than rejecting it.
    """
    if path.name.lower().endswith(HLA_SUFFIX):
        return path.name[: -len(HLA_SUFFIX)]
    return path.stem


def output_path_for(match_path: Path, output_dir: Path) -> Path:
    # Preserve the complete dataset/chunk name and the table number.
    return output_dir / "{}_with_hla.csv".format(match_path.stem)


def parse_manual_mappings(values: Sequence[str]) -> Dict[str, str]:
    mappings: Dict[str, str] = {}
    for value in values:
        if "=" not in value:
            raise SystemExit(
                "Invalid --file-map {!r}; expected MATCH_DATASET=HLA_DATASET".format(
                    value
                )
            )
        left, right = value.split("=", 1)
        left_key = normalise_dataset_name(left)
        right_key = normalise_dataset_name(right)
        if not left_key or not right_key:
            raise SystemExit(
                "Invalid --file-map {!r}; neither side may be blank".format(value)
            )
        previous = mappings.get(left_key)
        if previous is not None and previous != right_key:
            raise SystemExit(
                "Conflicting --file-map entries for {!r}".format(left)
            )
        mappings[left_key] = right_key
    return mappings


def pairing_score(match_name: str, hla_name: str) -> Optional[int]:
    """Return a filename-pairing score, with higher scores preferred."""
    match_norm = normalise_dataset_name(match_name)
    hla_norm = normalise_dataset_name(hla_name)
    match_base = strip_terminal_chunk_range(match_norm)
    hla_base = strip_terminal_chunk_range(hla_norm)

    if match_norm == hla_norm:
        return 500
    if match_base == hla_norm or match_base == hla_base:
        return 450
    if without_taa_prefix(match_norm) == without_taa_prefix(hla_norm):
        return 400
    if without_taa_prefix(match_base) == without_taa_prefix(hla_base):
        return 350
    if canonical_dataset_name(match_norm) == canonical_dataset_name(hla_norm):
        return 300
    return None


def discover_hla_files(hla_dir: Path, recursive: bool) -> List[Path]:
    pattern = "*{}".format(HLA_SUFFIX)
    iterator = hla_dir.rglob(pattern) if recursive else hla_dir.glob(pattern)
    return sorted(path for path in iterator if path.is_file())


def match_selected_by_prefix(
    path: Path,
    match_name: MatchName,
    prefixes: Sequence[str],
) -> bool:
    """Return True when a match-table file passes --match-prefix filters."""
    if not prefixes:
        return True

    dataset_key = normalise_dataset_name(match_name.dataset)
    filename_key = normalise_dataset_name(path.name)

    for raw_prefix in prefixes:
        prefix = normalise_dataset_name(raw_prefix)
        if not prefix:
            raise SystemExit("--match-prefix may not be blank")
        if (
            dataset_key == prefix
            or dataset_key.startswith(prefix + "_")
            or filename_key == prefix
            or filename_key.startswith(prefix + "_")
        ):
            return True
    return False


def resolve_explicit_hla_file(value: Path, hla_dir: Path) -> Path:
    """Resolve --hla-file as an absolute path or relative to --hla-dir."""
    candidates: List[Path]
    if value.is_absolute():
        candidates = [value]
    else:
        candidates = [hla_dir / value, value]

    checked: List[Path] = []
    for candidate in candidates:
        candidate = candidate.expanduser()
        if candidate in checked:
            continue
        checked.append(candidate)
        if candidate.is_file():
            return candidate.resolve()

    raise SystemExit(
        "The explicit --hla-file was not found. Checked: {}".format(
            ", ".join(str(path) for path in checked)
        )
    )


def discover_pairs(
    match_dir: Path,
    hla_dir: Path,
    recursive_hla: bool,
    manual_mappings: Dict[str, str],
    strict: bool,
    match_prefixes: Sequence[str],
    explicit_hla_path: Optional[Path],
    allow_hla_file_for_all: bool,
) -> Tuple[List[Tuple[Path, MatchName, Path]], List[str]]:
    """Discover modified match-table/HLA pairs.

    When ``explicit_hla_path`` is supplied, it is paired with every selected
    modified match-table file. Otherwise the normal filename/alias matching
    rules are used.
    """
    match_files: List[Tuple[Path, MatchName]] = []
    for path in sorted(match_dir.glob("*.csv")):
        if not path.is_file():
            continue
        parsed = parse_match_filename(path)
        if parsed is None:
            continue
        if match_selected_by_prefix(path, parsed, match_prefixes):
            match_files.append((path, parsed))

    if not match_files:
        prefix_text = (
            " after applying --match-prefix {}".format(
                ", ".join(repr(value) for value in match_prefixes)
            )
            if match_prefixes
            else ""
        )
        raise SystemExit(
            "No *_match_tables0.csv or *_match_tables1.csv files "
            "were found in {}{}".format(match_dir, prefix_text)
        )

    if explicit_hla_path is not None:
        families = sorted(
            {canonical_dataset_name(match_name.dataset) for _, match_name in match_files}
        )
        if len(families) > 1 and not allow_hla_file_for_all:
            raise SystemExit(
                "The selected modified match tables cover more than one dataset "
                "family ({}), but one --hla-file was supplied. Add one or more "
                "--match-prefix options to select a single family, or add "
                "--allow-hla-file-for-all if applying the same HLA file to all "
                "of them is intentional.".format(", ".join(families))
            )
        return [
            (match_path, match_name, explicit_hla_path)
            for match_path, match_name in match_files
        ], []

    hla_files = discover_hla_files(hla_dir, recursive_hla)
    if not hla_files:
        location_text = "recursively below" if recursive_hla else "in"
        raise SystemExit(
            "No *{} files found {} {}".format(
                HLA_SUFFIX, location_text, hla_dir
            )
        )

    hla_by_normalised_name: Dict[str, List[Path]] = defaultdict(list)
    for path in hla_files:
        key = normalise_dataset_name(dataset_name_from_hla(path))
        hla_by_normalised_name[key].append(path)

    pairs: List[Tuple[Path, MatchName, Path]] = []
    messages: List[str] = []

    for match_path, match_name in match_files:
        dataset_key = normalise_dataset_name(match_name.dataset)
        base_key = strip_terminal_chunk_range(dataset_key)

        mapping_key = next(
            (key for key in (dataset_key, base_key) if key in manual_mappings),
            None,
        )

        if mapping_key is not None:
            target_key = manual_mappings[mapping_key]
            candidates = list(hla_by_normalised_name.get(target_key, []))
            if not candidates:
                candidates = [
                    path
                    for path in hla_files
                    if without_taa_prefix(dataset_name_from_hla(path))
                    == without_taa_prefix(target_key)
                ]

            if len(candidates) == 1:
                pairs.append((match_path, match_name, candidates[0]))
                continue

            if not candidates:
                message = (
                    "UNPAIRED {}: manual mapping {!r} -> {!r}, but no "
                    "matching HLA file exists"
                ).format(match_path.name, mapping_key, target_key)
            else:
                message = (
                    "AMBIGUOUS {}: manual mapping {!r} -> {!r} matches: {}"
                ).format(
                    match_path.name,
                    mapping_key,
                    target_key,
                    ", ".join(str(path) for path in candidates),
                )
            if strict:
                raise SystemExit(message)
            messages.append(message)
            continue

        scored: List[Tuple[int, Path]] = []
        for hla_path in hla_files:
            score = pairing_score(
                match_name.dataset,
                dataset_name_from_hla(hla_path),
            )
            if score is not None:
                scored.append((score, hla_path))

        if not scored:
            message = "UNPAIRED {}: no compatible HLA filename".format(
                match_path.name
            )
            if strict:
                raise SystemExit(message)
            messages.append(message)
            continue

        best_score = max(score for score, _ in scored)
        best = [path for score, path in scored if score == best_score]
        if len(best) != 1:
            message = "AMBIGUOUS {}: equally good HLA files: {}".format(
                match_path.name,
                ", ".join(str(path) for path in best),
            )
            if strict:
                raise SystemExit(message)
            messages.append(message)
            continue

        pairs.append((match_path, match_name, best[0]))

    return pairs, messages

def make_unique_headers(raw_headers: Sequence[str]) -> Tuple[str, ...]:
    headers: List[str] = []
    used: Set[str] = set()

    for position, raw in enumerate(raw_headers, start=1):
        base = raw.strip()
        if not base:
            base = "unnamed_column_{}".format(position)

        candidate = base
        counter = 2
        while candidate in used:
            candidate = "{}__{}".format(base, counter)
            counter += 1

        headers.append(candidate)
        used.add(candidate)

    return tuple(headers)


def sniff_dialect(path: Path) -> object:
    """Detect the delimiter while enforcing normal RFC-style CSV quoting.

    ``csv.Sniffer`` can incorrectly infer ``doublequote=False`` when a field
    contains an escaped literal quote such as::

        three double quotes + In vitro T-cell analysis , flow cytometry
        + three double quotes

    With that incorrect setting, the comma inside the quoted text is treated
    as a delimiter and the row appears to have one extra column.  The project
    files use conventional double-quote escaping, so we use Sniffer only to
    identify the delimiter and force ``quotechar='\"'`` and
    ``doublequote=True``.
    """
    with path.open(
        "r", encoding="utf-8-sig", errors="replace", newline=""
    ) as handle:
        sample = handle.read(131072)

    if not sample:
        raise ValueError("CSV is empty: {}".format(path))

    try:
        sniffed = csv.Sniffer().sniff(sample, delimiters=",\t;")
        delimiter = sniffed.delimiter
        skipinitialspace = sniffed.skipinitialspace
    except csv.Error:
        delimiter = ","
        skipinitialspace = False

    # csv.reader accepts a Dialect class.  Constructing one explicitly avoids
    # retaining unreliable quote/doublequote guesses made by csv.Sniffer.
    return type(
        "ProjectCsvDialect",
        (csv.Dialect,),
        {
            "delimiter": delimiter,
            "quotechar": '"',
            "escapechar": None,
            "doublequote": True,
            "skipinitialspace": skipinitialspace,
            "lineterminator": "\n",
            "quoting": csv.QUOTE_MINIMAL,
            "strict": False,
        },
    )


def read_csv_spec(path: Path) -> CsvSpec:
    dialect = sniff_dialect(path)
    with path.open(
        "r", encoding="utf-8-sig", errors="replace", newline=""
    ) as handle:
        reader = csv.reader(handle, dialect=dialect)
        raw_header = next(reader, None)

    if raw_header is None:
        raise ValueError("CSV has no header row: {}".format(path))
    if len(raw_header) == 0:
        raise ValueError("CSV header is empty: {}".format(path))

    return CsvSpec(make_unique_headers(raw_header), dialect)


def iter_csv_rows(path: Path, spec: CsvSpec) -> Iterator[Tuple[int, List[str]]]:
    """Yield (1-based CSV record number, values); the header is record 1."""
    expected = len(spec.fieldnames)

    with path.open(
        "r", encoding="utf-8-sig", errors="replace", newline=""
    ) as handle:
        reader = csv.reader(handle, dialect=spec.dialect)
        next(reader, None)

        for record_number, values in enumerate(reader, start=2):
            if len(values) < expected:
                values = values + [""] * (expected - len(values))
            elif len(values) > expected:
                raise ValueError(
                    "{} record {} has {} values but the header has {} columns".format(
                        path,
                        record_number,
                        len(values),
                        expected,
                    )
                )
            yield record_number, values


def normalise_peptide_key(value: str) -> str:
    """Normalise a peptide sequence for joining."""
    text = "".join(str(value).strip().split())
    if text.startswith(">"):
        text = text[1:]

    query_match = QUERY_ID_RE.fullmatch(text)
    if query_match is not None:
        text = query_match.group("peptide")

    return text.upper()


def explicit_column_index(
    fieldnames: Sequence[str], value: str, role: str
) -> int:
    value = value.strip()

    if value.isdigit():
        index = int(value) - 1
        if not 0 <= index < len(fieldnames):
            raise ValueError(
                "{} column number {} is outside 1..{}".format(
                    role, value, len(fieldnames)
                )
            )
        return index

    exact = [index for index, name in enumerate(fieldnames) if name == value]
    if len(exact) == 1:
        return exact[0]

    target = normalise_header(value)
    normalised = [
        index
        for index, name in enumerate(fieldnames)
        if normalise_header(name) == target
    ]
    if len(normalised) == 1:
        return normalised[0]

    raise ValueError(
        "Could not identify {} column {!r}. Available columns: {}".format(
            role,
            value,
            ", ".join(fieldnames),
        )
    )


def first_named_column(
    fieldnames: Sequence[str], candidates: Sequence[str]
) -> Optional[int]:
    by_normalised: Dict[str, List[int]] = defaultdict(list)
    for index, name in enumerate(fieldnames):
        by_normalised[normalise_header(name)].append(index)

    for candidate in candidates:
        indexes = by_normalised.get(normalise_header(candidate), [])
        if len(indexes) == 1:
            return indexes[0]
    return None


def is_plausible_peptide(value: str) -> bool:
    peptide = normalise_peptide_key(value)
    return 5 <= len(peptide) <= 100 and PEPTIDE_RE.fullmatch(peptide) is not None


def detect_peptide_column(
    path: Path,
    spec: CsvSpec,
    explicit: Optional[str],
    role: str,
) -> int:
    if explicit:
        return explicit_column_index(spec.fieldnames, explicit, role)

    if role == "match-table":
        direct_names = (
            "Peptide Sequence",
            "peptide_sequence",
            "peptide sequence",
            "peptide_seq",
            "query_peptide",
            "query sequence",
            "query_sequence",
            "peptide",
            "qseq",
        )
    else:
        direct_names = (
            "peptide",
            "Peptide",
            "peptide_sequence",
            "peptide sequence",
            "peptide_seq",
            "query_peptide",
            "query_sequence",
            "query sequence",
            "source_peptide",
            "tumour_peptide",
            "tumor_peptide",
            "qseq",
        )

    named = first_named_column(spec.fieldnames, direct_names)
    if named is not None:
        return named

    # Header-based detection is expected to resolve the project files. This
    # fallback samples up to 250 records if an unusual header is encountered.
    nonblank = [0] * len(spec.fieldnames)
    peptide_hits = [0] * len(spec.fieldnames)

    for sampled, (_, values) in enumerate(iter_csv_rows(path, spec), start=1):
        for index, value in enumerate(values):
            if not str(value).strip():
                continue
            nonblank[index] += 1
            if is_plausible_peptide(value):
                peptide_hits[index] += 1
        if sampled >= 250:
            break

    candidates: List[Tuple[float, int]] = []
    for index, count in enumerate(nonblank):
        if count == 0:
            continue
        ratio = peptide_hits[index] / float(count)
        if ratio < 0.80:
            continue

        header = normalise_header(spec.fieldnames[index])
        if any(
            token in header
            for token in (
                "subject",
                "sseq",
                "accession",
                "title",
                "protein_name",
                "gene_name",
                "organism",
            )
        ):
            continue
        candidates.append((ratio, index))

    if candidates:
        best_ratio = max(ratio for ratio, _ in candidates)
        best = [index for ratio, index in candidates if ratio == best_ratio]
        if len(best) == 1:
            return best[0]

    option = "--match-peptide-column" if role == "match-table" else "--hla-peptide-column"
    raise ValueError(
        "Could not unambiguously detect the {} peptide column in {}. "
        "Available columns: {}. Supply {} explicitly.".format(
            role,
            path,
            ", ".join(spec.fieldnames),
            option,
        )
    )


def detect_allele_column(
    spec: CsvSpec,
    explicit: Optional[str],
) -> Optional[int]:
    if explicit:
        return explicit_column_index(spec.fieldnames, explicit, "HLA allele")

    return first_named_column(
        spec.fieldnames,
        (
            "allele",
            "hla_allele",
            "hla allele",
            "mhcflurry_allele",
            "netmhcpan_allele",
            "mhc",
            "hla",
        ),
    )


def collect_needed_peptides(
    jobs: Sequence[MatchJob],
    progress_every: int,
    strict_peptides: bool,
) -> Tuple[Set[str], Dict[Path, Dict[str, int]]]:
    needed: Set[str] = set()
    stats: Dict[Path, Dict[str, int]] = {}

    for job in jobs:
        file_stats = {
            "rows": 0,
            "blank_peptide_rows": 0,
            "new_unique_peptides": 0,
        }
        before = len(needed)

        print("  Reading peptide keys from {}".format(job.match_path.name))
        for record_number, values in iter_csv_rows(job.match_path, job.match_spec):
            file_stats["rows"] += 1
            peptide = normalise_peptide_key(values[job.match_peptide_index])
            if not peptide:
                file_stats["blank_peptide_rows"] += 1
                if strict_peptides:
                    raise ValueError(
                        "{} record {} has a blank peptide in column {!r}".format(
                            job.match_path,
                            record_number,
                            job.match_spec.fieldnames[job.match_peptide_index],
                        )
                    )
            else:
                needed.add(peptide)

            if file_stats["rows"] % progress_every == 0:
                print(
                    "    {:,} rows scanned; {:,} unique peptide keys so far".format(
                        file_stats["rows"], len(needed)
                    )
                )

        file_stats["new_unique_peptides"] = len(needed) - before
        stats[job.match_path] = file_stats
        print(
            "    {:,} rows; {:,} newly added unique peptides; {:,} blank "
            "peptide rows".format(
                file_stats["rows"],
                file_stats["new_unique_peptides"],
                file_stats["blank_peptide_rows"],
            )
        )

    return needed, stats


def create_temp_database(temp_dir: Path, stem: str) -> Tuple[sqlite3.Connection, Path]:
    temp_dir.mkdir(parents=True, exist_ok=True)
    handle = tempfile.NamedTemporaryFile(
        prefix="{}_".format(normalise_dataset_name(stem) or "hla"),
        suffix=".sqlite3",
        dir=str(temp_dir),
        delete=False,
    )
    database_path = Path(handle.name)
    handle.close()

    connection = sqlite3.connect(str(database_path))
    connection.execute("PRAGMA journal_mode=OFF")
    connection.execute("PRAGMA synchronous=OFF")
    connection.execute("PRAGMA temp_store=FILE")
    connection.execute("PRAGMA cache_size=-131072")  # approximately 128 MiB
    connection.execute("PRAGMA locking_mode=EXCLUSIVE")
    connection.execute(
        """
        CREATE TABLE hla_rows (
            peptide_key TEXT NOT NULL,
            row_hash TEXT NOT NULL,
            allele_key TEXT NOT NULL,
            source_record INTEGER NOT NULL,
            payload TEXT NOT NULL,
            PRIMARY KEY (peptide_key, row_hash)
        ) WITHOUT ROWID
        """
    )
    return connection, database_path


def prepare_hla_database(
    hla_path: Path,
    needed_peptides: Set[str],
    hla_peptide_column: Optional[str],
    hla_allele_column: Optional[str],
    temp_dir: Path,
    batch_size: int,
    progress_every: int,
    strict_peptides: bool,
) -> HlaDatabaseInfo:
    spec = read_csv_spec(hla_path)
    peptide_index = detect_peptide_column(
        hla_path,
        spec,
        hla_peptide_column,
        "HLA",
    )
    allele_index = detect_allele_column(spec, hla_allele_column)

    connection, database_path = create_temp_database(
        temp_dir,
        dataset_name_from_hla(hla_path),
    )

    scanned = 0
    relevant = 0
    inserted = 0
    blank_peptide = 0
    batch: List[Tuple[str, str, str, int, str]] = []

    insert_sql = (
        "INSERT OR IGNORE INTO hla_rows "
        "(peptide_key, row_hash, allele_key, source_record, payload) "
        "VALUES (?, ?, ?, ?, ?)"
    )

    print("  Indexing relevant HLA rows from {}".format(hla_path))
    print(
        "    HLA peptide column: {!r} (column {})".format(
            spec.fieldnames[peptide_index],
            peptide_index + 1,
        )
    )
    if allele_index is not None:
        print(
            "    HLA allele column: {!r} (column {})".format(
                spec.fieldnames[allele_index],
                allele_index + 1,
            )
        )
    else:
        print("    No HLA allele column detected; source-record order will be used")

    try:
        connection.execute("BEGIN")

        if needed_peptides:
            for source_record, values in iter_csv_rows(hla_path, spec):
                scanned += 1
                peptide = normalise_peptide_key(values[peptide_index])
                if not peptide:
                    blank_peptide += 1
                    if strict_peptides:
                        raise ValueError(
                            "{} record {} has a blank HLA peptide in column {!r}".format(
                                hla_path,
                                source_record,
                                spec.fieldnames[peptide_index],
                            )
                        )
                    continue

                if peptide not in needed_peptides:
                    if scanned % progress_every == 0:
                        print(
                            "    {:,} HLA rows scanned; {:,} relevant; {:,} "
                            "unique indexed".format(scanned, relevant, inserted)
                        )
                    continue

                relevant += 1
                allele_key = (
                    values[allele_index].strip().upper()
                    if allele_index is not None
                    else ""
                )
                payload = json.dumps(
                    values,
                    ensure_ascii=False,
                    separators=(",", ":"),
                )
                row_hash = hashlib.blake2b(
                    payload.encode("utf-8", errors="replace"),
                    digest_size=16,
                ).hexdigest()
                batch.append(
                    (peptide, row_hash, allele_key, source_record, payload)
                )

                if len(batch) >= batch_size:
                    before = connection.total_changes
                    connection.executemany(insert_sql, batch)
                    inserted += connection.total_changes - before
                    batch = []

                if scanned % progress_every == 0:
                    print(
                        "    {:,} HLA rows scanned; {:,} relevant; {:,} "
                        "unique indexed".format(scanned, relevant, inserted)
                    )
        else:
            print("    No nonblank match-table peptide keys; HLA data scan skipped")

        if batch:
            before = connection.total_changes
            connection.executemany(insert_sql, batch)
            inserted += connection.total_changes - before

        connection.commit()
        connection.execute(
            "CREATE INDEX hla_rows_lookup_order "
            "ON hla_rows(peptide_key, allele_key, source_record)"
        )
        connection.commit()

    except Exception:
        connection.close()
        try:
            database_path.unlink()
        except OSError:
            pass
        raise

    duplicates = relevant - inserted
    print(
        "    Finished: {:,} HLA rows scanned; {:,} relevant; {:,} unique HLA "
        "rows; {:,} exact duplicates removed; {:,} blank HLA peptide rows".format(
            scanned,
            relevant,
            inserted,
            duplicates,
            blank_peptide,
        )
    )

    return HlaDatabaseInfo(
        connection=connection,
        database_path=database_path,
        spec=spec,
        peptide_index=peptide_index,
        allele_index=allele_index,
        scanned_rows=scanned,
        relevant_rows=relevant,
        unique_rows=inserted,
        duplicate_rows=duplicates,
        blank_peptide_rows=blank_peptide,
    )


def unique_output_name(base: str, used: Set[str]) -> str:
    candidate = base
    counter = 2
    while candidate in used:
        candidate = "{}__{}".format(base, counter)
        counter += 1
    used.add(candidate)
    return candidate


def output_header(match_spec: CsvSpec, hla_spec: CsvSpec) -> List[str]:
    used: Set[str] = set()
    header: List[str] = []

    for fixed in ("match_source_file", "match_source_record"):
        header.append(unique_output_name(fixed, used))

    for name in match_spec.fieldnames:
        header.append(unique_output_name(name, used))

    for fixed in ("hla_source_file", "hla_source_record"):
        header.append(unique_output_name(fixed, used))

    for name in hla_spec.fieldnames:
        header.append(unique_output_name("hla_{}".format(name), used))

    return header


def merge_one_job(
    job: MatchJob,
    hla_info: HlaDatabaseInfo,
    cache_size: int,
    progress_every: int,
    strict_peptides: bool,
) -> Dict[str, int]:
    output_path = job.output_path
    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary_output = output_path.with_name(output_path.name + ".tmp")

    header = output_header(job.match_spec, hla_info.spec)
    connection = hla_info.connection

    @lru_cache(maxsize=cache_size)
    def lookup(peptide: str) -> Tuple[Tuple[int, Tuple[str, ...]], ...]:
        rows = connection.execute(
            """
            SELECT source_record, payload
            FROM hla_rows
            WHERE peptide_key = ?
            ORDER BY allele_key, source_record
            """,
            (peptide,),
        ).fetchall()
        return tuple(
            (int(source_record), tuple(json.loads(payload)))
            for source_record, payload in rows
        )

    input_rows = 0
    matched_input_rows = 0
    unmatched_input_rows = 0
    output_rows = 0
    blank_match_peptide_rows = 0

    print("  Writing {} (HLA-matched peptides only)".format(output_path))

    try:
        with temporary_output.open(
            "w", encoding="utf-8", newline=""
        ) as output_handle:
            writer = csv.writer(output_handle, lineterminator="\n")
            writer.writerow(header)

            for match_record, match_values in iter_csv_rows(
                job.match_path,
                job.match_spec,
            ):
                input_rows += 1
                peptide = normalise_peptide_key(
                    match_values[job.match_peptide_index]
                )
                if not peptide:
                    blank_match_peptide_rows += 1
                    if strict_peptides:
                        raise ValueError(
                            "{} record {} has a blank peptide in column {!r}".format(
                                job.match_path,
                                match_record,
                                job.match_spec.fieldnames[job.match_peptide_index],
                            )
                        )
                    hla_rows: Tuple[Tuple[int, Tuple[str, ...]], ...] = ()
                else:
                    hla_rows = lookup(peptide)

                fixed_left = [job.match_path.name, match_record]
                fixed_right = [job.hla_path.name]

                if hla_rows:
                    matched_input_rows += 1
                    for hla_record, hla_values in hla_rows:
                        writer.writerow(
                            fixed_left
                            + list(match_values)
                            + fixed_right
                            + [hla_record]
                            + list(hla_values)
                        )
                        output_rows += 1
                else:
                    # Matched-peptides-only behaviour: do not write this row.
                    # The paired HLA file is the authoritative list of retained
                    # (high-binding) peptides.
                    unmatched_input_rows += 1

                if input_rows % progress_every == 0:
                    print(
                        "    {:,} match-table rows processed; {:,} output rows "
                        "written".format(input_rows, output_rows)
                    )

        os.replace(str(temporary_output), str(output_path))

    except Exception:
        try:
            temporary_output.unlink()
        except OSError:
            pass
        raise

    cache_info = lookup.cache_info()
    print(
        "    Finished: {:,} input rows; {:,} matched input rows; {:,} "
        "unmatched rows filtered out; {:,} output rows".format(
            input_rows,
            matched_input_rows,
            unmatched_input_rows,
            output_rows,
        )
    )
    print(
        "    HLA lookup cache: {:,} hits, {:,} misses".format(
            cache_info.hits,
            cache_info.misses,
        )
    )

    return {
        "input_rows": input_rows,
        "matched_input_rows": matched_input_rows,
        "unmatched_input_rows": unmatched_input_rows,
        "filtered_out_rows": unmatched_input_rows,
        "output_rows": output_rows,
        "blank_match_peptide_rows": blank_match_peptide_rows,
    }


def validate_inputs(args: argparse.Namespace) -> Optional[Path]:
    """Validate directories and resolve an optional exact --hla-file."""
    if not args.match_dir.is_dir():
        raise SystemExit("Match directory does not exist: {}".format(args.match_dir))

    if args.hla_file is not None:
        if args.file_map:
            raise SystemExit(
                "Use either --hla-file or --file-map, not both. --hla-file "
                "already states the exact HLA results file."
            )
        return resolve_explicit_hla_file(args.hla_file, args.hla_dir)

    if not args.hla_dir.is_dir():
        raise SystemExit("HLA directory does not exist: {}".format(args.hla_dir))
    return None


def write_summary(summary_rows: List[Dict[str, object]], output_dir: Path) -> None:
    if not summary_rows:
        return

    summary_path = output_dir / "combine_match_tables01_with_high_binding_hla_summary.csv"
    summary_tmp = summary_path.with_name(summary_path.name + ".tmp")
    fieldnames = list(summary_rows[0].keys())

    with summary_tmp.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(summary_rows)

    os.replace(str(summary_tmp), str(summary_path))
    print("\nSummary written to {}".format(summary_path))


def main() -> int:
    args = parse_args()
    configure_csv_field_limit()
    explicit_hla_path = validate_inputs(args)

    manual_mappings = parse_manual_mappings(args.file_map)
    pairs, pairing_messages = discover_pairs(
        match_dir=args.match_dir,
        hla_dir=args.hla_dir,
        recursive_hla=args.recursive_hla,
        manual_mappings=manual_mappings,
        strict=args.strict_file_matching,
        match_prefixes=args.match_prefix,
        explicit_hla_path=explicit_hla_path,
        allow_hla_file_for_all=args.allow_hla_file_for_all,
    )

    print("Input filename rule: *_match_tables0_mod.csv or *_match_tables1_mod.csv or *_match_tables0.csv or *_match_tables1.csv")
    if args.match_prefix:
        print("Selected match prefix(es): {}".format(", ".join(args.match_prefix)))
    if explicit_hla_path is not None:
        print("Explicit HLA file: {}".format(explicit_hla_path))
    else:
        print("HLA selection: automatic filename/alias pairing under {}".format(args.hla_dir))
    print("Output directory: {}".format(args.output_dir))
    print("Retention rule: keep only peptides present in the paired HLA-results file")
    print("Found {:,} usable input/HLA pair(s):".format(len(pairs)))
    table_counts = {0: 0, 1: 0}
    for match_path, match_name, hla_path in pairs:
        table_counts[match_name.table_number] += 1
        print(
            "  {}  <->  {}  ->  {}".format(
                match_path.name,
                hla_path,
                output_path_for(match_path, args.output_dir).name,
            )
        )
    print(
        "  Pair counts: match_tables0={}, match_tables1={}".format(
            table_counts[0],
            table_counts[1],
        )
    )

    if pairing_messages:
        print("\nFiles not processed:")
        for message in pairing_messages:
            print("  {}".format(message))

    if args.dry_run:
        print("\nDry run complete; no CSV data were processed.")
        return 0

    if not pairs:
        raise SystemExit("No file pairs are available to process.")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    temp_dir = args.temp_dir or args.output_dir
    temp_dir.mkdir(parents=True, exist_ok=True)

    jobs: List[MatchJob] = []
    skipped_existing = 0

    for match_path, match_name, hla_path in pairs:
        output_path = output_path_for(match_path, args.output_dir)
        if output_path.exists() and not args.overwrite:
            print("SKIP existing output: {}".format(output_path))
            skipped_existing += 1
            continue

        match_spec = read_csv_spec(match_path)
        match_peptide_index = detect_peptide_column(
            match_path,
            match_spec,
            args.match_peptide_column,
            "match-table",
        )
        print(
            "Peptide source for {}: {!r} (column {})".format(
                match_path.name,
                match_spec.fieldnames[match_peptide_index],
                match_peptide_index + 1,
            )
        )

        jobs.append(
            MatchJob(
                match_path=match_path,
                hla_path=hla_path,
                output_path=output_path,
                match_name=match_name,
                match_spec=match_spec,
                match_peptide_index=match_peptide_index,
            )
        )

    if not jobs:
        print(
            "Nothing to process: all {:,} paired outputs already exist.".format(
                skipped_existing
            )
        )
        return 0

    jobs_by_hla: Dict[Path, List[MatchJob]] = defaultdict(list)
    for job in jobs:
        jobs_by_hla[job.hla_path].append(job)

    summary_rows: List[Dict[str, object]] = []

    for group_number, hla_path in enumerate(sorted(jobs_by_hla), start=1):
        group_jobs = sorted(
            jobs_by_hla[hla_path],
            key=lambda job: job.match_path.name,
        )
        print(
            "\n[HLA group {}/{}] {} paired with {} match-table file(s)".format(
                group_number,
                len(jobs_by_hla),
                hla_path,
                len(group_jobs),
            )
        )

        needed_peptides, key_stats = collect_needed_peptides(
            group_jobs,
            args.progress_every,
            args.strict_peptides,
        )
        print(
            "  Union contains {:,} unique nonblank peptide key(s)".format(
                len(needed_peptides)
            )
        )

        hla_info: Optional[HlaDatabaseInfo] = None
        try:
            hla_info = prepare_hla_database(
                hla_path=hla_path,
                needed_peptides=needed_peptides,
                hla_peptide_column=args.hla_peptide_column,
                hla_allele_column=args.hla_allele_column,
                temp_dir=temp_dir,
                batch_size=args.batch_size,
                progress_every=args.progress_every,
                strict_peptides=args.strict_peptides,
            )

            del needed_peptides

            for job in group_jobs:
                merge_stats = merge_one_job(
                    job=job,
                    hla_info=hla_info,
                    cache_size=args.cache_peptides,
                    progress_every=args.progress_every,
                    strict_peptides=args.strict_peptides,
                )
                scan_stats = key_stats[job.match_path]
                summary_rows.append(
                    {
                        "match_file": job.match_path.name,
                        "match_dataset": job.match_name.dataset,
                        "match_dataset_base": strip_terminal_chunk_range(
                            job.match_name.dataset
                        ),
                        "match_table_number": job.match_name.table_number,
                        "hla_file": str(job.hla_path),
                        "output_file": str(job.output_path),
                        "join": "inner",
                        "retention_rule": "peptide must occur in paired HLA file",
                        "match_peptide_column": job.match_spec.fieldnames[
                            job.match_peptide_index
                        ],
                        "hla_peptide_column": hla_info.spec.fieldnames[
                            hla_info.peptide_index
                        ],
                        "key_scan_rows": scan_stats["rows"],
                        "key_scan_new_unique_peptides": scan_stats[
                            "new_unique_peptides"
                        ],
                        "key_scan_blank_peptide_rows": scan_stats[
                            "blank_peptide_rows"
                        ],
                        "hla_rows_scanned": hla_info.scanned_rows,
                        "hla_relevant_rows": hla_info.relevant_rows,
                        "hla_unique_rows": hla_info.unique_rows,
                        "hla_duplicate_rows_removed": hla_info.duplicate_rows,
                        "hla_blank_peptide_rows": hla_info.blank_peptide_rows,
                        **merge_stats,
                    }
                )

        finally:
            if hla_info is not None:
                try:
                    hla_info.connection.close()
                finally:
                    try:
                        hla_info.database_path.unlink()
                    except OSError:
                        print(
                            "WARNING: could not remove temporary database {}".format(
                                hla_info.database_path
                            ),
                            file=sys.stderr,
                        )

    write_summary(summary_rows, args.output_dir)

    completed0 = sum(
        1 for row in summary_rows if int(row["match_table_number"]) == 0
    )
    completed1 = sum(
        1 for row in summary_rows if int(row["match_table_number"]) == 1
    )
    print(
        "Completed {:,} separate output file(s) in {}: {:,} from "
        "match_tables0 and {:,} from match_tables1. Skipped {:,} existing "
        "output file(s).".format(
            len(summary_rows),
            args.output_dir,
            completed0,
            completed1,
            skipped_existing,
        )
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        print("\nInterrupted by user.", file=sys.stderr)
        raise SystemExit(130)
    except (OSError, ValueError, csv.Error, sqlite3.Error) as exc:
        print("ERROR: {}".format(exc), file=sys.stderr)
        raise SystemExit(1)
