"""The parameter sheet: what it holds, and the reader that builds it.

    SheetLabels      the exact first-column labels the reader looks for, in one place
    Tile, Sample, IndexGroup
                     small immutable objects, one per tile, sample and group of samples sharing an index
    ParameterSheet   everything tileseq_demux needs from the sheet
    SheetReader      turns a CSV file into a ParameterSheet, or raises SheetError listing every problem
    SheetError       what SheetReader raises

The sheet is the CSV that `tsm csv2json` reads, plus one section that says which samples share a sequencing index:

    ,,,,
    Index Samples,,,,
    Shared Sample ID,Unique Sample ID,,,
    97,193,,,

The shared ID is the sample whose FASTQ files hold both tiles; the unique ID names the second tile's sample. Both
must also be listed in the ordinary sample table, which gives each one its tile. The row of empty cells above the
section is required: without it `tsm csv2json` reads the section as part of the table above.

There are two ways to check something, and the difference is whether the code that follows depends on the result:
    _record_problem_if(condition, message)   only records a problem; nothing below depends on it
    _check(rule_set, subject)                applies ordered rules and returns a Verdict to act on:
                                             stop, skip a row, fall back to a default
Use the first for a single check. Use a RuleSet (see rules.py) when a thing has several related checks, or when the
code must branch on the outcome.

Every problem found is collected and reported together, so a sheet can be fixed in one pass. Standard library only.
Nothing in this package imports TileSeqMut.
"""
import csv
import io
import itertools
import math
import re
from collections import Counter
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

from TileSeqDemux.rules import REPORTS, WARNS, Facts, Rule, RuleSet

ID_PATTERN = re.compile(r"^[^_\s/\\]+$")      # the sample ID is the first "_" token of a FASTQ name, so no "_"


class SheetLabels:
    """The exact first-column labels in the sheet. The Index Samples labels are provisional until the postdoc agrees."""
    SEQUENCE = "Sequence:"
    CDS_START = "CDS start:"
    CDS_END = "CDS end:"
    MIN_COVER = "Minimum coverage:"
    TILES = "Tile Number"                 # first cell of the tile table's header row
    SAMPLES = "Sample ID"                 # first cell of the sample table's header row
    SAMPLES_SECOND_COLUMN = "Tile ID"
    INDEX = "Index Samples"               # title row of the new section
    INDEX_HEADER = ("Shared Sample ID", "Unique Sample ID")
    DEFAULT_MIN_COVER = 0.6               # what tsm csv2json uses when the sheet has no Minimum coverage row
    ALLOWED_BASES = frozenset("ACGTN")


L = SheetLabels


class SheetError(Exception):
    """The sheet is not usable. `problems` holds every problem found."""

    def __init__(self, problems, path=None):
        self.problems = list(problems)
        self.path = path

    def __str__(self):
        count = len(self.problems)
        noun = "problem" if count == 1 else "problems"
        where = f" in {self.path}" if self.path else ""
        bullets = "\n".join(f"  - {problem}" for problem in self.problems)
        return f"{count} {noun}{where}:\n{bullets}"


@dataclass(frozen=True)
class Tile:
    number: str
    start_aa: int
    end_aa: int

    @property
    def begins(self):
        """First CDS position of the tile, as in tileseq_mut (count_mut.py)."""
        return self.start_aa * 3 - 2

    @property
    def ends(self):
        return self.end_aa * 3

    def min_map(self, min_cover):
        """Minimum mapped length, computed with the same float arithmetic as tileseq_mut."""
        return math.ceil((self.ends - self.begins) * min_cover)


@dataclass(frozen=True)
class Sample:
    sample_id: str
    tile: str
    condition: str = ""
    timepoint: str = ""
    replicate: str = ""


@dataclass(frozen=True)
class IndexGroup:
    """The samples that share one sequencing index: the shared sample first, then the unique ones."""
    shared: str
    unique: Tuple[str, ...]

    @property
    def members(self):
        return (self.shared,) + self.unique


@dataclass
class ParameterSheet:
    path: str
    sequence: str
    cds_start: int
    cds_end: Optional[int]
    min_cover: float
    tiles: Dict[str, Tile]
    samples: Dict[str, Sample]
    groups: List[IndexGroup]
    warnings: List[str] = field(default_factory=list)

    @property
    def shared_ids(self):
        """The sample IDs whose FASTQ files hold the pooled reads, one per index group."""
        return [group.shared for group in self.groups]

    @property
    def unique_ids(self):
        """The sample IDs that exist only after demultiplexing, from every index group."""
        return [unique for group in self.groups for unique in group.unique]

    def tile_of(self, sample_id):
        return self.tiles[self.samples[sample_id].tile]

    def describe(self, min_cover=None):
        """Human-readable lines: each index group, its tiles and their geometry."""
        cover = self.min_cover if min_cover is None else min_cover
        lines = [f"minimum coverage {cover}, CDS start {self.cds_start}"]
        for group in self.groups:
            lines += self._describe_group(group, cover)
        return lines

    def _describe_group(self, group, cover):
        lines = [f"index shared by {len(group.members)} samples:"]
        lines += [self._describe_member(group, sample_id, cover) for sample_id in group.members]
        return lines

    def _describe_member(self, group, sample_id, cover):
        tile = self.tile_of(sample_id)
        role = "  (shared: holds the FASTQ files)" if sample_id == group.shared else ""
        return (f"  sample {sample_id:<8} tile {tile.number:<3} CDS {tile.begins}-{tile.ends}, "
                f"minimum mapped length {tile.min_map(cover)}{role}")


def _rows(text):
    """(cells, line number) for every CSV record; cells are stripped."""
    reader = csv.reader(io.StringIO(text, newline=""))
    return [([c.strip() for c in row], reader.line_num) for row in reader]

# -- row rules ------------------------------------------------------------------------------------------------
BAD_ID_MESSAGE = ("sample ID {bad_id!r} must not be empty or contain '_', spaces or slashes "
                  "(the sample ID is the first '_' token of the FASTQ file name)")


def _whole_number(text):
    """The integer in `text`, or None when it is not a whole number."""
    try:
        return int(text)
    except ValueError:
        return None


@dataclass(frozen=True)
class TileRow:
    """One row of the tile table. The context its rules receive is the dict of tiles read so far."""
    number: str
    start_text: str
    end_text: str

    @property
    def start(self):
        return _whole_number(self.start_text)

    @property
    def end(self):
        return _whole_number(self.end_text)

    def details(self, known_tiles):
        return {"number": self.number, "start_text": self.start_text, "end_text": self.end_text, "start": self.start, "end": self.end}


TILE_RULES = RuleSet(
    details=("number", "start_text", "end_text", "start", "end"),
    rules=[
        Rule("not_numbers", "start and end amino acids must be whole numbers, found {start_text!r} and {end_text!r}",
             lambda row, known_tiles: row.start is None or row.end is None),
        Rule("bad_range", "tile {number} has start {start} and end {end}; need 1 <= start <= end",
             lambda row, known_tiles: row.start < 1 or row.end < row.start),
        Rule("duplicate", "tile {number} is listed twice",
             lambda row, known_tiles: row.number in known_tiles),
    ])


@dataclass(frozen=True)
class SampleContext:
    """What the sample rules need to know about the rest of the sheet."""
    known_samples: dict
    tiles: dict

    @property
    def tile_list(self):
        return ", ".join(sorted(self.tiles, key=lambda number: (len(number), number)))


@dataclass(frozen=True)
class SampleRow:
    """One row of the sample table."""
    sample_id: str
    tile: str
    condition: str
    timepoint: str
    replicate: str

    def details(self, context):
        return {"bad_id": self.sample_id, "sample_id": self.sample_id, "tile": self.tile, "known_tiles": context.tile_list}


SAMPLE_RULES = RuleSet(
    details=("bad_id", "sample_id", "tile", "known_tiles"),
    rules=[
        Rule("bad_id", BAD_ID_MESSAGE,
             lambda row, context: not ID_PATTERN.match(row.sample_id)),
        Rule("duplicate", "sample {sample_id} is listed twice",
             lambda row, context: row.sample_id in context.known_samples),
        Rule("unknown_tile", "sample {sample_id} uses tile {tile!r}, which is not in the tile table (tiles: {known_tiles})",
             lambda row, context: bool(context.tiles) and row.tile not in context.tiles,
             blocks_row=False),                       # reported, but the sample is kept
    ])


@dataclass(frozen=True)
class PairContext:
    """What the pair rules need to know about the rows accepted so far."""
    seen: set          # the (shared, unique) pairs
    owners: dict       # unique ID -> (shared ID, line number) of the row that used it first


@dataclass(frozen=True)
class PairRow:
    """One row of the Index Samples section."""
    shared: str
    unique: str

    @property
    def bad_id(self):
        """The first of the two IDs that cannot be a sample ID, or None."""
        return next((sample_id for sample_id in (self.shared, self.unique) if not ID_PATTERN.match(sample_id)), None)

    def details(self, context):
        first_shared, first_line = context.owners.get(self.unique, ("", ""))
        return {"shared": self.shared, "unique": self.unique, "bad_id": self.bad_id, "first_shared": first_shared, "first_line": first_line}


PAIR_RULES = RuleSet(
    details=("shared", "unique", "bad_id", "first_shared", "first_line"),
    rules=[
        Rule("no_unique", "shared sample {shared} has no unique sample ID",
             lambda row, context: row.unique == ""),
        Rule("bad_id", BAD_ID_MESSAGE,
             lambda row, context: row.bad_id is not None),
        Rule("same_id", "shared and unique sample ID are both {shared}",
             lambda row, context: row.shared == row.unique),
        Rule("pair_twice", "the pair {shared},{unique} is listed twice",
             lambda row, context: (row.shared, row.unique) in context.seen),
        Rule("unique_reused", "unique sample {unique} is already paired with shared sample {first_shared} on line {first_line}",
             lambda row, context: row.unique in context.owners),
    ])

def _decimal_number(text):
    """The float in `text`, or None when it is not a number."""
    try:
        return float(text)
    except ValueError:
        return None

# -- rules for everything else: single values, tables and sections --------------------------------------------
REPEATED_LABEL = Rule("repeated", "{label!r} appears {count} times (lines {lines})",
                      lambda facts, context: facts.count > 1)

LOOKUP_RULES = RuleSet(              # a labelled row that must appear exactly once
    details=("label", "what", "count", "lines", "hint"),
    rules=[
        Rule("missing", "no {what} labelled {label!r}{hint}", lambda facts, context: facts.count == 0),
        REPEATED_LABEL,
    ])

INDEX_LOOKUP_RULES = RuleSet(        # the title row of the Index Samples section
    details=("label", "count", "lines", "hint"),
    rules=[
        Rule("missing", "no {label!r} section: nothing to demultiplex{hint}", lambda facts, context: facts.count == 0),
        REPEATED_LABEL,
    ])

INTEGER_RULES = RuleSet(
    details=("label", "text", "number", "minimum"),
    rules=[
        Rule("not_whole", "{label!r} should be a whole number, found {text!r}",
             lambda facts, context: facts.number is None),
        Rule("too_small", "{label!r} should be at least {minimum}, found {number}",
             lambda facts, context: facts.number < facts.minimum),
    ])

COVER_RULES = RuleSet(
    details=("label", "count", "text", "number", "default"),
    rules=[
        Rule("repeated", "{label!r} appears {count} times",
             lambda facts, context: facts.count > 1),
        Rule("no_value", "no {label!r} value: using {default}, as tsm csv2json does",
             lambda facts, context: facts.text == "", effect=WARNS),
        Rule("not_a_number", "{label!r} should be a number, found {text!r}",
             lambda facts, context: facts.text != "" and facts.number is None),
        Rule("out_of_range", "{label!r} should be above 0 and at most 1, found {text}",
             lambda facts, context: facts.number is not None and not 0 < facts.number <= 1, effect=REPORTS),
    ])

INDEX_SECTION_RULES = RuleSet(       # the rows around the Index Samples title
    details=("label", "title_line", "above", "above_has_content", "has_header", "header_line", "found", "header_ok", "expected"),
    rules=[
        Rule("row_above_not_empty", "line {title_line}: the row above {label!r} must be a row of empty cells (',,,,'), not {above!r}; "
                                    "otherwise tsm csv2json reads the section as part of the table above",
             lambda facts, context: facts.above_has_content, effect=REPORTS),
        Rule("no_header_row", "{label!r} has no header row (expected {expected})",
             lambda facts, context: not facts.has_header),
        Rule("wrong_header", "line {header_line}: the header under {label!r} should be {expected}, found {found}",
             lambda facts, context: not facts.header_ok),
    ])

CDS_RULES = RuleSet(                 # the CDS against the sequence
    details=("start_label", "end_label", "start", "end", "sequence_length", "cds_length"),
    rules=[
        Rule("start_beyond_sequence", "{start_label!r} is {start} but the sequence has only {sequence_length} bases",
             lambda facts, context: facts.start > facts.sequence_length),
        Rule("end_before_start", "{end_label!r} ({end}) is before {start_label!r} ({start})",
             lambda facts, context: facts.end is not None and facts.end < facts.start),
        Rule("end_beyond_sequence", "{end_label!r} is {end} but the sequence has only {sequence_length} bases",
             lambda facts, context: facts.end is not None and facts.end > facts.sequence_length),
        Rule("length_not_multiple_of_3", "the CDS ({start} to {end}) is {cds_length} bases long, not a multiple of 3",
             lambda facts, context: facts.end is not None and facts.cds_length % 3 != 0, effect=WARNS),
    ])


class _SheetParse:
    def __init__(self, rows, path):
        self.rows = rows
        self.path = path
        self.problems = []
        self.warnings = []

    # -- helpers ---------------------------------------------------------------------------------------------
    def _check(self, rule_set, subject, context=None, where=None):
      """Apply `rule_set` to `subject` and record what it finds. The Verdict says whether the subject is usable."""
      verdict = rule_set.check(subject, context)
      prefix = f"{where}: " if where else ""
      self.problems += [prefix + message for message in verdict.problems]
      self.warnings += [prefix + message for message in verdict.warnings]
      return verdict

    def _record_problem_if(self, condition, message, where=None):
        """Record `message` as a problem when `condition` is true.

        This only records. It applies no rule and decides nothing. To act on a problem (stop, skip a row, use a default),
        use a RuleSet with _check, which returns a Verdict.
        """
        if condition:
            self.problems.append(f"{where}: {message}" if where else message)

    def _find(self, label):
        """Indexes of the rows whose first cell is exactly `label`."""
        return [index for index, (cells, _line_number) in enumerate(self.rows) if cells and cells[0] == label]

    def _hint(self, label):
        """A hint for a label that is present but written with different capitals, such as 'index samples'."""
        near_misses = [(cells[0], line_number) for cells, line_number in self.rows
                       if cells and cells[0].lower() == label.lower() and cells[0] != label]
        if not near_misses:
            return ""
        found, line_number = near_misses[0]
        return f" (found {found!r} on line {line_number}; labels must match exactly)"

    def _once(self, label, what):
        """The index of the one row labelled `label`, or None after recording that it is missing or repeated."""
        hits = self._find(label)
        lookup = Facts(label=label, what=what, count=len(hits), lines=self._lines(hits), hint=self._hint(label))
        return None if self._check(LOOKUP_RULES, lookup).blocked else hits[0]

    def _scalar(self, label, required=True):
        index = self._once(label, "row") if required else (self._find(label) or [None])[0]
        if index is None:
            return None
        _cells, line_number = self.rows[index]
        value = self._second_cell(index)
        self._record_problem_if(required and value == "", f"{label!r} on line {line_number} has no value")
        return value or None

    def _int(self, label, text, minimum=1):
        number = _whole_number(text)
        verdict = self._check(INTEGER_RULES, Facts(label=label, text=text, number=number, minimum=minimum))
        return None if verdict.blocked else number

    def _second_cell(self, row_index):
      cells, _line_number = self.rows[row_index]
      return cells[1] if len(cells) > 1 else ""

    @staticmethod
    def _padded(cells, width):
        """The first 'width' cells of a row, padded with empty strings if the row is shorter."""
        return (cells + [""] *width)[:width]

    @staticmethod
    def _has_first_cell(row):
        cells, _line_number = row
        return bool(cells) and cells[0] != ""

    def _table(self, header_index):
        """Rows after a header row, up to the first row with an empty first cell (where tsm csv2json ends a table)."""
        return list(itertools.takewhile(self._has_first_cell, self.rows[header_index + 1:]))

    # -- sections --------------------------------------------------------------------------------------------
    def sequence(self):
        value = self._scalar(L.SEQUENCE)
        if value is None:
            return None
        sequence = value.upper()
        bad = sorted(set(sequence) - L.ALLOWED_BASES)
        position = next((index for index, base in enumerate(sequence) if base in bad), 0) + 1
        shown = " ".join(repr(base) for base in bad[:5])
        self._record_problem_if(bool(bad), f"{L.SEQUENCE!r} contains characters other than A, C, G, T, N: {shown} (first at position {position})")
        return None if bad else sequence

    def min_cover(self):
        rows_found = self._find(L.MIN_COVER)
        text = self._second_cell(rows_found[0]) if rows_found else ""
        facts = Facts(label=L.MIN_COVER, count=len(rows_found), text=text, number=_decimal_number(text), default=L.DEFAULT_MIN_COVER)
        verdict = self._check(COVER_RULES, facts)
        return L.DEFAULT_MIN_COVER if verdict.blocked or facts.number is None else facts.number

    def tiles(self):
        header_index = self._once(L.TILES, "tile table header row")
        if header_index is None:
            return {}
        rows = self._table(header_index)
        self._record_problem_if(not rows, f"the tile table under {L.TILES!r} is empty")
        tiles = {}
        for cells, line_number in rows:
            tile = self._read_tile(cells, line_number, tiles)
            if tile is not None:
                tiles[tile.number] = tile
        return tiles

    def _read_tile(self, cells, line_number, known_tiles):
        """One Tile from one row of the tile table, or None after recording what is wrong with the row."""
        row = TileRow(*self._padded(cells, 3))
        verdict = self._check(TILE_RULES, row, known_tiles, where=f"tile table line {line_number}")
        return None if verdict.blocked else Tile(row.number, row.start, row.end)

    def samples(self, tiles):
        header_index = self._once(L.SAMPLES, "sample table header row")
        if header_index is None:
            return {}
        header_cells, _line_number = self.rows[header_index]
        found = self._padded(header_cells, 2)[1]
        self._record_problem_if(found != L.SAMPLES_SECOND_COLUMN, f"sample table: the second column should be {L.SAMPLES_SECOND_COLUMN!r}, found {found!r}")
        samples = {}
        for cells, line_number in self._table(header_index):
            sample = self._read_sample(cells, line_number, samples, tiles)
            if sample is not None:
                samples[sample.sample_id] = sample
        self._record_problem_if(not samples, f"the sample table under {L.SAMPLES!r} is empty")
        return samples

    def _read_sample(self, cells, line_number, known_samples, tiles):
      """One Sample from one row of the sample table, or None after recording what is wrong with the row."""
      row = SampleRow(*self._padded(cells, 5))
      verdict = self._check(SAMPLE_RULES, row, SampleContext(known_samples, tiles), where=f"sample table line {line_number}")
      return None if verdict.blocked else Sample(row.sample_id, row.tile, row.condition, row.timepoint, row.replicate)
    
    def groups(self, samples):
        rows = self._index_rows()
        if rows is None:
            return []
        self._record_problem_if(not rows, f"the {L.INDEX!r} section lists no samples")
        pairs = self._read_pairs(rows)
        self._check_samples_exist(pairs, samples)
        groups = self._build_groups(pairs)
        self._check_group_tiles(groups, samples)
        return groups

    def _index_rows(self):
        """The data rows of the Index Samples section, or None after recording why the section cannot be read."""
        titles = self._find(L.INDEX)
        lookup = Facts(label=L.INDEX, count=len(titles), lines=self._lines(titles), hint=self._hint(L.INDEX))
        if self._check(INDEX_LOOKUP_RULES, lookup).blocked:
            return None
        if self._check(INDEX_SECTION_RULES, self._index_layout(titles[0])).blocked:
            return None
        return self._table(titles[0] + 1)

    def _index_layout(self, title_index):
      """What the section rules need to know about the rows around the section title."""
      cells_above = self.rows[title_index - 1][0] if title_index > 0 else []
      header = self.rows[title_index + 1] if title_index + 1 < len(self.rows) else None
      header_cells, header_line = header if header else ([], None)
      return Facts(label=L.INDEX, title_line=self.rows[title_index][1],
                    above=",".join(cells_above), above_has_content=any(cells_above),
                    has_header=header is not None, header_line=header_line,
                    found=",".join(header_cells[:2]) or "an empty row", header_ok=tuple(header_cells[:2]) == L.INDEX_HEADER,
                    expected=",".join(L.INDEX_HEADER))

    def _read_pairs(self, rows):
        """(shared, unique, line_number) for every usable row of the section. Rows with a problem are reported and left out."""
        pairs, context = [], PairContext(seen=set(), owners={})
        for cells, line_number in rows:
            row = PairRow(*self._padded(cells, 2))
            verdict = self._check(PAIR_RULES, row, context, where=f"line {line_number}")
            if not verdict.blocked:
                pairs.append((row.shared, row.unique, line_number))
                context.seen.add((row.shared, row.unique))
                context.owners[row.unique] = (row.shared, line_number)
        return pairs

    @staticmethod
    def _mentions(pairs):
        """Every (sample_id, role, line_number) the section mentions, row by row."""
        for shared, unique, line_number in pairs:
            yield shared, "shared", line_number
            yield unique, "unique", line_number

    def _check_samples_exist(self, pairs, samples):
        """Every ID in the section must also be in the sample table, which gives it its tile."""
        for sample_id, role, line_number in self._mentions(pairs):
            self._record_problem_if(sample_id not in samples, f"{role} sample {sample_id} is not in the sample table", f"line {line_number}")

    def _build_groups(self, pairs):
        uniques_by_shared = {}
        for shared, unique, _line_number in pairs:
            uniques_by_shared.setdefault(shared, []).append(unique)
        for _shared, unique, line_number in pairs:
            self._record_problem_if(unique in uniques_by_shared, f"sample {unique} is a unique sample here and a shared sample elsewhere",
                             f"line {line_number}")
        return [IndexGroup(shared, tuple(uniques)) for shared, uniques in uniques_by_shared.items()]

    def _check_group_tiles(self, groups, samples):
        for group in groups:
            tiles_used = [samples[member].tile for member in group.members if member in samples]
            repeated = sorted(tile for tile, count in Counter(tiles_used).items() if count > 1)
            members, repeated_text = ", ".join(group.members), ", ".join(repeated)
            self._record_problem_if(bool(repeated), f"samples {members} share an index but tile {repeated_text} appears more than once; "
                                             "tiles in one index must differ")

    def run(self):
        sequence = self.sequence()
        start_value = self._scalar(L.CDS_START)
        cds_start = self._int(L.CDS_START, start_value, 1) if start_value is not None else None
        end_value = self._scalar(L.CDS_END, required=False)
        cds_end = self._int(L.CDS_END, end_value, 1) if end_value else None
        cover = self.min_cover()
        tiles = self.tiles()
        samples = self.samples(tiles)
        groups = self.groups(samples)
        self._cross_check(sequence, cds_start, cds_end, tiles)
        if self.problems:
            raise SheetError(self.problems, self.path)
        return ParameterSheet(self.path, sequence, cds_start, cds_end, cover, tiles, samples, groups, self.warnings)

    def _cross_check(self, sequence, cds_start, cds_end, tiles):
        """Checks that need several values at once: the CDS must fit in the sequence, and every tile must fit in the CDS."""
        if sequence is None or cds_start is None:
            return
        last_base = cds_end if cds_end is not None else len(sequence)
        cds = Facts(start_label=L.CDS_START, end_label=L.CDS_END, start=cds_start, end=cds_end,
                    sequence_length=len(sequence), cds_length=last_base - cds_start + 1)
        if self._check(CDS_RULES, cds).blocked:
            return
        for tile in tiles.values():
            self._record_problem_if(tile.ends > cds.cds_length,
                             f"tile {tile.number} ends at CDS position {tile.ends} (amino acid {tile.end_aa}) but the CDS has only {cds.cds_length} bases")


class SheetReader:
    """Turns a parameter sheet CSV into a ParameterSheet. It keeps no state, so one reader can read many sheets."""

    def read(self, path):
        """Read the file at `path`. Raises OSError if it cannot be read and SheetError if it is not usable."""
        with open(path, newline="", encoding="utf-8-sig") as fh:
            text = fh.read()
        return self.parse(text, path)

    def parse(self, text, path="<string>"):
        """Parse sheet text that is already in memory. Raises SheetError listing every problem."""
        return _SheetParse(_rows(text), path).run()