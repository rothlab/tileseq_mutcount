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

Every problem found is collected and reported together, so a sheet can be fixed in one pass. Standard library only.
Nothing in this package imports TileSeqMut.
"""
import csv
import io
import math
import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

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


class _SheetParse:
    def __init__(self, rows, path):
        self.rows = rows
        self.path = path
        self.problems = []
        self.warnings = []

    # -- helpers ---------------------------------------------------------------------------------------------
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
        hits = self._find(label)
        if not hits:
            self.problems.append("no %s labelled %r%s" % (what, label, self._hint(label)))
            return None
        if len(hits) > 1:
            self.problems.append("%r appears %d times (lines %s)" % (label, len(hits), ", ".join(str(self.rows[i][1]) for i in hits)))
            return None
        return hits[0]

    def _scalar(self, label, required=True):
        i = self._once(label, "row") if required else (self._find(label) or [None])[0]
        if i is None:
            return None
        row, line = self.rows[i]
        value = row[1] if len(row) > 1 else ""
        if value == "":
            if required:
                self.problems.append("%r on line %d has no value" % (label, line))
            return None
        return value

    def _int(self, label, value, minimum=1):
        try:
            number = int(value)
        except ValueError:
            self.problems.append("%r should be a whole number, found %r" % (label, value))
            return None
        if number < minimum:
            self.problems.append("%r should be at least %d, found %d" % (label, minimum, number))
            return None
        return number

    def _table(self, header_index):
        """Rows after a header row, up to the first row with an empty first cell (where tsm csv2json ends a table)."""
        out = []
        for r, line in self.rows[header_index + 1:]:
            if not r or r[0] == "":
                break
            out.append((r, line))
        return out

    def _bad_id(self, value, where):
        if not ID_PATTERN.match(value):
            self.problems.append("%s: sample ID %r must not be empty or contain '_', spaces or slashes "
                                 "(the sample ID is the first '_' token of the FASTQ file name)" % (where, value))
            return True
        return False

    # -- sections --------------------------------------------------------------------------------------------
    def sequence(self):
        value = self._scalar(L.SEQUENCE)
        if value is None:
            return None
        seq = value.upper()
        bad = sorted(set(seq) - L.ALLOWED_BASES)
        if bad:
            first = next(i for i, ch in enumerate(seq) if ch in bad)
            self.problems.append("%r contains characters other than A, C, G, T, N: %s (first at position %d)"
                                 % (L.SEQUENCE, " ".join(repr(b) for b in bad[:5]), first + 1))
            return None
        return seq

    def min_cover(self):
        hits = self._find(L.MIN_COVER)
        if len(hits) > 1:
            self.problems.append("%r appears %d times" % (L.MIN_COVER, len(hits)))
            return L.DEFAULT_MIN_COVER
        value = self.rows[hits[0]][0][1] if hits and len(self.rows[hits[0]][0]) > 1 else ""
        if not hits or value == "":
            self.warnings.append("no %r value: using %s, as tsm csv2json does" % (L.MIN_COVER, L.DEFAULT_MIN_COVER))
            return L.DEFAULT_MIN_COVER
        try:
            cover = float(value)
        except ValueError:
            self.problems.append("%r should be a number, found %r" % (L.MIN_COVER, value))
            return L.DEFAULT_MIN_COVER
        if not 0 < cover <= 1:
            self.problems.append("%r should be above 0 and at most 1, found %s" % (L.MIN_COVER, value))
        return cover

    def tiles(self):
        i = self._once(L.TILES, "tile table header row")
        tiles = {}
        if i is None:
            return tiles
        rows = self._table(i)
        if not rows:
            self.problems.append("the tile table under %r is empty" % L.TILES)
        for r, line in rows:
            r = r + [""] * (3 - len(r))
            start, end = r[1], r[2]
            try:
                s_aa, e_aa = int(start), int(end)
            except ValueError:
                self.problems.append("tile table line %d: start and end amino acids must be whole numbers, found %r and %r" % (line, start, end))
                continue
            if s_aa < 1 or e_aa < s_aa:
                self.problems.append("tile table line %d: tile %s has start %d and end %d; need 1 <= start <= end" % (line, r[0], s_aa, e_aa))
                continue
            if r[0] in tiles:
                self.problems.append("tile table line %d: tile %s is listed twice" % (line, r[0]))
                continue
            tiles[r[0]] = Tile(r[0], s_aa, e_aa)
        return tiles

    def samples(self, tiles):
        i = self._once(L.SAMPLES, "sample table header row")
        samples = {}
        if i is None:
            return samples
        header = self.rows[i][0] + [""] * 2
        if header[1] != L.SAMPLES_SECOND_COLUMN:
            self.problems.append("sample table: the second column should be %r, found %r" % (L.SAMPLES_SECOND_COLUMN, header[1]))
        for r, line in self._table(i):
            r = r + [""] * (5 - len(r))
            sid, tile = r[0], r[1]
            if self._bad_id(sid, "sample table line %d" % line):
                continue
            if sid in samples:
                self.problems.append("sample table line %d: sample %s is listed twice" % (line, sid))
                continue
            if tiles and tile not in tiles:
                self.problems.append("sample table line %d: sample %s uses tile %r, which is not in the tile table (tiles: %s)"
                                     % (line, sid, tile, ", ".join(sorted(tiles, key=lambda t: (len(t), t)))))
            samples[sid] = Sample(sid, tile, r[2], r[3], r[4])
        if not samples:
            self.problems.append("the sample table under %r is empty" % L.SAMPLES)
        return samples

    def groups(self, samples):
        titles = self._find(L.INDEX)
        if not titles:
            self.problems.append("no %r section: nothing to demultiplex%s" % (L.INDEX, self._hint(L.INDEX)))
            return []
        if len(titles) > 1:
            self.problems.append("%r appears %d times (lines %s)" % (L.INDEX, len(titles), ", ".join(str(self.rows[i][1]) for i in titles)))
            return []
        t = titles[0]
        if t > 0 and any(self.rows[t - 1][0]):
            self.problems.append("line %d: the row above %r must be a row of empty cells (',,,,'), not %r; otherwise tsm csv2json "
                                 "reads the section as part of the table above" % (self.rows[t][1], L.INDEX, ",".join(self.rows[t - 1][0])))
        if t + 1 >= len(self.rows):
            self.problems.append("%r has no header row (expected %s)" % (L.INDEX, ",".join(L.INDEX_HEADER)))
            return []
        header, hline = self.rows[t + 1]
        if tuple(header[:2]) != L.INDEX_HEADER:
            self.problems.append("line %d: the header under %r should be %s, found %s"
                                 % (hline, L.INDEX, ",".join(L.INDEX_HEADER), ",".join(header[:2]) or "an empty row"))
            return []
        by_shared, seen_pairs, uniques = {}, set(), {}
        for r, line in self._table(t + 1):
            r = r + [""] * (2 - len(r))
            shared, unique = r[0], r[1]
            if unique == "":
                self.problems.append("line %d: shared sample %s has no unique sample ID" % (line, shared))
                continue
            bad = self._bad_id(shared, "line %d" % line) | self._bad_id(unique, "line %d" % line)
            if bad:
                continue
            if shared == unique:
                self.problems.append("line %d: shared and unique sample ID are both %s" % (line, shared))
                continue
            if (shared, unique) in seen_pairs:
                self.problems.append("line %d: the pair %s,%s is listed twice" % (line, shared, unique))
                continue
            seen_pairs.add((shared, unique))
            for sid, role in ((shared, "shared"), (unique, "unique")):
                if sid not in samples:
                    self.problems.append("line %d: %s sample %s is not in the sample table" % (line, role, sid))
            if unique in uniques:
                self.problems.append("line %d: unique sample %s is already paired with shared sample %s on line %d"
                                     % (line, unique, uniques[unique][0], uniques[unique][1]))
                continue
            uniques[unique] = (shared, line)
            by_shared.setdefault(shared, []).append(unique)
        for unique, (shared, line) in uniques.items():
            if unique in by_shared:
                self.problems.append("line %d: sample %s is a unique sample here and a shared sample elsewhere" % (line, unique))
        groups = [IndexGroup(s, tuple(u)) for s, u in by_shared.items()]
        for g in groups:
            tiles_in = [samples[m].tile for m in g.members if m in samples]
            dupes = sorted({x for x in tiles_in if tiles_in.count(x) > 1})
            if dupes:
                self.problems.append("samples %s share an index but tile %s appears more than once; tiles in one index must differ"
                                     % (", ".join(g.members), ", ".join(dupes)))
        if not groups and not any("section" in p for p in self.problems):
            self.problems.append("the %r section lists no samples" % L.INDEX)
        return groups

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
        if sequence is None or cds_start is None:
            return
        if cds_start > len(sequence):
            self.problems.append("%r is %d but the sequence has only %d bases" % (L.CDS_START, cds_start, len(sequence)))
            return
        if cds_end is not None:
            if cds_end < cds_start:
                self.problems.append("%r (%d) is before %r (%d)" % (L.CDS_END, cds_end, L.CDS_START, cds_start))
                return
            if cds_end > len(sequence):
                self.problems.append("%r is %d but the sequence has only %d bases" % (L.CDS_END, cds_end, len(sequence)))
                return
            if (cds_end - cds_start + 1) % 3:
                self.warnings.append("the CDS (%d to %d) is %d bases long, not a multiple of 3" % (cds_start, cds_end, cds_end - cds_start + 1))
        cds_len = (cds_end - cds_start + 1) if cds_end is not None else len(sequence) - cds_start + 1
        for t in tiles.values():
            if t.ends > cds_len:
                self.problems.append("tile %s ends at CDS position %d (amino acid %d) but the CDS has only %d bases"
                                     % (t.number, t.ends, t.end_aa, cds_len))


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