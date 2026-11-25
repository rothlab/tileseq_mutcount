"""Deciding which tile a read pair belongs to (slice C of ticket 1.4).

    TileRule         what one tile asks of a mate: tileseq_mut's test, plus the overlap check pooled tiles need
    PairClassifier   puts a pair into one category, using one TileRule per tile of the job
    Assignment       the answer: a category, and the tiles whose files the pair is written to
    Tally            counts the pairs of a job by category and by tile

tileseq_mut keeps a read if it ends far enough into a tile. That test looks only at where a read ends, so a read from a
tile further downstream also passes it. Two tiles that share one index make that a real problem, so a mate must also
overlap the tile by at least the minimum mapped length. The arithmetic is tileseq_mut's (count_mut.py): a tile covers
CDS positions begins to ends, and a mate has to cover min_map of them, where min_map is the minimum coverage times the
tile's length, rounded up.

Categories: assigned:<tile>, r1_unmapped, r2_unmapped, both_unmapped, off_tile (no tile takes both mates) and
ambiguous (more than one does; dropped by default, or written to every tile that takes it).

Standard library only; tiles and alignments are used by their attributes, so nothing here imports another module.
"""
from collections import Counter
from typing import NamedTuple, Tuple

ASSIGNED = "assigned:"
R1_UNMAPPED = "r1_unmapped"
R2_UNMAPPED = "r2_unmapped"
BOTH_UNMAPPED = "both_unmapped"
OFF_TILE = "off_tile"
AMBIGUOUS = "ambiguous"

POLICIES = {"drop": lambda hits: (), "both": tuple}      # what to write when more than one tile takes a pair


class Assignment(NamedTuple):
    """Immutable, so the answers that come up for most pairs are made once and shared."""
    category: str
    tiles: Tuple[str, ...] = ()     # the tile numbers whose files the pair is written to

    @property
    def is_assigned(self):
        return self.category.startswith(ASSIGNED)


class UnmappedMates(NamedTuple):
    """Which mates of a pair have no alignment."""
    r1: bool
    r2: bool


_UNMAPPED = {
    UnmappedMates(r1=False, r2=False): None,                       # both mates aligned, so the pair is classified further
    UnmappedMates(r1=True, r2=False): Assignment(R1_UNMAPPED),
    UnmappedMates(r1=False, r2=True): Assignment(R2_UNMAPPED),
    UnmappedMates(r1=True, r2=True): Assignment(BOTH_UNMAPPED),
}
_OFF_TILE = Assignment(OFF_TILE)


class TileRule:
    """What one tile asks of a mate. The tile's numbers are worked out once, because this runs for every mate."""

    def __init__(self, tile, cds_start, min_cover):
        self.number = tile.number
        self._begins = tile.begins
        self._ends = tile.ends
        self._minimum = tile.min_map(min_cover)
        self._cds_start = cds_start

    def accepts(self, alignment):
        """True if the mate ends far enough into the tile (tileseq_mut's test) and overlaps it by the minimum."""
        reaches = alignment.position + alignment.length - self._cds_start >= self._begins + self._minimum
        first = alignment.position - self._cds_start + 1          # the CDS position of the first aligned base
        overlap = min(first + alignment.span - 1, self._ends) - max(first, self._begins) + 1
        return reaches and overlap >= self._minimum


class PairClassifier:
    """Puts read pairs into categories. `tiles` are the tiles of one job, each with number, begins, ends and
    min_map(cover); `ambiguous` is "drop" or "both"."""

    def __init__(self, tiles, cds_start, min_cover, ambiguous="drop"):
        if ambiguous not in POLICIES:
            raise ValueError(f"ambiguous must be one of {sorted(POLICIES)}, not {ambiguous!r}")
        self._rules = [TileRule(tile, cds_start, min_cover) for tile in tiles]
        self._assignments = {rule.number: Assignment(ASSIGNED + rule.number, (rule.number,)) for rule in self._rules}
        self._policy = POLICIES[ambiguous]
        self._outcomes = (self._off_tile, self._assigned, self._ambiguous)       # for no tile, one tile, several tiles

    def classify(self, first, second):
        """The Assignment for a pair, given the Alignment of each mate or None where a mate did not align."""
        unmapped = _UNMAPPED[UnmappedMates(r1=first is None, r2=second is None)]
        if unmapped:
            return unmapped
        hits = [rule.number for rule in self._rules if rule.accepts(first) and rule.accepts(second)]
        return self._outcomes[min(len(hits), 2)](hits)

    @staticmethod
    def _off_tile(hits):
        return _OFF_TILE

    def _assigned(self, hits):
        return self._assignments[hits[0]]

    def _ambiguous(self, hits):
        return Assignment(AMBIGUOUS, tuple(self._policy(hits)))


class Tally:
    """The pairs of one job counted by category and by tile. `tiles` are the tile numbers of the job."""

    def __init__(self, tiles):
        self.categories = Counter()
        self.written = {tile: 0 for tile in tiles}       # pairs written per tile; an ambiguous pair can count for two
        self.pairs_in = 0

    def add(self, assignment):
        self.pairs_in += 1
        self.categories[assignment.category] += 1
        for tile in assignment.tiles:
            self.written[tile] += 1

    @property
    def assigned(self):
        return sum(count for category, count in self.categories.items() if category.startswith(ASSIGNED))

    @property
    def unassigned(self):
        return sum(count for category, count in self.categories.items() if not category.startswith(ASSIGNED))

    def summary(self):
        """Everything the run record needs, as plain numbers."""
        return {"pairs_in": self.pairs_in, "assigned": self.assigned, "unassigned": self.unassigned,
                "categories": dict(sorted(self.categories.items())), "written": dict(self.written)}