# These tests were generated with Claude (Anthropic) and then reviewed, edited and run by the maintainers.
"""
Tests for deciding which tile a read pair belongs to (slice C of ticket 1.4): the tile rule, every category, the
ambiguity policies, the counting, and agreement with the prototype of ticket 1.2 on random alignments.

Run with:  python -m pytest tests/test_demux_classify.py -v     (from the repository root)
Without pytest: python -m unittest discover -s tests -p "test_demux_*.py"
"""
import math
import random
import unittest

from TileSeqDemux import align as A
from TileSeqDemux import classify as C
from TileSeqDemux.sheet import SheetReader, Tile
from test_demux_sheet import make_sheet

TILE_1 = Tile("1", 2, 13)         # CDS 4 to 39, and 32 bases to cover at minimum coverage 0.9
TILE_12 = Tile("12", 373, 396)    # CDS 1117 to 1188, and 64 bases to cover at 0.9
INSIDE_1 = A.Alignment(position=4, length=41, span=41)
INSIDE_12 = A.Alignment(position=1117, length=72, span=72)


class TileRuleTests(unittest.TestCase):
    def rule(self, cds_start=1, min_cover=0.9, tile=TILE_1):
        return C.TileRule(tile, cds_start, min_cover)

    def test_the_numbers_are_the_tiles(self):
        self.assertEqual((TILE_1.begins, TILE_1.ends, TILE_1.min_map(0.9)), (4, 39, 32))
        self.assertEqual((TILE_12.begins, TILE_12.ends, TILE_12.min_map(0.9)), (1117, 1188, 64))

    def test_a_mate_inside_the_tile_is_accepted(self):
        self.assertTrue(self.rule().accepts(INSIDE_1))

    def test_a_mate_from_a_tile_further_downstream_passes_the_end_test_but_not_the_overlap_check(self):
        """tileseq_mut's test looks only at where a read ends, so this read passes it; the overlap check is the addition."""
        downstream = A.Alignment(position=100, length=50, span=50)
        self.assertTrue(downstream.position + downstream.length - 1 >= TILE_1.begins + 32)        # the end test alone says yes
        self.assertFalse(self.rule().accepts(downstream))

    def test_a_mate_that_ends_before_the_tile_is_far_enough_in_is_rejected(self):
        self.assertFalse(self.rule().accepts(A.Alignment(position=1, length=3, span=3)))

    def test_the_overlap_must_reach_the_minimum_exactly(self):
        """A mate starting inside the tile and running past its end overlaps it by (ends - start + 1) bases."""
        self.assertTrue(self.rule().accepts(A.Alignment(position=8, length=60, span=60)))      # 32 bases of overlap
        self.assertFalse(self.rule().accepts(A.Alignment(position=9, length=60, span=60)))     # 31

    def test_the_minimum_depends_on_the_minimum_coverage(self):
        mate = A.Alignment(position=15, length=30, span=30)
        self.assertTrue(self.rule(min_cover=0.6).accepts(mate))
        self.assertFalse(self.rule(min_cover=0.9).accepts(mate))

    def test_positions_are_relative_to_where_the_cds_starts(self):
        for cds_start in (1, 111, 500):
            with self.subTest(cds_start=cds_start):
                shifted = A.Alignment(INSIDE_1.position + cds_start - 1, INSIDE_1.length, INSIDE_1.span)
                self.assertTrue(self.rule(cds_start=cds_start).accepts(shifted))

    def test_the_overlap_uses_the_template_the_read_covers_not_the_length_of_the_read(self):
        """A long insertion makes a read longer than the template it covers. Its length still decides where it ends."""
        self.assertTrue(self.rule().accepts(A.Alignment(position=5, length=60, span=60)))        # covers CDS 5 to 64: 35 in the tile
        self.assertFalse(self.rule().accepts(A.Alignment(position=5, length=60, span=20)))       # covers CDS 5 to 24: 20 in the tile

    def test_the_rule_knows_its_tile(self):
        self.assertEqual(self.rule(tile=TILE_12).number, "12")


class PairClassifierTests(unittest.TestCase):
    def classifier(self, tiles=(TILE_1, TILE_12), ambiguous="drop"):
        return C.PairClassifier(tiles, cds_start=1, min_cover=0.9, ambiguous=ambiguous)

    def test_a_pair_inside_one_tile_is_assigned_to_it(self):
        self.assertEqual(self.classifier().classify(INSIDE_1, INSIDE_1), C.Assignment("assigned:1", ("1",)))
        self.assertEqual(self.classifier().classify(INSIDE_12, INSIDE_12), C.Assignment("assigned:12", ("12",)))

    def test_mates_in_different_tiles_are_off_tile(self):
        self.assertEqual(self.classifier().classify(INSIDE_1, INSIDE_12), C.Assignment("off_tile", ()))

    def test_a_pair_in_a_region_of_neither_tile_is_off_tile(self):
        elsewhere = A.Alignment(position=481, length=60, span=60)
        self.assertEqual(self.classifier().classify(elsewhere, elsewhere).category, "off_tile")

    def test_each_kind_of_unmapped_pair(self):
        classifier = self.classifier()
        self.assertEqual(classifier.classify(None, INSIDE_1), C.Assignment("r1_unmapped"))
        self.assertEqual(classifier.classify(INSIDE_1, None), C.Assignment("r2_unmapped"))
        self.assertEqual(classifier.classify(None, None), C.Assignment("both_unmapped"))

    def test_an_unmapped_pair_goes_to_no_tile(self):
        self.assertEqual(self.classifier().classify(None, None).tiles, ())

    def test_assigned_is_the_only_assigned_category(self):
        classifier = self.classifier()
        self.assertTrue(classifier.classify(INSIDE_1, INSIDE_1).is_assigned)
        for assignment in (classifier.classify(None, None), classifier.classify(INSIDE_1, INSIDE_12)):
            self.assertFalse(assignment.is_assigned)


class AmbiguityTests(unittest.TestCase):
    TWINS = (Tile("1", 2, 13), Tile("2", 2, 13))        # two tiles that both take the same reads

    def classify(self, ambiguous):
        return C.PairClassifier(self.TWINS, cds_start=1, min_cover=0.9, ambiguous=ambiguous).classify(INSIDE_1, INSIDE_1)

    def test_a_pair_that_two_tiles_take_is_ambiguous(self):
        self.assertEqual(self.classify("drop").category, "ambiguous")

    def test_by_default_an_ambiguous_pair_is_written_nowhere(self):
        self.assertEqual(self.classify("drop").tiles, ())

    def test_with_both_an_ambiguous_pair_goes_to_every_tile_that_takes_it_in_the_order_given(self):
        self.assertEqual(self.classify("both").tiles, ("1", "2"))

    def test_an_unknown_policy_is_refused_with_the_choices(self):
        with self.assertRaisesRegex(ValueError, r"\['both', 'drop'\].*'sometimes'"):
            C.PairClassifier(self.TWINS, 1, 0.9, ambiguous="sometimes")

    def test_tiles_can_be_given_as_any_iterable(self):
        classifier = C.PairClassifier((tile for tile in (TILE_1, TILE_12)), 1, 0.9)
        self.assertEqual(classifier.classify(INSIDE_12, INSIDE_12).category, "assigned:12")


class TallyTests(unittest.TestCase):
    def fill(self, tally, *assignments):
        for assignment in assignments:
            tally.add(assignment)
        return tally

    def test_pairs_are_counted_by_category_and_by_tile(self):
        tally = self.fill(C.Tally(["1", "12"]), C.Assignment("assigned:1", ("1",)), C.Assignment("assigned:1", ("1",)),
                          C.Assignment("assigned:12", ("12",)), C.Assignment("off_tile"), C.Assignment("both_unmapped"))
        self.assertEqual(tally.categories, {"assigned:1": 2, "assigned:12": 1, "off_tile": 1, "both_unmapped": 1})
        self.assertEqual(tally.written, {"1": 2, "12": 1})

    def test_assigned_and_unassigned_add_up_to_the_pairs_in(self):
        tally = self.fill(C.Tally(["1", "12"]), C.Assignment("assigned:1", ("1",)), C.Assignment("off_tile"),
                          C.Assignment("ambiguous"), C.Assignment("r1_unmapped"), C.Assignment("assigned:12", ("12",)))
        self.assertEqual((tally.pairs_in, tally.assigned, tally.unassigned), (5, 2, 3))
        self.assertEqual(tally.assigned + tally.unassigned, tally.pairs_in)

    def test_an_ambiguous_pair_written_to_both_counts_once_as_a_pair_and_once_for_each_tile(self):
        tally = self.fill(C.Tally(["1", "2"]), C.Assignment("ambiguous", ("1", "2")))
        self.assertEqual((tally.pairs_in, tally.unassigned, tally.written), (1, 1, {"1": 1, "2": 1}))

    def test_a_tile_nobody_was_written_to_is_still_listed_with_zero(self):
        self.assertEqual(C.Tally(["1", "12"]).written, {"1": 0, "12": 0})

    def test_a_tile_that_is_not_the_jobs_is_an_error(self):
        with self.assertRaises(KeyError):
            C.Tally(["1"]).add(C.Assignment("assigned:7", ("7",)))

    def test_the_summary_is_plain_numbers_in_a_stable_order(self):
        tally = self.fill(C.Tally(["1"]), C.Assignment("off_tile"), C.Assignment("assigned:1", ("1",)))
        self.assertEqual(tally.summary(), {"pairs_in": 2, "assigned": 1, "unassigned": 1,
                                           "categories": {"assigned:1": 1, "off_tile": 1}, "written": {"1": 1}})


class PassOverSamTests(unittest.TestCase):
    """The pieces together on hand-made SAM lines: the join of two streams to a list of reads, as the job will do it."""

    @staticmethod
    def line(name, position, length):
        fields = [name, "0", "gene", str(position), "42", f"{length}M", "*", "0", "0", "A" * length, "I" * length, "AS:i:1"]
        return ("\t".join(fields) + "\n").encode()

    def test_a_pass_gives_every_category_and_uses_up_both_streams(self):
        ids = [b"p1", b"p2", b"p3", b"p4", b"p5", b"p6", b"p7"]
        r1 = [self.line("p1", 4, 41), self.line("p2", 1117, 72), self.line("p4", 481, 60), self.line("p5", 4, 41), self.line("p7", 4, 41)]
        r2 = [self.line("p1", 4, 41), self.line("p2", 1117, 72), self.line("p3", 4, 41), self.line("p4", 481, 60), self.line("p5", 1117, 72)]
        first, second = A.AlignmentStream(r1), A.AlignmentStream(r2)
        classifier = C.PairClassifier([TILE_1, TILE_12], 1, 0.9)
        tally = C.Tally(["1", "12"])
        for read_id in ids:
            tally.add(classifier.classify(first.take(read_id), second.take(read_id)))
        self.assertEqual(dict(tally.categories), {"assigned:1": 1, "assigned:12": 1, "r1_unmapped": 1, "r2_unmapped": 1,
                                                  "both_unmapped": 1, "off_tile": 2})
        self.assertEqual((tally.pairs_in, tally.assigned, tally.unassigned), (7, 2, 5))
        self.assertEqual((first.leftover, second.leftover), (None, None))
        self.assertEqual((first.records_read, second.records_read), (5, 5))


# ---------------------------------------------------------------------------------------------------------------
# The prototype of ticket 1.2 (prototype_gene.py), copied here exactly as it was. The port must agree with it.
def reference_geometry(start_aa, end_aa, min_cover):
    begins, ends = start_aa * 3 - 2, end_aa * 3
    return {"begins": begins, "ends": ends, "min_map": math.ceil((ends - begins) * min_cover)}


def reference_mate_ok(aln, cds_start, tile):
    pos, seq_len, span = aln
    reaches = pos + seq_len - cds_start >= tile["begins"] + tile["min_map"]
    first = pos - cds_start + 1
    overlap = min(first + span - 1, tile["ends"]) - max(first, tile["begins"]) + 1
    return reaches and overlap >= tile["min_map"]


REFERENCE_UNMAPPED = {(False, False): None, (True, False): "r1_unmapped", (False, True): "r2_unmapped", (True, True): "both_unmapped"}
REFERENCE_AMBIGUOUS = {"drop": lambda hits: [], "both": lambda hits: hits}


def reference_classify(a1, a2, geom, cds_start, ambiguous):
    kind = REFERENCE_UNMAPPED[(a1 is None, a2 is None)]
    if kind:
        return kind, []
    hits = [t for t, g in geom.items() if reference_mate_ok(a1, cds_start, g) and reference_mate_ok(a2, cds_start, g)]
    if len(hits) == 1:
        return "assigned:" + hits[0], hits
    if not hits:
        return "off_tile", []
    return "ambiguous", REFERENCE_AMBIGUOUS[ambiguous](hits)
# ---------------------------------------------------------------------------------------------------------------


class AgreementWithThePrototypeTests(unittest.TestCase):
    """The port must give the prototype's answer for every pair. Pairs are random but placed near the tiles, so every
    category turns up, and the test checks that they did."""

    @classmethod
    def setUpClass(cls):
        cls.sheet = SheetReader().parse(make_sheet())
        cls.cds_start = cls.sheet.cds_start

    def random_mate(self, rng, tiles):
        if rng.random() < 0.05:
            return None
        tile = rng.choice(tiles)
        cds_position = tile.begins + rng.randint(-80, 60)                # near the tile, before it, inside it or just past it
        length = rng.randint(30, 150)
        return A.Alignment(self.cds_start + cds_position - 1, length, max(10, length + rng.randint(-8, 8)))

    def compare(self, tiles, min_cover, ambiguous, pairs=6000, seed=0):
        rng = random.Random(seed)
        reference_geom = {tile.number: reference_geometry(tile.start_aa, tile.end_aa, min_cover) for tile in tiles}
        classifier = C.PairClassifier(tiles, self.cds_start, min_cover, ambiguous)
        categories = set()
        for _ in range(pairs):
            first, second = self.random_mate(rng, tiles), self.random_mate(rng, tiles)
            # most pairs have both mates in the same region, as real pairs do
            if first is not None and second is not None and rng.random() < 0.7:
                second = A.Alignment(first.position + rng.randint(-30, 30), second.length, second.span)
            expected = reference_classify(first, second, reference_geom, self.cds_start, ambiguous)
            got = classifier.classify(first, second)
            self.assertEqual((got.category, list(got.tiles)), expected, f"pair {first} {second}")
            categories.add(got.category.split(":")[0])
        return categories

    def test_agreement_on_real_tiles_far_apart(self):
        tiles = [self.sheet.tiles["1"], self.sheet.tiles["12"]]
        for min_cover in (0.6, 0.9, 1.0):
            self.compare(tiles, min_cover, "drop")

    def test_agreement_on_neighbouring_tiles(self):
        tiles = [self.sheet.tiles["6"], self.sheet.tiles["7"]]
        for min_cover in (0.6, 0.9, 1.0):
            self.compare(tiles, min_cover, "drop")

    def test_agreement_on_overlapping_tiles_with_both_policies(self):
        tiles = [Tile("a", 2, 13), Tile("b", 6, 20), Tile("c", 3, 18)]
        for ambiguous in ("drop", "both"):
            for min_cover in (0.5, 0.9):
                categories = self.compare(tiles, min_cover, ambiguous, seed=ord(ambiguous[0]))
                self.assertIn("ambiguous", categories)

    def test_every_category_turns_up_so_the_agreement_is_not_vacuous(self):
        tiles = [Tile("a", 2, 13), Tile("b", 6, 20)]
        seen = set()
        for seed in range(3):
            seen |= self.compare(tiles, 0.6, "drop", seed=seed)
        self.assertEqual(seen, {"assigned", "off_tile", "ambiguous", "r1_unmapped", "r2_unmapped", "both_unmapped"})


if __name__ == "__main__":
    unittest.main()