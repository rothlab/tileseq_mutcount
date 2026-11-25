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
                                           "categories": {"assigned:1": 1, "off_tile": 1}, "written": {"1": 1},
                                           "off_tile_regions": {}})

    def test_off_tile_pairs_are_counted_by_region_and_only_those(self):
        tally = self.fill(C.Tally(["1", "12"]), C.Assignment("off_tile", (), "6"), C.Assignment("off_tile", (), "6"),
                          C.Assignment("off_tile", (), "mixed"), C.Assignment("off_tile", (), "12"),
                          C.Assignment("assigned:1", ("1",)), C.Assignment("r1_unmapped"))
        self.assertEqual(tally.off_tile_regions, {"6": 2, "mixed": 1, "12": 1})
        self.assertEqual(sum(tally.off_tile_regions.values()), tally.categories["off_tile"])

    def test_a_region_adds_nothing_to_the_tile_counts(self):
        tally = self.fill(C.Tally(["1", "12"]), C.Assignment("off_tile", (), "6"))
        self.assertEqual((tally.written, tally.assigned, tally.unassigned), ({"1": 0, "12": 0}, 0, 1))

    def test_the_summary_lists_tiles_in_numeric_order_then_the_other_labels(self):
        tally = self.fill(C.Tally([]), *[C.Assignment("off_tile", (), label) for label in ("outside_tiles", "10", "mixed", "2", "1")])
        self.assertEqual(list(tally.summary()["off_tile_regions"]), ["1", "2", "10", "mixed", "outside_tiles"])


ALL_TILES = [Tile(str(number), start, end) for number, (start, end) in enumerate(
    [(2, 13), (14, 49), (50, 87), (88, 123), (124, 160), (161, 195), (196, 229), (230, 264), (265, 301), (302, 337), (338, 372), (373, 396)], 1)]


class RegionMapTests(unittest.TestCase):
    """Where on the gene a mate lies, using the twelve tiles of the sheet (tile 3 is CDS 148 to 261)."""

    def regions(self, cds_start=1, tiles=ALL_TILES):
        return C.RegionMap(tiles, cds_start)

    def test_the_tile_table_is_the_one_in_the_test_sheet(self):
        sheet_tiles = SheetReader().parse(make_sheet()).tiles
        self.assertEqual([(tile.number, tile.begins, tile.ends) for tile in ALL_TILES],
                         [(number, sheet_tiles[number].begins, sheet_tiles[number].ends) for number in sorted(sheet_tiles, key=int)])

    def test_a_mate_inside_a_tile_is_in_that_tile(self):
        self.assertEqual(self.regions().of_mate(A.Alignment(200, 60, 60)), "3")

    def test_a_mate_over_two_tiles_is_in_the_one_it_overlaps_more(self):
        self.assertEqual(self.regions().of_mate(A.Alignment(130, 60, 60)), "3")        # CDS 130 to 189: 18 bases in tile 2, 42 in tile 3

    def test_equal_overlap_goes_to_the_lower_tile_whatever_the_order_the_tiles_were_given_in(self):
        mate = A.Alignment(118, 60, 60)                                               # CDS 118 to 177: 30 bases in tile 2, 30 in tile 3
        self.assertEqual(self.regions().of_mate(mate), "2")
        self.assertEqual(self.regions(tiles=list(reversed(ALL_TILES))).of_mate(mate), "2")

    def test_a_mate_beyond_either_end_is_outside_the_tiles(self):
        self.assertEqual(self.regions().of_mate(A.Alignment(1300, 80, 80)), "outside_tiles")
        self.assertEqual(self.regions().of_mate(A.Alignment(-50, 40, 40)), "outside_tiles")

    def test_a_mate_that_ends_one_base_before_the_first_tile_or_starts_one_base_after_the_last_is_outside(self):
        self.assertEqual(self.regions().of_mate(A.Alignment(-36, 40, 40)), "outside_tiles")        # CDS -36 to 3; tile 1 starts at 4
        self.assertEqual(self.regions().of_mate(A.Alignment(1189, 40, 40)), "outside_tiles")       # tile 12 ends at 1188

    def test_a_mate_that_only_just_reaches_a_tile_is_in_it(self):
        self.assertEqual(self.regions().of_mate(A.Alignment(1, 4, 4)), "1")           # CDS 1 to 4: one base in tile 1

    def test_the_cds_start_is_taken_into_account(self):
        self.assertEqual(self.regions(cds_start=111).of_mate(A.Alignment(310, 60, 60)), "3")        # CDS 200 to 259

    def test_the_span_and_not_the_length_decides(self):
        self.assertEqual(self.regions().of_mate(A.Alignment(200, 150, 20)), "3")      # a long insertion: only 20 bases of template

    def test_two_mates_in_one_place_give_that_region(self):
        self.assertEqual(self.regions().of_pair(A.Alignment(200, 60, 60), A.Alignment(210, 50, 50)), "3")

    def test_two_mates_in_different_places_are_mixed(self):
        self.assertEqual(self.regions().of_pair(A.Alignment(200, 60, 60), A.Alignment(500, 60, 60)), "mixed")

    def test_a_mate_outside_and_a_mate_inside_are_mixed(self):
        self.assertEqual(self.regions().of_pair(A.Alignment(200, 60, 60), A.Alignment(1300, 60, 60)), "mixed")

    def test_two_mates_outside_are_outside_tiles(self):
        self.assertEqual(self.regions().of_pair(A.Alignment(1300, 60, 60), A.Alignment(1320, 60, 60)), "outside_tiles")

    def test_with_no_tiles_everything_is_outside(self):
        self.assertEqual(self.regions(tiles=[]).of_mate(A.Alignment(200, 60, 60)), "outside_tiles")


class OffTileRegionTests(unittest.TestCase):
    """The classifier with the sheet's tiles given for labelling. The job's tiles are 1 and 12; positions are CDS positions."""

    def setUp(self):
        self.sheet = SheetReader().parse(make_sheet())
        self.classifier = C.PairClassifier([TILE_1, TILE_12], 1, 0.9, regions=list(self.sheet.tiles.values()))

    def test_an_off_tile_pair_is_labelled_with_where_its_mates_lie(self):
        mate = A.Alignment(481, 60, 60)                                               # inside tile 6
        self.assertEqual(self.classifier.classify(mate, mate), C.Assignment("off_tile", (), "6"))

    def test_mates_in_different_tiles_are_labelled_mixed(self):
        self.assertEqual(self.classifier.classify(INSIDE_1, INSIDE_12), C.Assignment("off_tile", (), "mixed"))

    def test_a_pair_too_short_inside_this_jobs_own_tile_is_labelled_with_that_tile(self):
        short = A.Alignment(1125, 45, 45)                                             # 45 of the 64 bases tile 12 needs
        self.assertEqual(self.classifier.classify(short, short), C.Assignment("off_tile", (), "12"))

    def test_a_pair_beyond_the_gene_is_labelled_outside_tiles(self):
        beyond = A.Alignment(1300, 60, 60)
        self.assertEqual(self.classifier.classify(beyond, beyond), C.Assignment("off_tile", (), "outside_tiles"))

    def test_an_off_tile_pair_is_still_written_to_no_tile(self):
        mate = A.Alignment(481, 60, 60)
        self.assertEqual(self.classifier.classify(mate, mate).tiles, ())

    def test_only_off_tile_pairs_get_a_region(self):
        inside = self.classifier.classify(INSIDE_1, INSIDE_1)
        unmapped = self.classifier.classify(None, INSIDE_1)
        twins = C.PairClassifier([Tile("1", 2, 13), Tile("2", 2, 13)], 1, 0.9, regions=ALL_TILES).classify(INSIDE_1, INSIDE_1)
        self.assertEqual((inside.region, unmapped.region, twins.category, twins.region), ("", "", "ambiguous", ""))

    def test_without_the_sheets_tiles_an_off_tile_pair_has_no_region(self):
        mate = A.Alignment(481, 60, 60)
        plain = C.PairClassifier([TILE_1, TILE_12], 1, 0.9)
        self.assertEqual(plain.classify(mate, mate), C.Assignment("off_tile"))
        self.assertEqual(plain.classify(mate, mate).region, "")


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

    def random_pair(self, rng, tiles):
        first, second = self.random_mate(rng, tiles), self.random_mate(rng, tiles)
        if first is not None and second is not None and rng.random() < 0.7:      # most pairs have both mates in one place
            second = A.Alignment(first.position + rng.randint(-30, 30), second.length, second.span)
        return first, second

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
            first, second = self.random_pair(rng, tiles)
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

    def test_regions_never_change_the_category_or_the_tiles(self):
        """The region is a label on off-tile pairs and nothing more: with all the sheet's tiles given for labelling, every
        pair gets exactly the category and the tiles it gets without them, and only off-tile pairs carry a region."""
        job_tiles = [self.sheet.tiles["1"], self.sheet.tiles["12"]]
        everywhere = list(self.sheet.tiles.values())
        plain = C.PairClassifier(job_tiles, self.cds_start, 0.9)
        labelled = C.PairClassifier(job_tiles, self.cds_start, 0.9, regions=everywhere)
        rng = random.Random(3)
        seen = set()
        for _ in range(8000):
            first, second = self.random_pair(rng, everywhere)
            without, with_regions = plain.classify(first, second), labelled.classify(first, second)
            self.assertEqual((with_regions.category, with_regions.tiles), (without.category, without.tiles))
            self.assertEqual(bool(with_regions.region), with_regions.category == "off_tile")
            seen.add(with_regions.region)
        self.assertTrue({"", "mixed", "outside_tiles"} <= seen and len(seen) > 6, f"regions seen: {sorted(seen)}")

    def test_every_category_turns_up_so_the_agreement_is_not_vacuous(self):
        tiles = [Tile("a", 2, 13), Tile("b", 6, 20)]
        seen = set()
        for seed in range(3):
            seen |= self.compare(tiles, 0.6, "drop", seed=seed)
        self.assertEqual(seen, {"assigned", "off_tile", "ambiguous", "r1_unmapped", "r2_unmapped", "both_unmapped"})


if __name__ == "__main__":
    unittest.main()