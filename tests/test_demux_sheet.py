# These tests were generated with Claude (Anthropic) and then reviewed, edited and run by the maintainers.
"""
Tests for the tileseq_demux sheet reader (ticket 1.3): reading the parameter sheet and its Index Samples section, and
reporting every problem at once. The FASTQ folder and tool checks are in test_demux_validation.py.

Run with:  python -m pytest tests/test_demux_sheet.py -v     (from the repository root)
Without pytest (for example in the cluster's py37 environment):
           python -m unittest discover -s tests -p "test_demux_*.py"
Set TILESEQ_TEST_SHEET to a real sheet to also check that it parses. The sheets here are built in code, because
.gitignore excludes *.csv.
"""
import os
import re
import random
import tempfile
import unittest

from TileSeqDemux import sheet as S
from TileSeqDemux.rules import Facts, RuleSet

TILES = [(1, 2, 13), (2, 14, 49), (3, 50, 87), (4, 88, 123), (5, 124, 160), (6, 161, 195), (7, 196, 229),
         (8, 230, 264), (9, 265, 301), (10, 302, 337), (11, 338, 372), (12, 373, 396)]


def make_sheet(sequence=None, index=True, extra_before_index=",,,,"):
    """A sheet shaped like the lab's AAT.csv template: 4 samples, two sharing an index each."""
    rng = random.Random(1)
    seq = sequence if sequence is not None else "".join(rng.choice("acgt") for _ in range(1400))
    lines = ["Project name:,test,,,", ",,,,", "Template construct,,,,", "Gene name:,GENE,,,", "Sequence:,%s,,," % seq,
             "CDS start:,111,,,", "CDS end:,1328,,,", ",,,,", "Sequencing Tiles,,,,", "Tile Number,Start AA,End AA,,"]
    lines += ["%d,%d,%d,," % t for t in TILES]
    lines += [",,,,", "Sequencing samples,,,,", "Sample ID,Tile ID,Condition,Time point,Replicate",
              "97,1,QC,1,1", "193,12,QC,1,1", "98,2,QC,1,1", "194,11,QC,1,1", ",,,,"]
    if index:
        lines += [extra_before_index, "Index Samples,,,,", "Shared Sample ID,Unique Sample ID,,,", "97,193,,,", "98,194,,,", ",,,,"]
    lines += ["Variant calling parameters,,,,", "Posterior threshold:,0.9,,,", "Minimum coverage:,0.9,,,", ",,,,"]
    return "\n".join(l for l in lines if l is not None) + "\n"


def swap(text, old, new):
    assert old in text, "test setup: %r not in the sheet" % old
    return text.replace(old, new, 1)


class ReadTests(unittest.TestCase):
    def test_reads_the_whole_sheet(self):
        sh = S.SheetReader().parse(make_sheet())
        self.assertEqual(len(sh.sequence), 1400)
        self.assertTrue(sh.sequence.isupper())
        self.assertEqual((sh.cds_start, sh.cds_end, sh.min_cover), (111, 1328, 0.9))
        self.assertEqual(sorted(sh.samples, key=int), ["97", "98", "193", "194"])
        self.assertEqual([(g.shared, g.unique) for g in sh.groups], [("97", ("193",)), ("98", ("194",))])
        self.assertEqual(sh.samples["193"].tile, "12")
        self.assertEqual(sh.warnings, [])

    def test_tile_geometry_matches_the_runs(self):
        """Numbers printed by prototype_gene.py on the real sheet (minimum coverage 0.9, CDS start 111)."""
        sh = S.SheetReader().parse(make_sheet())
        got = {n: (t.begins, t.ends, t.min_map(0.9)) for n, t in sh.tiles.items()}
        self.assertEqual(got["1"], (4, 39, 32))
        self.assertEqual(got["12"], (1117, 1188, 64))
        self.assertEqual(got["6"], (481, 585, 94))
        self.assertEqual(got["7"], (586, 687, 91))
        self.assertEqual(got["2"], (40, 147, 97))
        self.assertEqual(got["11"], (1012, 1116, 94))
        self.assertEqual(got["9"], (793, 903, 99))

    def test_default_minimum_coverage_when_the_row_is_missing(self):
        sh = S.SheetReader().parse(swap(make_sheet(), "Minimum coverage:,0.9,,,\n", ""))
        self.assertEqual(sh.min_cover, 0.6)
        self.assertTrue(any("0.6" in w for w in sh.warnings))

    def test_crlf_bom_quotes_and_unpadded_rows(self):
        text = make_sheet().replace("\n", "\r\n")
        text = swap(text, '"Index Samples,,,,', '"Index Samples,,,,') if '"Index' in text else text
        text = text.replace("97,193,,,", '"97", 193')               # quoted cell, spaces, short row
        with tempfile.NamedTemporaryFile("w", suffix=".csv", delete=False, newline="", encoding="utf-8-sig") as fh:
            fh.write(text)
        try:
            sh = S.SheetReader().read(fh.name)
        finally:
            os.unlink(fh.name)
        self.assertEqual(sh.groups[0].unique, ("193",))

    def test_a_sheet_without_the_section_has_nothing_to_demultiplex(self):
        sheet = S.SheetReader().parse(make_sheet(index=False))
        self.assertEqual(sheet.groups, [])
        self.assertFalse(sheet.needs_demultiplexing)
        self.assertEqual((sheet.shared_ids, sheet.unique_ids), ([], []))
        self.assertTrue(any("nothing to demultiplex" in warning for warning in sheet.warnings), sheet.warnings)
        self.assertEqual(len(sheet.tiles), 12)                      # the rest of the sheet is still read and checked
        self.assertIn("nothing to demultiplex", " ".join(sheet.describe()))

    def test_a_sheet_with_the_section_needs_demultiplexing(self):
        self.assertTrue(S.SheetReader().parse(make_sheet()).needs_demultiplexing)

    def test_index_section_can_sit_at_the_end_of_the_file(self):
        text = make_sheet(index=False) + ",,,,\nIndex Samples,,,,\nShared Sample ID,Unique Sample ID,,,\n97,193,,,\n98,194"
        sh = S.SheetReader().parse(text)
        self.assertEqual(len(sh.groups), 2)

    def test_three_tiles_can_share_one_index(self):
        text = swap(make_sheet(), "97,193,,,", "97,193,,,\n97,195,,,")
        text = swap(text, "194,11,QC,1,1", "194,11,QC,1,1\n195,5,QC,1,1")
        sh = S.SheetReader().parse(text)
        self.assertEqual(sh.groups[0].members, ("97", "193", "195"))

    def test_shared_and_unique_ids(self):
        sh = S.SheetReader().parse(make_sheet())
        self.assertEqual(sh.shared_ids, ["97", "98"])
        self.assertEqual(sh.unique_ids, ["193", "194"])

    def test_describe_lists_every_member(self):
        lines = S.SheetReader().parse(make_sheet()).describe()
        self.assertIn("minimum coverage 0.9, CDS start 111", lines[0])
        self.assertTrue(any("sample 193" in l and "CDS 1117-1188" in l and "64" in l for l in lines))
        self.assertEqual(len(S.SheetReader().parse(make_sheet()).describe(0.6)), len(lines))


class ProblemTests(unittest.TestCase):
    def problems(self, text):
        with self.assertRaises(S.SheetError) as ctx:
            S.SheetReader().parse(text)
        return ctx.exception.problems

    def assertProblem(self, text, fragment):
        found = self.problems(text)
        self.assertTrue(any(fragment in p for p in found), "no problem containing %r in %s" % (fragment, found))

    def test_missing_sequence(self):
        self.assertProblem(swap(make_sheet(), "Sequence:", "Seq:"), "'Sequence:'")

    def test_label_with_the_wrong_case_is_explained(self):
        self.assertProblem(swap(make_sheet(), "Index Samples,,,,", "index samples,,,,"), "labels must match exactly")

    def test_bad_bases(self):
        self.assertProblem(make_sheet(sequence="ACGTXACGT" * 200), "other than A, C, G, T, N")

    def test_cds_start_beyond_the_sequence(self):
        self.assertProblem(swap(make_sheet(), "CDS start:,111", "CDS start:,5000"), "has only 1400 bases")

    def test_cds_end_before_start(self):
        self.assertProblem(swap(make_sheet(), "CDS end:,1328", "CDS end:,100"), "is before")

    def test_tile_beyond_the_cds(self):
        text = swap(make_sheet(), "12,373,396,,", "12,373,500,,")
        self.assertProblem(text, "tile 12 ends at CDS position 1500")

    def test_tile_with_end_before_start(self):
        self.assertProblem(swap(make_sheet(), "3,50,87,,", "3,87,50,,"), "need 1 <= start <= end")

    def test_duplicate_tile(self):
        self.assertProblem(swap(make_sheet(), "4,88,123,,", "3,88,123,,"), "tile 3 is listed twice")

    def test_sample_with_an_unknown_tile(self):
        self.assertProblem(swap(make_sheet(), "98,2,QC,1,1", "98,99,QC,1,1"), "tile '99', which is not in the tile table")

    def test_duplicate_sample(self):
        self.assertProblem(swap(make_sheet(), "194,11,QC,1,1", "98,11,QC,1,1"), "sample 98 is listed twice")

    def test_underscore_in_a_sample_id(self):
        self.assertProblem(swap(make_sheet(), "194,11,QC,1,1", "19_4,11,QC,1,1"), "must not be empty or contain '_'")

    def test_a_missing_section_does_not_hide_other_problems(self):
        problems = self.problems(swap(make_sheet(index=False), "CDS start:,111,,,", "CDS start:,,,,"))
        self.assertTrue(any("has no value" in p for p in problems))
        self.assertFalse(any("Index Samples" in p for p in problems))

    def test_a_misspelled_section_label_is_an_error_not_a_no_op(self):
        problems = self.problems(swap(make_sheet(), "Index Samples,,,,", "index samples,,,,"))
        self.assertTrue(any(p.startswith("no 'Index Samples' section (found 'index samples' on line") for p in problems), problems)
        self.assertFalse(any("nothing to demultiplex" in p for p in problems))

    def test_row_above_the_section_must_be_empty_cells(self):
        text = make_sheet(extra_before_index="Replicate notes,x,,,")
        self.assertProblem(text, "must be a row of empty cells")

    def test_wrong_index_header(self):
        self.assertProblem(swap(make_sheet(), "Shared Sample ID,Unique Sample ID", "Physical,Logical"), "the header under 'Index Samples'")

    def test_index_samples_missing_from_the_sample_table(self):
        self.assertProblem(swap(make_sheet(), "98,194,,,", "98,555,,,"), "unique sample 555 is not in the sample table")
        self.assertProblem(swap(make_sheet(), "97,193,,,", "999,193,,,"), "shared sample 999 is not in the sample table")

    def test_same_sample_twice_in_a_pair(self):
        self.assertProblem(swap(make_sheet(), "98,194,,,", "98,98,,,"), "are both 98")

    def test_pair_listed_twice(self):
        self.assertProblem(swap(make_sheet(), "98,194,,,", "97,193,,,"), "is listed twice")

    def test_unique_sample_paired_with_two_shared_samples(self):
        self.assertProblem(swap(make_sheet(), "98,194,,,", "98,193,,,"), "is already paired with shared sample 97")

    def test_unique_sample_that_is_also_shared(self):
        text = swap(make_sheet(), "98,194,,,", "193,194,,,")
        self.assertProblem(text, "unique sample here and a shared sample elsewhere")

    def test_two_samples_in_one_index_with_the_same_tile(self):
        text = swap(make_sheet(), "193,12,QC,1,1", "193,1,QC,1,1")
        self.assertProblem(text, "tiles in one index must differ")

    def test_empty_unique_cell(self):
        self.assertProblem(swap(make_sheet(), "98,194,,,", "98,,,,"), "has no unique sample ID")

    def test_minimum_coverage_out_of_range_or_not_a_number(self):
        self.assertProblem(swap(make_sheet(), "Minimum coverage:,0.9", "Minimum coverage:,1.5"), "above 0 and at most 1")
        self.assertProblem(swap(make_sheet(), "Minimum coverage:,0.9", "Minimum coverage:,high"), "should be a number")

    def test_required_value_missing(self):
        self.assertProblem(swap(make_sheet(), "CDS start:,111,,,", "CDS start:,,,,"), "'CDS start:' on line")
        self.assertProblem(swap(make_sheet(), "CDS start:,111,,,", "CDS start:,,,,"), "has no value")

    def test_empty_tile_table(self):
        text = re.sub(r"(?m)^\d+,\d+,\d+,,\n", "", make_sheet())
        self.assertProblem(text, "the tile table under 'Tile Number' is empty")

    def test_empty_sample_table(self):
        text = make_sheet()
        for row in ("97,1,QC,1,1\n", "193,12,QC,1,1\n", "98,2,QC,1,1\n", "194,11,QC,1,1\n"):
            text = swap(text, row, "")
        self.assertProblem(text, "the sample table under 'Sample ID' is empty")

    def test_wrong_second_column_in_the_sample_table(self):
        self.assertProblem(swap(make_sheet(), "Sample ID,Tile ID,", "Sample ID,Tile,"), "the second column should be 'Tile ID', found 'Tile'")

    def test_index_section_with_a_header_but_no_rows(self):
        self.assertProblem(swap(make_sheet(), "97,193,,,\n98,194,,,\n", ""), "the 'Index Samples' section lists no samples")

    def test_all_problems_are_reported_together(self):
        text = swap(make_sheet(), "Sequence:", "Seq:")
        text = swap(text, "98,194,,,", "98,555,,,")
        text = swap(text, "3,50,87,,", "3,87,50,,")
        self.assertGreaterEqual(len(self.problems(text)), 3)

    def test_error_message_lists_the_problems(self):
        with self.assertRaises(S.SheetError) as ctx:
            S.SheetReader().parse(swap(make_sheet(), "CDS start:,111,,,", "CDS start:,,,,"), path="x.csv")
        self.assertIn("in x.csv", str(ctx.exception))
        self.assertIn("  - ", str(ctx.exception))


def lookup(count, hint=""):
    return Facts(label="Sequence:", what="row", count=count, lines="3, 9", hint=hint)


def integer(text, minimum=1):
    return Facts(label="CDS start:", text=text, number=S._whole_number(text), minimum=minimum)


def cover(text, count=1):
    return Facts(label="Minimum coverage:", count=count, text=text, number=S._decimal_number(text), default=0.6)


def layout(**changes):
    facts = dict(label="Index Samples", title_line=30, above="", above_has_content=False, has_header=True, header_line=31,
                 found="Shared Sample ID,Unique Sample ID", header_ok=True, expected="Shared Sample ID,Unique Sample ID")
    facts.update(changes)
    return Facts(**facts)


def cds(start, end=None, sequence_length=1500, cds_length=1000):
    return Facts(start_label="CDS start:", end_label="CDS end:", start=start, end=end, sequence_length=sequence_length, cds_length=cds_length)


class RuleTests(unittest.TestCase):
    """Every rule fires on a subject built to break it, with its own name and a message that formats."""

    def cases(self):
        """(rule set, subject, context, the rules it breaks, whether the subject is blocked)"""
        sample_context = S.SampleContext(known_samples={"97": None}, tiles={"1": S.Tile("1", 2, 13)})
        pair_context = S.PairContext(seen={("97", "193")}, owners={"193": ("96", 4)})
        no_row_above = dict(above_has_content=False)
        return [
            (S.TILE_RULES, S.TileRow("1", "x", "13"), {}, ["not_numbers"], True),
            (S.TILE_RULES, S.TileRow("1", "5", "3"), {}, ["bad_range"], True),
            (S.TILE_RULES, S.TileRow("1", "2", "13"), {"1": None}, ["duplicate"], True),
            (S.TILE_RULES, S.TileRow("2", "2", "13"), {"1": None}, [], False),
            (S.SAMPLE_RULES, S.SampleRow("a_b", "1", "", "", ""), sample_context, ["bad_id"], True),
            (S.SAMPLE_RULES, S.SampleRow("97", "1", "", "", ""), sample_context, ["duplicate"], True),
            (S.SAMPLE_RULES, S.SampleRow("98", "9", "", "", ""), sample_context, ["unknown_tile"], False),
            (S.SAMPLE_RULES, S.SampleRow("98", "1", "", "", ""), sample_context, [], False),
            (S.PAIR_RULES, S.PairRow("97", ""), pair_context, ["no_unique"], True),
            (S.PAIR_RULES, S.PairRow("9_7", "194"), pair_context, ["bad_id"], True),
            (S.PAIR_RULES, S.PairRow("97", "97"), pair_context, ["same_id"], True),
            (S.PAIR_RULES, S.PairRow("97", "193"), pair_context, ["pair_twice"], True),
            (S.PAIR_RULES, S.PairRow("98", "193"), pair_context, ["unique_reused"], True),
            (S.PAIR_RULES, S.PairRow("98", "194"), pair_context, [], False),
            (S.LOOKUP_RULES, lookup(0), None, ["missing"], True),
            (S.LOOKUP_RULES, lookup(2), None, ["repeated"], True),
            (S.LOOKUP_RULES, lookup(1), None, [], False),
            (S.INDEX_LOOKUP_RULES, lookup(0), None, ["absent"], True),
            (S.INDEX_LOOKUP_RULES, lookup(0, hint=" (found 'index samples' on line 31; labels must match exactly)"), None, ["misspelled"], True),
            (S.INDEX_LOOKUP_RULES, lookup(2), None, ["repeated"], True),
            (S.INDEX_LOOKUP_RULES, lookup(1), None, [], False),
            (S.INTEGER_RULES, integer("abc"), None, ["not_whole"], True),
            (S.INTEGER_RULES, integer("0"), None, ["too_small"], True),
            (S.INTEGER_RULES, integer("111"), None, [], False),
            (S.COVER_RULES, cover("0.9", count=2), None, ["repeated"], True),
            (S.COVER_RULES, cover(""), None, ["no_value"], False),
            (S.COVER_RULES, cover("high"), None, ["not_a_number"], True),
            (S.COVER_RULES, cover("1.5"), None, ["out_of_range"], False),
            (S.COVER_RULES, cover("0.9"), None, [], False),
            (S.INDEX_SECTION_RULES, layout(above="x,y", above_has_content=True), None, ["row_above_not_empty"], False),
            (S.INDEX_SECTION_RULES, layout(has_header=False, header_line=None, found="an empty row", header_ok=False), None, ["no_header_row"], True),
            (S.INDEX_SECTION_RULES, layout(header_ok=False, found="a,b"), None, ["wrong_header"], True),
            (S.INDEX_SECTION_RULES, layout(**no_row_above), None, [], False),
            (S.CDS_RULES, cds(start=2000), None, ["start_beyond_sequence"], True),
            (S.CDS_RULES, cds(start=100, end=50), None, ["end_before_start"], True),
            (S.CDS_RULES, cds(start=100, end=2000), None, ["end_beyond_sequence"], True),
            (S.CDS_RULES, cds(start=100, end=1300, cds_length=1201), None, ["length_not_multiple_of_3"], False),
            (S.CDS_RULES, cds(start=100, end=1299, cds_length=1200), None, [], False),
            (S.CDS_RULES, cds(start=100, end=None, cds_length=1401), None, [], False),
        ]

    def test_each_subject_breaks_exactly_the_rules_it_was_built_to_break(self):
        for rule_set, subject, context, names, blocked in self.cases():
            with self.subTest(subject=subject):
                verdict = rule_set.check(subject, context)
                self.assertEqual((verdict.names, verdict.blocked), (names, blocked))
                self.assertEqual(len(verdict.problems) + len(verdict.warnings), len(names))

    def test_every_rule_in_every_table_is_exercised_by_a_case_above(self):
        """Adding a Rule, or a whole RuleSet, without adding a case for it fails here."""
        defined = {(value, rule.name) for value in vars(S).values() if isinstance(value, RuleSet) for rule in value.rules}
        exercised = {(rule_set, name) for rule_set, _subject, _context, names, _blocked in self.cases() for name in names}
        self.assertEqual(exercised, defined)
        self.assertGreater(len({rule_set for rule_set, _name in defined}), 5)

    def test_the_unknown_tile_message_lists_the_tiles_in_numeric_order(self):
        context = S.SampleContext(known_samples={}, tiles={"10": None, "2": None, "1": None})
        message = S.SAMPLE_RULES.check(S.SampleRow("98", "9", "", "", ""), context).problems[0]
        self.assertTrue(message.endswith("(tiles: 1, 2, 10)"))

    def test_a_pair_with_two_bad_ids_reports_the_first_one_only(self):
        verdict = S.PAIR_RULES.check(S.PairRow("9_7", "1_93"), S.PairContext(seen=set(), owners={}))
        self.assertEqual(len(verdict.problems), 1)
        self.assertIn("'9_7'", verdict.problems[0])

    def test_the_reused_unique_message_names_the_first_row_that_used_it(self):
        context = S.PairContext(seen=set(), owners={"193": ("96", 4)})
        message = S.PAIR_RULES.check(S.PairRow("98", "193"), context).problems[0]
        self.assertEqual(message, "unique sample 193 is already paired with shared sample 96 on line 4")

    def test_an_absent_section_is_a_warning_that_stops_the_reading_not_a_problem(self):
        verdict = S.INDEX_LOOKUP_RULES.check(lookup(0))
        self.assertEqual((verdict.problems, verdict.blocked), ([], True))
        self.assertTrue(verdict.warnings[0].endswith("no samples share an index, so there is nothing to demultiplex"))

    def test_a_misspelled_section_label_is_a_problem_with_a_hint(self):
        verdict = S.INDEX_LOOKUP_RULES.check(lookup(0, hint=" (found 'index samples' on line 31; labels must match exactly)"))
        self.assertEqual((verdict.warnings, verdict.blocked), ([], True))
        self.assertEqual(verdict.problems, ["no 'Sequence:' section (found 'index samples' on line 31; labels must match exactly)"])

    def test_the_cds_length_check_is_a_warning_not_a_problem(self):
        verdict = S.CDS_RULES.check(cds(start=100, end=1300, cds_length=1201))
        self.assertEqual(verdict.problems, [])
        self.assertEqual(verdict.warnings, ["the CDS (100 to 1300) is 1201 bases long, not a multiple of 3"])

    def test_a_missing_minimum_coverage_is_a_warning_that_names_the_default(self):
        verdict = S.COVER_RULES.check(cover(""))
        self.assertEqual((verdict.problems, verdict.blocked), ([], False))
        self.assertEqual(verdict.warnings, ["no 'Minimum coverage:' value: using 0.6, as tsm csv2json does"])


@unittest.skipUnless(os.environ.get("TILESEQ_TEST_SHEET"), "set TILESEQ_TEST_SHEET to check a real sheet")
class RealSheetTest(unittest.TestCase):
    def test_real_sheet_parses(self):
        sh = S.SheetReader().read(os.environ["TILESEQ_TEST_SHEET"])
        self.assertTrue(sh.groups)


if __name__ == "__main__":
    unittest.main()