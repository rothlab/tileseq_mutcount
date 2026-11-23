# These tests were generated with Claude (Anthropic) and then reviewed, edited and run by the maintainers.
"""
Tests for the tileseq_demux input validation (ticket 1.3): ToolChecker, FastqFolderChecker and the
DemuxInputsValidator that runs them together with the sheet reader.

Run with:  python -m pytest tests/test_demux_validation.py -v     (from the repository root)
Without pytest: python -m unittest discover -s tests -p "test_demux_*.py"
"""
import contextlib
import io
import os
import stat
import subprocess
import tempfile
import unittest
from unittest import mock

from TileSeqDemux import validation as V
from TileSeqDemux.sheet import SheetReader
from test_demux_sheet import make_sheet      # tests/ is on sys.path under both pytest and unittest


def done(returncode=0, out=""):
    return subprocess.CompletedProcess([], returncode, stdout=out)


def fake_run(table):
    """table maps a command tuple to a CompletedProcess or an exception to raise."""
    def run(cmd, timeout=60):
        result = table[tuple(cmd)]
        if isinstance(result, Exception):
            raise result
        return result
    return run


def which_in(*present):
    return lambda tool: "/usr/bin/" + tool if tool in present else None


GOOD_TOOLS = {("bowtie2", "--version"): done(0, "bowtie2-align-s version 2.5.5\n"),
              ("bowtie2-build", "--version"): done(0, "bowtie2-build-s version 2.5.5\n")}


def good_tool_checker():
    return V.ToolChecker(run=fake_run(GOOD_TOOLS), which=which_in("bowtie2", "bowtie2-build"))


def missing_tool_checker():
    return V.ToolChecker(run=fake_run(GOOD_TOOLS), which=which_in())


class ToolCheckerTests(unittest.TestCase):
    def test_everything_present(self):
        self.assertEqual(good_tool_checker().check(), [])

    def test_missing_tools_are_all_listed(self):
        problems = missing_tool_checker().check()
        self.assertEqual(len(problems), 2)
        self.assertIn("bowtie2 was not found on PATH", problems[0])

    def test_tool_that_fails_when_asked_for_its_version(self):
        table = dict(GOOD_TOOLS)
        table[("bowtie2", "--version")] = done(127, "error while loading shared libraries: libtbb.so.2\n")
        problems = V.ToolChecker(run=fake_run(table), which=which_in("bowtie2", "bowtie2-build")).check()
        self.assertEqual(len(problems), 1)
        self.assertIn("libtbb", problems[0])

    def test_tool_that_cannot_be_executed(self):
        table = dict(GOOD_TOOLS)
        table[("bowtie2", "--version")] = OSError("Exec format error")
        problems = V.ToolChecker(run=fake_run(table), which=which_in("bowtie2", "bowtie2-build")).check()
        self.assertIn("Exec format error", problems[0])

    def test_tool_that_times_out(self):
        table = dict(GOOD_TOOLS)
        table[("bowtie2", "--version")] = subprocess.TimeoutExpired("bowtie2", 60)
        self.assertEqual(len(V.ToolChecker(run=fake_run(table), which=which_in("bowtie2", "bowtie2-build")).check()), 1)

    def test_the_last_line_of_the_output_is_the_one_reported(self):
        """Whatever warnings come first, on any machine, the error is usually the last line."""
        table = dict(GOOD_TOOLS)
        table[("bowtie2", "--version")] = done(1, "some warning\n\nanother warning\nreal error here\n\n")
        problem = V.ToolChecker(run=fake_run(table), which=which_in("bowtie2", "bowtie2-build")).check()[0]
        self.assertTrue(problem.endswith("real error here"))
        self.assertNotIn("some warning", problem)

    def test_a_tool_that_fails_silently_is_described(self):
        table = dict(GOOD_TOOLS)
        table[("bowtie2", "--version")] = done(2, "")
        problem = V.ToolChecker(run=fake_run(table), which=which_in("bowtie2", "bowtie2-build")).check()[0]
        self.assertTrue(problem.endswith("no output"))

    def test_a_different_list_of_tools(self):
        self.assertEqual(V.ToolChecker(tools=("samtools",), run=fake_run({}), which=which_in()).check(),
                         ["samtools was not found on PATH."])

    def test_real_executables_on_a_path(self):
        """No fakes: scripts in a temporary directory that is the only thing on PATH."""
        with tempfile.TemporaryDirectory() as d:
            for name in ("bowtie2", "bowtie2-build"):
                path = os.path.join(d, name)
                with open(path, "w") as fh:
                    fh.write("#!/bin/sh\necho '%s version 2.5.5'\n" % name)
                os.chmod(path, os.stat(path).st_mode | stat.S_IEXEC)
            with mock.patch.dict(os.environ, {"PATH": d}):
                self.assertEqual(V.ToolChecker().check(), [])
            with mock.patch.dict(os.environ, {"PATH": os.path.join(d, "empty")}):
                self.assertEqual(len(V.ToolChecker().check()), 2)


class ToolRuleTests(unittest.TestCase):
    """Every tool rule fires on a probe built to break it, and the probe keeps what the tool printed."""
    FOUND = "/usr/bin/bowtie2"

    def cases(self):
        return [
            (V.ToolProbe("bowtie2", None), ["not_found"]),
            (V.ToolProbe("bowtie2", self.FOUND, error="Exec format error"), ["cannot_run"]),
            (V.ToolProbe("bowtie2", self.FOUND, status=127, output="a warning\nlibrary missing"), ["bad_exit"]),
            (V.ToolProbe("bowtie2", self.FOUND, status=0, output="bowtie2 version 2.5.5"), []),
        ]

    def test_each_probe_breaks_exactly_the_rules_it_was_built_to_break(self):
        for probe, names in self.cases():
            with self.subTest(probe=probe):
                verdict = V.TOOL_RULES.check(probe, None)
                self.assertEqual(verdict.names, names)
                self.assertEqual(verdict.blocked, bool(names))

    def test_every_tool_rule_is_exercised_by_a_case_above(self):
        """Adding a Rule without adding a case for it fails here."""
        exercised = {name for _probe, names in self.cases() for name in names}
        self.assertEqual(exercised, {rule.name for rule in V.TOOL_RULES.rules})

    def test_the_bad_exit_message_shows_the_last_line_of_the_output(self):
        message = V.TOOL_RULES.check(self.cases()[2][0]).problems[0]
        self.assertEqual(message, "bowtie2 (/usr/bin/bowtie2) exited with status 127 when asked for its version: library missing")

    def test_a_probe_keeps_the_version_the_tool_printed(self):
        """The run record in ticket 1.4 can take the version from the probe without running the tool again."""
        checker = V.ToolChecker(run=fake_run(GOOD_TOOLS), which=which_in("bowtie2", "bowtie2-build"))
        probe = checker._probe("bowtie2")
        self.assertEqual((probe.path, probe.status), ("/usr/bin/bowtie2", 0))
        self.assertIn("2.5.5", probe.output)

    def test_a_tool_that_is_not_on_path_is_not_run(self):
        checker = V.ToolChecker(run=fake_run({}), which=which_in())      # fake_run({}) would raise KeyError if called
        self.assertIsNone(checker._probe("bowtie2").path)


class FolderTests(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.TemporaryDirectory()
        self.sheet = SheetReader().parse(make_sheet())
        self.checker = V.FastqFolderChecker()

    def tearDown(self):
        self.dir.cleanup()

    def touch(self, *names):
        for n in names:
            open(os.path.join(self.dir.name, n), "w").close()

    def test_complete_folder_passes(self):
        self.touch("97_S97_R1_001.fastq.gz", "97_S97_R2_001.fastq.gz", "98_S98_R1_001.fastq.gz", "98_S98_R2_001.fastq.gz")
        self.assertEqual(self.checker.check(self.sheet, self.dir.name), [])

    def test_missing_mate_is_reported(self):
        self.touch("97_S97_R1_001.fastq.gz", "98_S98_R1_001.fastq.gz", "98_S98_R2_001.fastq.gz")
        problems = self.checker.check(self.sheet, self.dir.name)
        self.assertEqual(len(problems), 1)
        self.assertIn("shared sample 97: expected one R2 file", problems[0])

    def test_two_files_for_one_mate_are_reported(self):
        self.touch("97_S97_R1_001.fastq.gz", "97_S1_R1_001.fastq.gz", "97_S97_R2_001.fastq.gz", "98_S98_R1_001.fastq.gz", "98_S98_R2_001.fastq.gz")
        self.assertIn("found 2", self.checker.check(self.sheet, self.dir.name)[0])

    def test_existing_file_for_a_unique_sample_would_collide(self):
        self.touch("97_S97_R1_001.fastq.gz", "97_S97_R2_001.fastq.gz", "98_S98_R1_001.fastq.gz", "98_S98_R2_001.fastq.gz", "193_S193_R1_001.fastq.gz")
        self.assertIn("would collide", self.checker.check(self.sheet, self.dir.name)[0])

    def test_unrelated_files_are_ignored(self):
        self.touch("97_S97_R1_001.fastq.gz", "97_S97_R2_001.fastq.gz", "98_S98_R1_001.fastq.gz", "98_S98_R2_001.fastq.gz",
                   "Undetermined_S0_R1_001.fastq.gz", "README", "970_S970_R1_001.fastq.gz")
        self.assertEqual(self.checker.check(self.sheet, self.dir.name), [])

    def test_a_sample_id_is_the_whole_first_token(self):
        """97 must not match 970, and a file that is not a gzipped FASTQ does not count as a mate."""
        self.touch("970_S970_R1_001.fastq.gz", "970_S970_R2_001.fastq.gz", "97_S97_R1_001.fastq.gz", "97_S97_R2_001.fastq", "98_S98_R1_001.fastq.gz", "98_S98_R2_001.fastq.gz")
        problems = self.checker.check(self.sheet, self.dir.name)
        self.assertEqual(len(problems), 1)
        self.assertIn("shared sample 97: expected one R2 file", problems[0])

    def test_all_mate_problems_come_before_clashes(self):
        self.touch("193_S193_R1_001.fastq.gz")           # a clash, and no FASTQs at all for 97 and 98
        problems = self.checker.check(self.sheet, self.dir.name)
        self.assertEqual(len(problems), 5)
        self.assertTrue(all(p.startswith("shared sample") for p in problems[:4]))
        self.assertTrue(problems[4].startswith("unique sample 193"))

    def test_missing_folder(self):
        self.assertIn("does not exist", self.checker.check(self.sheet, os.path.join(self.dir.name, "nope"))[0])


class ValidatorTests(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.TemporaryDirectory()
        self.fastq = os.path.join(self.dir.name, "fastq")
        os.mkdir(self.fastq)
        self.sheet_path = os.path.join(self.dir.name, "sheet.csv")

    def tearDown(self):
        self.dir.cleanup()

    def write_sheet(self, text):
        with open(self.sheet_path, "w") as fh:
            fh.write(text)

    def write_fastq(self, *names):
        for n in names:
            open(os.path.join(self.fastq, n), "w").close()

    def good_inputs(self):
        self.write_sheet(make_sheet())
        self.write_fastq("97_S97_R1_001.fastq.gz", "97_S97_R2_001.fastq.gz", "98_S98_R1_001.fastq.gz", "98_S98_R2_001.fastq.gz")

    def validator(self, tools=None):
        return V.DemuxInputsValidator(self.sheet_path, self.fastq, tool_checker=tools or good_tool_checker())

    def headings(self, err):
        return [h for h, _ in err.sections]

    def test_everything_fine_returns_the_sheet(self):
        self.good_inputs()
        sheet = self.validator().validate()
        self.assertEqual([g.shared for g in sheet.groups], ["97", "98"])

    def test_missing_tools_are_reported_with_the_hint(self):
        self.good_inputs()
        with self.assertRaises(V.InputValidationError) as ctx:
            self.validator(missing_tool_checker()).validate()
        self.assertEqual(self.headings(ctx.exception), ["Tools"])
        self.assertIn("bowtie2 was not found", str(ctx.exception))
        self.assertIn("inherits the environment", str(ctx.exception))

    def test_tool_and_sheet_problems_arrive_together(self):
        self.write_sheet(make_sheet(index=False))
        with self.assertRaises(V.InputValidationError) as ctx:
            self.validator(missing_tool_checker()).validate()
        heads = self.headings(ctx.exception)
        self.assertEqual(heads[0], "Tools")
        self.assertTrue(heads[1].startswith("Parameter sheet"))
        self.assertIn("nothing to demultiplex", str(ctx.exception))
        self.assertIn("not checked until the parameter sheet is valid", str(ctx.exception))

    def test_tool_and_folder_problems_arrive_together(self):
        self.write_sheet(make_sheet())                      # good sheet, empty FASTQ folder
        with self.assertRaises(V.InputValidationError) as ctx:
            self.validator(missing_tool_checker()).validate()
        heads = self.headings(ctx.exception)
        self.assertEqual(heads[0], "Tools")
        self.assertTrue(heads[1].startswith("Input files"))
        self.assertIn("shared sample 97: expected one R1 file", str(ctx.exception))

    def test_folder_problems_alone(self):
        self.write_sheet(make_sheet())
        self.write_fastq("97_S97_R1_001.fastq.gz", "97_S97_R2_001.fastq.gz", "98_S98_R1_001.fastq.gz")
        with self.assertRaises(V.InputValidationError) as ctx:
            self.validator().validate()
        self.assertEqual(len(ctx.exception.sections), 1)
        self.assertIn("shared sample 98: expected one R2 file", str(ctx.exception))

    def test_missing_sheet_file(self):
        self.sheet_path = os.path.join(self.dir.name, "nope.csv")
        with self.assertRaises(V.InputValidationError) as ctx:
            self.validator().validate()
        self.assertIn("cannot read", str(ctx.exception))
        self.assertIn("nope.csv", str(ctx.exception))

    def test_the_message_has_a_heading_per_section(self):
        self.write_sheet(make_sheet(index=False))
        with self.assertRaises(V.InputValidationError) as ctx:
            self.validator(missing_tool_checker()).validate()
        text = str(ctx.exception)
        self.assertTrue(text.startswith("tileseq_demux cannot start."))
        self.assertIn("\nTools\n  - ", text)

    def test_warnings_survive_on_the_returned_sheet(self):
        self.write_sheet(make_sheet().replace("Minimum coverage:,0.9,,,\n", ""))
        self.write_fastq("97_S97_R1_001.fastq.gz", "97_S97_R2_001.fastq.gz", "98_S98_R1_001.fastq.gz", "98_S98_R2_001.fastq.gz")
        self.assertTrue(any("0.6" in w for w in self.validator().validate().warnings))

    def test_collaborators_can_be_replaced(self):
        """The validator only gathers what its parts report; fakes show it makes no checks of its own."""
        class FakeReader:
            def read(self, path):
                return SheetReader().parse(make_sheet(), path)

        class FakeFolder:
            def check(self, sheet, folder):
                return ["a fake folder problem"]

        validator = V.DemuxInputsValidator("x.csv", "y", tool_checker=good_tool_checker(), sheet_reader=FakeReader(), folder_checker=FakeFolder())
        with self.assertRaises(V.InputValidationError) as ctx:
            validator.validate()
        self.assertEqual(self.headings(ctx.exception), ["Input files (y)"])
        self.assertIn("a fake folder problem", str(ctx.exception))


class DeveloperCommandTests(unittest.TestCase):
    def test_exit_codes(self):
        with tempfile.TemporaryDirectory() as d:
            sheet_path = os.path.join(d, "sheet.csv")
            with open(sheet_path, "w") as fh:
                fh.write(make_sheet(index=False))
            with contextlib.redirect_stderr(io.StringIO()):          # the error text is not what is being tested
                self.assertEqual(V.main([sheet_path, d]), 1)         # no Index Samples section (and real tools may be absent too)
                self.assertEqual(V.main([os.path.join(d, "nope.csv"), d]), 1)


if __name__ == "__main__":
    unittest.main()