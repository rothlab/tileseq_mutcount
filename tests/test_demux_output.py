# These tests were generated with Claude (Anthropic) and then reviewed, edited and run by the maintainers.
"""
Tests for the output plan and the output folder (slice A of ticket 1.4): which files a run will create, what state a
folder is in, what a fresh run, --resume and --force each do with it, and how a failed or killed run is left.

Run with:  python -m pytest tests/test_demux_output.py -v     (from the repository root)
Without pytest: python -m unittest discover -s tests -p "test_demux_*.py"
"""
import json
import os
import tempfile
import unittest

from TileSeqDemux import output as O
from TileSeqDemux.sheet import SheetReader
from test_demux_sheet import make_sheet, swap      # tests/ is on sys.path under both pytest and unittest

RAW_NAMES = ["97_S97_R1_001.fastq.gz", "97_S97_R2_001.fastq.gz", "98_S98_R1_001.fastq.gz", "98_S98_R2_001.fastq.gz",
             "99_S99_R1_001.fastq.gz", "99_S99_R2_001.fastq.gz", "Undetermined_S0_R1_001.fastq.gz", "Undetermined_S0_R2_001.fastq.gz"]
TILE_97 = ["97_S97_R1_001.fastq.gz", "97_S97_R2_001.fastq.gz", "193_S193_R1_001.fastq.gz", "193_S193_R2_001.fastq.gz"]
TILE_98 = ["98_S98_R1_001.fastq.gz", "98_S98_R2_001.fastq.gz", "194_S194_R1_001.fastq.gz", "194_S194_R2_001.fastq.gz"]
LINKS = ["99_S99_R1_001.fastq.gz", "99_S99_R2_001.fastq.gz", "Undetermined_S0_R1_001.fastq.gz", "Undetermined_S0_R2_001.fastq.gz"]


class PlanTests(unittest.TestCase):
    def setUp(self):
        self.sheet = SheetReader().parse(make_sheet())

    def build(self, names=RAW_NAMES, raw_dir="/data/raw"):
        return O.OutputPlan.build(self.sheet, raw_dir, names)

    def test_every_member_of_a_group_gets_an_r1_and_an_r2_file(self):
        plan = self.build()
        self.assertEqual([entry.name for entry in plan.tile_files("97")], TILE_97)
        self.assertEqual([entry.name for entry in plan.tile_files("98")], TILE_98)

    def test_jobs_follow_the_sheet(self):
        self.assertEqual(self.build().jobs, ("97", "98"))

    def test_each_job_knows_its_raw_pair(self):
        self.assertEqual(self.build().raw_pairs["97"], ("97_S97_R1_001.fastq.gz", "97_S97_R2_001.fastq.gz"))

    def test_links_cover_every_raw_fastq_that_is_not_a_pooled_file(self):
        self.assertEqual([entry.name for entry in self.build().links], LINKS)

    def test_files_that_are_not_fastq_are_ignored(self):
        plan = self.build(RAW_NAMES + ["README", "notes.txt", "97_S97_R1_001.fastq"])
        self.assertEqual([entry.name for entry in plan.links], LINKS)

    def test_link_targets_are_absolute_paths_into_the_raw_folder(self):
        plan = self.build(raw_dir="raw_data")
        for entry in plan.links:
            self.assertTrue(os.path.isabs(entry.target))
            self.assertEqual(entry.target, os.path.join(os.path.abspath("raw_data"), entry.name))

    def test_tile_files_have_a_tmp_name_while_they_are_being_written(self):
        entry = self.build().tile_files("97")[0]
        self.assertEqual(entry.tmp_name, "97_S97_R1_001.fastq.gz.tmp")

    def test_a_raw_file_for_a_unique_sample_would_clash_with_an_output(self):
        with self.assertRaisesRegex(O.PlanError, "193_S193_R1_001.fastq.gz"):
            self.build(RAW_NAMES + ["193_S193_R1_001.fastq.gz"])

    def test_a_missing_mate_is_an_error(self):
        with self.assertRaisesRegex(O.PlanError, "shared sample 97: expected one R2 file.*found 0"):
            self.build([name for name in RAW_NAMES if name != "97_S97_R2_001.fastq.gz"])

    def test_two_files_for_one_mate_are_an_error(self):
        with self.assertRaisesRegex(O.PlanError, "shared sample 97: expected one R1 file.*found 2"):
            self.build(RAW_NAMES + ["97_S5_R1_001.fastq.gz"])

    def test_the_plan_describes_what_each_job_reads_and_writes(self):
        lines = self.build().describe()
        self.assertEqual(lines[0], "job 97: reads 97_S97_R1_001.fastq.gz and 97_S97_R2_001.fastq.gz")
        self.assertEqual(lines[1], "  writes " + ", ".join(TILE_97))
        self.assertEqual(lines[-1], "4 links to raw files (unshared samples and other FASTQ files)")

    def test_the_plan_survives_a_round_trip_through_json(self):
        plan = self.build()
        again = O.OutputPlan.from_dict(json.loads(json.dumps(plan.to_dict())))
        self.assertEqual((again.entries, again.raw_pairs, again.jobs), (plan.entries, plan.raw_pairs, plan.jobs))

    def test_a_sheet_that_needs_no_demultiplexing_has_no_jobs_and_links_everything(self):
        plan = O.OutputPlan.build(SheetReader().parse(make_sheet(index=False)), "/data/raw", RAW_NAMES)
        self.assertEqual((plan.jobs, plan.tile_entries), ((), ()))
        self.assertEqual(len(plan.links), len(RAW_NAMES))


class FolderTestCase(unittest.TestCase):
    """A raw folder with real files, an output path next to it, and helpers to put the folder in any state."""

    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.root = os.path.realpath(self.directory.name)
        self.raw = os.path.join(self.root, "raw")
        os.mkdir(self.raw)
        for name in RAW_NAMES:
            self.touch(os.path.join(self.raw, name))
        self.out = os.path.join(self.root, "demux")
        self.sheet = SheetReader().parse(make_sheet())

    def tearDown(self):
        self.directory.cleanup()

    @staticmethod
    def touch(path, text=""):
        with open(path, "w") as handle:
            handle.write(text)

    def plan(self, sheet=None):
        return O.OutputPlan.build(sheet or self.sheet, self.raw, sorted(os.listdir(self.raw)))

    def folder(self, sheet=None, path=None):
        return O.OutputFolder(path or self.out, self.plan(sheet))

    def write_tmp_files(self, folder, job):
        for entry in folder.plan.tile_files(job):
            self.touch(folder.path_of(entry.tmp_name), "reads")

    def write_job_record(self, folder, job):
        os.makedirs(os.path.dirname(folder.job_record_path(job)), exist_ok=True)
        self.touch(folder.job_record_path(job), "{}")

    def finished_job(self, folder, job):
        """Write the job's files and commit them, then write its record: what a job that succeeded leaves behind."""
        self.write_tmp_files(folder, job)
        folder.commit_job(job)
        self.write_job_record(folder, job)

    def complete_run(self):
        folder = self.folder()
        folder.prepare(O.RunMode.NEW)
        for job in folder.plan.jobs:
            self.finished_job(folder, job)
        folder.write_report(O.MANIFEST_NAME, "manifest\n")
        folder.finish()
        return folder

    def exists(self, folder, name):
        return os.path.lexists(folder.path_of(name))

    def names_in(self, folder):
        return sorted(os.listdir(folder.path))


class StateTests(FolderTestCase):
    def test_a_missing_folder_is_new(self):
        self.assertEqual(self.folder().state, O.FolderState.NEW)

    def test_an_existing_empty_folder_is_new(self):
        os.mkdir(self.out)
        self.assertEqual(self.folder().state, O.FolderState.NEW)

    def test_a_folder_with_other_files_and_no_run_record_is_foreign(self):
        os.mkdir(self.out)
        self.touch(os.path.join(self.out, "notes.txt"))
        self.assertEqual(self.folder().state, O.FolderState.FOREIGN)

    def test_a_file_in_place_of_the_folder_is_foreign(self):
        self.touch(self.out)
        self.assertEqual(self.folder().state, O.FolderState.FOREIGN)

    def test_a_prepared_folder_is_incomplete(self):
        folder = self.folder()
        folder.prepare(O.RunMode.NEW)
        self.assertEqual(folder.state, O.FolderState.INCOMPLETE)

    def test_a_finished_run_is_complete(self):
        self.assertEqual(self.complete_run().state, O.FolderState.COMPLETE)

    def test_status_complete_without_a_manifest_is_still_incomplete(self):
        folder = self.complete_run()
        os.remove(folder.manifest)
        self.assertEqual(folder.state, O.FolderState.INCOMPLETE)

    def test_a_marker_that_is_still_there_means_incomplete(self):
        folder = self.complete_run()
        self.touch(folder.marker, "left over")
        self.assertEqual(folder.state, O.FolderState.INCOMPLETE)


class ModeTests(FolderTestCase):
    """What each run mode does in each state: the table the cards promise."""

    EXPECTED = {
        (O.FolderState.NEW, O.RunMode.NEW): None, (O.FolderState.NEW, O.RunMode.RESUME): None, (O.FolderState.NEW, O.RunMode.FORCE): None,
        (O.FolderState.INCOMPLETE, O.RunMode.NEW): O.OutputFolderError, (O.FolderState.INCOMPLETE, O.RunMode.RESUME): None,
        (O.FolderState.INCOMPLETE, O.RunMode.FORCE): None,
        (O.FolderState.COMPLETE, O.RunMode.NEW): O.OutputFolderError, (O.FolderState.COMPLETE, O.RunMode.RESUME): O.AlreadyComplete,
        (O.FolderState.COMPLETE, O.RunMode.FORCE): None,
        (O.FolderState.FOREIGN, O.RunMode.NEW): O.OutputFolderError, (O.FolderState.FOREIGN, O.RunMode.RESUME): O.OutputFolderError,
        (O.FolderState.FOREIGN, O.RunMode.FORCE): O.OutputFolderError,
    }

    def put_in_state(self, state):
        if state is O.FolderState.INCOMPLETE:
            self.folder().prepare(O.RunMode.NEW)
        elif state is O.FolderState.COMPLETE:
            self.complete_run()
        elif state is O.FolderState.FOREIGN:
            os.mkdir(self.out)
            self.touch(os.path.join(self.out, "notes.txt"))

    def test_every_state_and_mode(self):
        for (state, mode), expected in self.EXPECTED.items():
            with self.subTest(state=state.value, mode=mode.value):
                self.tearDown()
                self.setUp()
                self.put_in_state(state)
                folder = self.folder()
                self.assertEqual(folder.state, state)
                if expected is None:
                    folder.prepare(mode)
                    return_state = folder.state
                    self.assertEqual(return_state, O.FolderState.INCOMPLETE)     # a started run is always unfinished
                else:
                    with self.assertRaises(expected) as caught:
                        folder.prepare(mode)
                    self.assertIs(type(caught.exception), expected)

    def test_an_unfinished_run_is_refused_with_the_count_and_both_ways_forward(self):
        self.folder().prepare(O.RunMode.NEW)
        with self.assertRaises(O.OutputFolderError) as caught:
            self.folder().prepare(O.RunMode.NEW)
        message = str(caught.exception)
        for expected in ("0 of 2 jobs complete", "--resume", "--force", O.MARKER_NAME):
            self.assertIn(expected, message)

    def test_a_foreign_folder_is_never_touched_not_even_with_force(self):
        self.put_in_state(O.FolderState.FOREIGN)
        folder = self.folder()
        for mode in O.RunMode:
            with self.assertRaises(O.OutputFolderError):
                folder.prepare(mode)
        self.assertEqual(self.names_in(folder), ["notes.txt"])

    def test_a_file_where_the_folder_should_be_is_refused(self):
        self.touch(self.out)
        with self.assertRaisesRegex(O.OutputFolderError, "not a folder"):
            self.folder().prepare(O.RunMode.NEW)

    def test_the_raw_folder_is_refused_as_an_output_before_anything_is_written(self):
        folder = self.folder(path=self.raw)
        with self.assertRaisesRegex(O.OutputFolderError, "raw FASTQ folder"):
            folder.prepare(O.RunMode.FORCE)
        self.assertEqual(sorted(os.listdir(self.raw)), sorted(RAW_NAMES))


class PrepareTests(FolderTestCase):
    def test_prepare_creates_the_folder_the_logs_the_marker_and_the_records(self):
        folder = self.folder()
        folder.prepare(O.RunMode.NEW)
        self.assertEqual(self.names_in(folder), sorted([O.LOGS_NAME, O.MARKER_NAME]))
        self.assertEqual(sorted(os.listdir(folder.logs)), sorted([O.PLAN_NAME, O.RUN_RECORD_NAME]))
        with open(folder.run_record) as handle:
            self.assertEqual(json.load(handle)["status"], "running")

    def test_the_marker_says_what_happened_and_what_to_do(self):
        folder = self.folder()
        folder.prepare(O.RunMode.NEW)
        with open(folder.marker) as handle:
            text = handle.read()
        for expected in ("did not finish", "demux.log", ".tmp", "rerun the same command with --resume", "--force"):
            self.assertIn(expected, text)

    def test_the_recorded_plan_is_the_plan_of_the_run(self):
        folder = self.folder()
        folder.prepare(O.RunMode.NEW)
        with open(folder.plan_record) as handle:
            self.assertEqual(O.OutputPlan.from_dict(json.load(handle)).entries, folder.plan.entries)

    def test_the_rerun_hint_is_the_one_the_command_gave(self):
        folder = O.OutputFolder(self.out, self.plan(), rerun_hint="tileseq_demux -p sheet.csv -f raw -o demux --resume")
        folder.prepare(O.RunMode.NEW)
        with open(folder.marker) as handle:
            self.assertIn("To finish: tileseq_demux -p sheet.csv -f raw -o demux --resume.", handle.read())

    def test_no_file_ends_up_beside_the_sheet_or_in_the_raw_folder(self):
        folder = self.folder()
        folder.prepare(O.RunMode.NEW)
        self.assertEqual(sorted(os.listdir(self.raw)), sorted(RAW_NAMES))


class CommitTests(FolderTestCase):
    def setUp(self):
        super().setUp()
        self.folder_under_test = self.folder()
        self.folder_under_test.prepare(O.RunMode.NEW)

    def test_commit_gives_a_jobs_files_their_final_names(self):
        folder = self.folder_under_test
        self.write_tmp_files(folder, "97")
        folder.commit_job("97")
        for name in TILE_97:
            self.assertTrue(self.exists(folder, name))
            self.assertFalse(self.exists(folder, name + ".tmp"))

    def test_commit_leaves_other_jobs_alone(self):
        folder = self.folder_under_test
        self.write_tmp_files(folder, "97")
        self.write_tmp_files(folder, "98")
        folder.commit_job("97")
        self.assertTrue(all(self.exists(folder, name + ".tmp") for name in TILE_98))
        self.assertFalse(any(self.exists(folder, name) for name in TILE_98))

    def test_commit_refuses_when_a_file_is_missing_and_changes_nothing(self):
        folder = self.folder_under_test
        self.write_tmp_files(folder, "97")
        os.remove(folder.path_of("193_S193_R2_001.fastq.gz.tmp"))
        with self.assertRaisesRegex(O.OutputFolderError, "193_S193_R2_001.fastq.gz.tmp"):
            folder.commit_job("97")
        self.assertFalse(any(self.exists(folder, name) for name in TILE_97))

    def test_a_job_counts_as_complete_only_once_it_has_a_record(self):
        folder = self.folder_under_test
        self.write_tmp_files(folder, "97")
        folder.commit_job("97")
        self.assertEqual(folder.completed_jobs(), set())            # renamed, but the record is the commit point
        self.write_job_record(folder, "97")
        self.assertEqual(folder.completed_jobs(), {"97"})


class FailureTests(FolderTestCase):
    def setUp(self):
        super().setUp()
        self.folder_under_test = self.folder()
        self.folder_under_test.prepare(O.RunMode.NEW)
        self.finished_job(self.folder_under_test, "97")
        self.write_tmp_files(self.folder_under_test, "98")
        self.touch(os.path.join(self.out, "notes.tmp"), "someone else's file")

    def test_a_caught_failure_deletes_the_failed_jobs_tmp_files_and_keeps_the_finished_job(self):
        folder = self.folder_under_test
        folder.mark_failed("98", "R1 and R2 are out of step at pair 12")
        self.assertFalse(any(self.exists(folder, name + ".tmp") for name in TILE_98))
        self.assertTrue(all(self.exists(folder, name) for name in TILE_97))

    def test_an_unrelated_tmp_file_survives(self):
        self.folder_under_test.mark_failed("98", "boom")
        self.assertTrue(self.exists(self.folder_under_test, "notes.tmp"))

    def test_keep_tmp_leaves_the_tmp_files_for_debugging(self):
        folder = self.folder_under_test
        folder.mark_failed("98", "boom", keep_tmp=True)
        self.assertTrue(all(self.exists(folder, name + ".tmp") for name in TILE_98))

    def test_final_named_files_of_the_failed_job_go_even_with_keep_tmp(self):
        """A job can fail part way through its renames. A final name must only ever mean a complete job."""
        folder = self.folder_under_test
        os.replace(folder.path_of("98_S98_R1_001.fastq.gz.tmp"), folder.path_of("98_S98_R1_001.fastq.gz"))
        folder.mark_failed("98", "died during the renames", keep_tmp=True)
        self.assertFalse(self.exists(folder, "98_S98_R1_001.fastq.gz"))

    def test_the_run_record_and_the_marker_say_which_job_failed_and_why(self):
        folder = self.folder_under_test
        folder.mark_failed("98", "R1 and R2 are out of step at pair 12")
        with open(folder.run_record) as handle:
            record = json.load(handle)
        self.assertEqual((record["status"], record["failed_job"], record["jobs_complete"], record["jobs_total"]), ("failed", "98", 1, 2))
        with open(folder.marker) as handle:
            text = handle.read()
        for expected in ("Job 98 failed: R1 and R2 are out of step at pair 12", "1 of 2 jobs finished"):
            self.assertIn(expected, text)

    def test_a_failed_folder_is_incomplete_and_has_no_links_or_manifest(self):
        folder = self.folder_under_test
        folder.mark_failed("98", "boom")
        self.assertEqual(folder.state, O.FolderState.INCOMPLETE)
        self.assertFalse(any(self.exists(folder, name) for name in LINKS + [O.MANIFEST_NAME]))

    def test_a_hard_kill_leaves_tmp_files_and_a_running_status_and_the_next_run_refuses(self):
        """Nothing runs after a kill: build a fresh OutputFolder, as the next command would."""
        again = self.folder()
        self.assertEqual(again.state, O.FolderState.INCOMPLETE)
        with open(again.run_record) as handle:
            self.assertEqual(json.load(handle)["status"], "running")
        self.assertTrue(all(self.exists(again, name + ".tmp") for name in TILE_98))
        with self.assertRaisesRegex(O.OutputFolderError, "1 of 2 jobs complete"):
            again.prepare(O.RunMode.NEW)


class ResumeTests(FolderTestCase):
    def setUp(self):
        super().setUp()
        first = self.folder()
        first.prepare(O.RunMode.NEW)
        self.finished_job(first, "97")
        self.write_tmp_files(first, "98")
        first.mark_failed("98", "boom", keep_tmp=True)
        self.second = self.folder()

    def test_resume_clears_only_the_unfinished_jobs_files(self):
        self.second.prepare(O.RunMode.RESUME)
        self.assertTrue(all(self.exists(self.second, name) for name in TILE_97))
        self.assertFalse(any(self.exists(self.second, name + ".tmp") for name in TILE_98))

    def test_resume_also_clears_final_named_files_of_a_job_with_no_record(self):
        """Renamed, then killed before the record was written: the job has to be redone."""
        os.replace(self.second.path_of("98_S98_R1_001.fastq.gz.tmp"), self.second.path_of("98_S98_R1_001.fastq.gz"))
        self.second.prepare(O.RunMode.RESUME)
        self.assertFalse(self.exists(self.second, "98_S98_R1_001.fastq.gz"))

    def test_resume_keeps_the_record_of_the_finished_job(self):
        self.second.prepare(O.RunMode.RESUME)
        self.assertEqual(self.second.completed_jobs(), {"97"})

    def test_resume_sets_the_status_back_to_running_and_counts_the_resumes(self):
        self.second.prepare(O.RunMode.RESUME)
        with open(self.second.run_record) as handle:
            record = json.load(handle)
        self.assertEqual((record["status"], record["resumes"], record["failed_job"]), ("running", 1, None))

    def test_resume_writes_the_marker_again_and_an_unrelated_file_survives(self):
        self.touch(os.path.join(self.out, "notes.tmp"))
        self.second.prepare(O.RunMode.RESUME)
        self.assertTrue(self.exists(self.second, O.MARKER_NAME))
        self.assertTrue(self.exists(self.second, "notes.tmp"))


class ForceTests(FolderTestCase):
    def test_force_removes_what_the_run_made_and_starts_again(self):
        folder = self.complete_run()
        self.folder().prepare(O.RunMode.FORCE)
        after = self.names_in(folder)
        self.assertEqual(after, sorted([O.LOGS_NAME, O.MARKER_NAME]))

    def test_force_removes_the_files_of_the_recorded_plan_even_if_the_sheet_changed(self):
        self.complete_run()
        smaller = SheetReader().parse(swap(make_sheet(), "97,193,,,\n", ""))      # 97 and 193 are no longer a group
        folder = self.folder(sheet=smaller)
        folder.prepare(O.RunMode.FORCE)
        self.assertFalse(any(self.exists(folder, name) for name in TILE_97))
        self.assertEqual(folder.plan.jobs, ("98",))

    def test_force_can_clear_a_tile_file_whose_name_is_a_link_in_the_new_plan(self):
        """A sample that was split last time is unshared now: its name was a tile file, and is a link in the new plan."""
        self.complete_run()
        smaller = SheetReader().parse(swap(make_sheet(), "97,193,,,\n", ""))
        folder = self.folder(sheet=smaller)
        folder.prepare(O.RunMode.FORCE)
        self.assertFalse(self.exists(folder, "97_S97_R1_001.fastq.gz"))
        folder.link_raw()
        self.assertTrue(os.path.islink(folder.path_of("97_S97_R1_001.fastq.gz")))

    def test_force_leaves_files_that_are_not_in_any_plan(self):
        self.complete_run()
        self.touch(os.path.join(self.out, "notes.tmp"))
        self.touch(os.path.join(self.out, "keep_me.txt"))
        folder = self.folder()
        folder.prepare(O.RunMode.FORCE)
        self.assertTrue(self.exists(folder, "notes.tmp") and self.exists(folder, "keep_me.txt"))

    def test_force_refuses_and_deletes_nothing_if_a_link_name_is_held_by_a_regular_file(self):
        folder = self.complete_run()
        os.remove(folder.path_of(LINKS[0]))
        self.touch(folder.path_of(LINKS[0]), "a real file")
        with self.assertRaisesRegex(O.OutputFolderError, "regular files"):
            self.folder().prepare(O.RunMode.FORCE)
        self.assertTrue(all(self.exists(folder, name) for name in TILE_97 + TILE_98))

    def test_force_deletes_links_but_never_what_they_point_to(self):
        self.complete_run()
        self.folder().prepare(O.RunMode.FORCE)
        self.assertEqual(sorted(os.listdir(self.raw)), sorted(RAW_NAMES))


class FinishTests(FolderTestCase):
    def test_finish_needs_a_record_for_every_job(self):
        folder = self.folder()
        folder.prepare(O.RunMode.NEW)
        self.finished_job(folder, "97")
        with self.assertRaisesRegex(O.OutputFolderError, "no record for job.*98"):
            folder.finish()

    def test_finish_links_the_raw_files_with_absolute_targets(self):
        folder = self.complete_run()
        for name in LINKS:
            path = folder.path_of(name)
            self.assertTrue(os.path.islink(path))
            self.assertEqual(os.readlink(path), os.path.join(self.raw, name))

    def test_finish_does_not_link_the_pooled_files_of_the_split_samples(self):
        folder = self.complete_run()
        self.assertFalse(os.path.islink(folder.path_of("97_S97_R1_001.fastq.gz")))     # that name is the tile file

    def test_finish_marks_the_run_complete_and_removes_the_marker_last(self):
        folder = self.complete_run()
        self.assertFalse(self.exists(folder, O.MARKER_NAME))
        with open(folder.run_record) as handle:
            record = json.load(handle)
        self.assertEqual((record["status"], record["jobs_complete"], record["jobs_total"]), ("complete", 2, 2))

    def test_linking_twice_is_fine_but_a_different_file_in_the_way_is_an_error(self):
        folder = self.complete_run()
        folder.link_raw()
        os.remove(folder.path_of(LINKS[0]))
        self.touch(folder.path_of(LINKS[0]))
        with self.assertRaisesRegex(O.OutputFolderError, "cannot link"):
            folder.link_raw()

    def test_a_tile_file_is_a_real_file_and_a_link_is_a_link(self):
        folder = self.complete_run()
        self.assertTrue(os.path.isfile(folder.path_of(TILE_97[0])) and not os.path.islink(folder.path_of(TILE_97[0])))

    def test_the_output_folder_holds_exactly_the_planned_files_and_the_reports(self):
        folder = self.complete_run()
        expected = sorted(TILE_97 + TILE_98 + LINKS + [O.LOGS_NAME, O.MANIFEST_NAME])
        self.assertEqual(self.names_in(folder), expected)


class ReportTests(FolderTestCase):
    def test_a_report_is_written_whole_with_no_tmp_file_left_behind(self):
        folder = self.folder()
        folder.prepare(O.RunMode.NEW)
        folder.write_report(O.MANIFEST_NAME, "a\tb\n")
        with open(folder.manifest) as handle:
            self.assertEqual(handle.read(), "a\tb\n")
        self.assertNotIn(O.MANIFEST_NAME + ".tmp", self.names_in(folder))

    def test_update_run_record_adds_fields_and_keeps_the_others(self):
        folder = self.folder()
        folder.prepare(O.RunMode.NEW)
        folder.update_run_record(version="0.0.1")
        with open(folder.run_record) as handle:
            record = json.load(handle)
        self.assertEqual((record["status"], record["version"]), ("running", "0.0.1"))


if __name__ == "__main__":
    unittest.main()