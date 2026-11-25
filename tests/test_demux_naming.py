# These tests were generated with Claude (Anthropic) and then reviewed, edited and run by the maintainers.
"""
Tests for the FASTQ naming contract that tileseq_demux shares with tileseq_mut (slice A of ticket 1.4).

Run with:  python -m pytest tests/test_demux_naming.py -v     (from the repository root)
Without pytest: python -m unittest discover -s tests -p "test_demux_*.py"
"""
import unittest

from TileSeqDemux import naming


class SampleIdTests(unittest.TestCase):
    def test_the_id_is_the_first_underscore_token(self):
        self.assertEqual(naming.sample_id("97_S97_R1_001.fastq.gz"), "97")
        self.assertEqual(naming.sample_id("Undetermined_S0_R1_001.fastq.gz"), "Undetermined")

    def test_a_name_without_an_underscore_has_no_id(self):
        self.assertIsNone(naming.sample_id("README"))
        self.assertIsNone(naming.sample_id("97.fastq.gz"))


class MateFileTests(unittest.TestCase):
    def test_the_standard_illumina_name(self):
        self.assertTrue(naming.is_mate_file("97_S97_R1_001.fastq.gz", "97", 1))
        self.assertTrue(naming.is_mate_file("97_S97_R2_001.fastq.gz", "97", 2))

    def test_r1_is_not_r2(self):
        self.assertFalse(naming.is_mate_file("97_S97_R1_001.fastq.gz", "97", 2))

    def test_the_id_must_be_the_whole_first_token(self):
        self.assertFalse(naming.is_mate_file("970_S970_R1_001.fastq.gz", "97", 1))
        self.assertFalse(naming.is_mate_file("1097_S1_R1_001.fastq.gz", "97", 1))

    def test_the_file_must_end_in_fastq_gz(self):
        self.assertFalse(naming.is_mate_file("97_S97_R1_001.fastq", "97", 1))
        self.assertFalse(naming.is_mate_file("97_S97_R1_001.fastq.gz.tmp", "97", 1))

    def test_the_mate_marker_must_come_after_the_id(self):
        """tileseq_mut looks for _R1_ in the name, and the one in '97_R1_001' sits right after the ID, so it does not count."""
        self.assertFalse(naming.is_mate_file("97_R1_001.fastq.gz", "97", 1))

    def test_a_lane_number_is_allowed(self):
        self.assertTrue(naming.is_mate_file("97_S5_L001_R1_001.fastq.gz", "97", 1))


class TileFileNameTests(unittest.TestCase):
    def test_the_name_tileseq_demux_writes(self):
        self.assertEqual(naming.tile_file_name("193", 1), "193_S193_R1_001.fastq.gz")
        self.assertEqual(naming.tile_file_name("193", 2), "193_S193_R2_001.fastq.gz")

    def test_a_written_name_is_recognised_as_a_mate_file_of_its_sample(self):
        for mate in (1, 2):
            with self.subTest(mate=mate):
                self.assertTrue(naming.is_mate_file(naming.tile_file_name("193", mate), "193", mate))
                self.assertEqual(naming.sample_id(naming.tile_file_name("193", mate)), "193")


class GroupTests(unittest.TestCase):
    def test_names_are_grouped_by_id_in_the_order_given(self):
        groups = naming.group_by_sample_id(["97_S97_R1_001.fastq.gz", "98_S98_R1_001.fastq.gz", "97_S97_R2_001.fastq.gz"])
        self.assertEqual(groups, {"97": ["97_S97_R1_001.fastq.gz", "97_S97_R2_001.fastq.gz"], "98": ["98_S98_R1_001.fastq.gz"]})

    def test_names_without_an_underscore_are_left_out(self):
        self.assertEqual(naming.group_by_sample_id(["README", "notes.txt", "97_S97_R1_001.fastq.gz"]), {"97": ["97_S97_R1_001.fastq.gz"]})


if __name__ == "__main__":
    unittest.main()