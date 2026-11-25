# These tests were generated with Claude (Anthropic) and then reviewed, edited and run by the maintainers.
"""
Tests for the FASTQ reader and the output writer (slice B of ticket 1.4): pairs come out in order and exactly as they
were read, every kind of bad input is named with its file and record, the output is valid gzip even after a failure,
and memory does not grow with the number of reads.

Run with:  python -m pytest tests/test_demux_fastq.py -v     (from the repository root)
Without pytest: python -m unittest discover -s tests -p "test_demux_*.py"
"""
import gzip
import os
import tempfile
import tracemalloc
import unittest

from TileSeqDemux import fastq as F


def record(read_id, sequence=b"ACGTACGT", quality=None, comment=b" 1:N:0:ACGT"):
    """One FASTQ record as bytes."""
    quality = quality if quality is not None else b"I" * len(sequence)
    return b"@" + read_id + comment + b"\n" + sequence + b"\n+\n" + quality + b"\n"


class FastqTestCase(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.root = self.directory.name

    def tearDown(self):
        self.directory.cleanup()

    def path(self, name):
        return os.path.join(self.root, name)

    def write_gz(self, name, content):
        """A gzip file holding `content` (bytes, or a list of records)."""
        data = content if isinstance(content, bytes) else b"".join(content)
        with gzip.open(self.path(name), "wb") as handle:
            handle.write(data)
        return self.path(name)

    def pair_files(self, r1_records, r2_records):
        return self.write_gz("r1.fastq.gz", r1_records), self.write_gz("r2.fastq.gz", r2_records)

    def read_all(self, r1_path, r2_path):
        with F.FastqPairReader(r1_path, r2_path) as reader:
            return list(reader.pairs())

    def gunzip(self, path):
        with gzip.open(path, "rb") as handle:
            return handle.read()


class ReaderTests(FastqTestCase):
    def test_pairs_come_out_in_order_and_exactly_as_read(self):
        r1 = [record(b"read1"), record(b"read2", b"TTTT"), record(b"read3", b"GGGGGG")]
        r2 = [record(b"read1", b"AAAA"), record(b"read2"), record(b"read3", b"CC")]
        pairs = self.read_all(*self.pair_files(r1, r2))
        self.assertEqual([pair.read_id for pair in pairs], [b"read1", b"read2", b"read3"])
        self.assertEqual([pair.r1 for pair in pairs], r1)
        self.assertEqual([pair.r2 for pair in pairs], r2)

    def test_the_read_id_is_the_first_token_after_the_at_sign(self):
        pair = self.read_all(*self.pair_files([record(b"A00:123:HX:1:1101:1000:2000", comment=b" 1:N:0:ACGT+TTTT")],
                                              [record(b"A00:123:HX:1:1101:1000:2000", comment=b" 2:N:0:ACGT+TTTT")]))[0]
        self.assertEqual(pair.read_id, b"A00:123:HX:1:1101:1000:2000")

    def test_a_header_with_no_comment_works(self):
        pair = self.read_all(*self.pair_files([record(b"read1", comment=b"")], [record(b"read1", comment=b"")]))[0]
        self.assertEqual(pair.read_id, b"read1")

    def test_empty_files_have_no_pairs(self):
        self.assertEqual(self.read_all(*self.pair_files([], [])), [])

    def test_the_reader_is_lazy_and_reads_one_pair_at_a_time(self):
        r1_path, r2_path = self.pair_files([record(b"a"), record(b"b")], [record(b"a"), record(b"b")])
        with F.FastqPairReader(r1_path, r2_path) as reader:
            stream = reader.pairs()
            self.assertEqual(next(stream).read_id, b"a")
            self.assertEqual(next(stream).read_id, b"b")
            with self.assertRaises(StopIteration):
                next(stream)

    def test_a_last_line_without_a_line_end_is_accepted_and_kept_as_it_is(self):
        last = record(b"read2").rstrip(b"\n")
        pairs = self.read_all(*self.pair_files([record(b"read1"), last], [record(b"read1"), last]))
        self.assertEqual(pairs[1].r1, last)

    def test_windows_line_ends_are_kept_and_accepted(self):
        crlf = b"@read1 1:N:0:X\r\nACGT\r\n+\r\nIIII\r\n"
        pair = self.read_all(*self.pair_files(crlf, crlf))[0]
        self.assertEqual((pair.read_id, pair.r1), (b"read1", crlf))

    def test_a_record_that_is_cut_off_is_named(self):
        r1_path, r2_path = self.pair_files(record(b"read1") + b"@read2\nACGT\n", [record(b"read1"), record(b"read2")])
        with self.assertRaisesRegex(F.FastqFormatError, r"r1\.fastq\.gz, record 2: the record is cut off"):
            self.read_all(r1_path, r2_path)

    def test_a_first_line_that_does_not_start_with_at(self):
        r1_path, r2_path = self.pair_files(record(b"read1") + record(b"read2").replace(b"@", b">", 1), [record(b"read1"), record(b"read2")])
        with self.assertRaisesRegex(F.FastqFormatError, r"r1\.fastq\.gz, record 2: the first line does not start with @"):
            self.read_all(r1_path, r2_path)

    def test_a_third_line_that_does_not_start_with_plus(self):
        bad = b"@read2\nACGT\nIIII\nIIII\n"
        r1_path, r2_path = self.pair_files([record(b"read1"), bad], [record(b"read1"), record(b"read2")])
        with self.assertRaisesRegex(F.FastqFormatError, r"record 2: the third line does not start with \+"):
            self.read_all(r1_path, r2_path)

    def test_a_plus_line_that_repeats_the_header_is_fine_and_kept(self):
        repeated = b"@read1 1:N:0:X\nACGT\n+read1 1:N:0:X\nIIII\n"
        pair = self.read_all(*self.pair_files(repeated, repeated))[0]
        self.assertEqual(pair.r1, repeated)

    def test_a_sequence_and_a_quality_of_different_length(self):
        r1_path, r2_path = self.pair_files([record(b"read1", b"ACGT", quality=b"III")], [record(b"read1")])
        with self.assertRaisesRegex(F.FastqFormatError, r"record 1: the sequence and the quality line have different lengths"):
            self.read_all(r1_path, r2_path)

    def test_a_header_with_no_read_id(self):
        r1_path, r2_path = self.pair_files(b"@\nACGT\n+\nIIII\n", [record(b"read1")])
        with self.assertRaisesRegex(F.FastqFormatError, r"record 1: the first line has no read ID"):
            self.read_all(r1_path, r2_path)

    def test_r1_and_r2_with_different_ids_are_out_of_step(self):
        r1_path, r2_path = self.pair_files([record(b"read1"), record(b"read2")], [record(b"read1"), record(b"read9")])
        with self.assertRaisesRegex(F.FastqFormatError, r"out of step at pair 2.*b'read2' and b'read9'"):
            self.read_all(r1_path, r2_path)

    def test_r1_with_more_records_than_r2(self):
        r1_path, r2_path = self.pair_files([record(b"a"), record(b"b")], [record(b"a")])
        with self.assertRaisesRegex(F.FastqFormatError, r"r2\.fastq\.gz has fewer records than its mate: it ends after record 1"):
            self.read_all(r1_path, r2_path)

    def test_r2_with_more_records_than_r1(self):
        r1_path, r2_path = self.pair_files([record(b"a")], [record(b"a"), record(b"b")])
        with self.assertRaisesRegex(F.FastqFormatError, r"r1\.fastq\.gz has fewer records than its mate"):
            self.read_all(r1_path, r2_path)

    def test_pairs_before_the_problem_are_still_delivered(self):
        r1_path, r2_path = self.pair_files([record(b"a"), record(b"b")], [record(b"a"), record(b"x")])
        delivered = []
        with F.FastqPairReader(r1_path, r2_path) as reader:
            with self.assertRaises(F.FastqFormatError):
                for pair in reader.pairs():
                    delivered.append(pair.read_id)
        self.assertEqual(delivered, [b"a"])

    def test_a_blank_line_at_the_end_is_an_error_not_a_silent_skip(self):
        r1_path, r2_path = self.pair_files(record(b"a") + b"\n", record(b"a"))
        with self.assertRaisesRegex(F.FastqFormatError, r"record 2: the record is cut off"):
            self.read_all(r1_path, r2_path)

    def test_a_gzip_file_that_is_cut_off_is_named(self):
        full = self.write_gz("full.fastq.gz", [record(b"read%d" % i) for i in range(2000)])
        with open(full, "rb") as handle:
            data = handle.read()
        with open(self.path("r1.fastq.gz"), "wb") as handle:
            handle.write(data[: int(len(data) * 0.6)])
        r2_path = self.write_gz("r2.fastq.gz", [record(b"read%d" % i) for i in range(2000)])
        with self.assertRaisesRegex(F.FastqFormatError, r"r1\.fastq\.gz"):
            self.read_all(self.path("r1.fastq.gz"), r2_path)

    def test_a_file_that_is_not_gzip_is_named(self):
        with open(self.path("r1.fastq.gz"), "wb") as handle:
            handle.write(b"@read1\nACGT\n+\nIIII\n")
        r2_path = self.write_gz("r2.fastq.gz", [record(b"read1")])
        with self.assertRaisesRegex(F.FastqFormatError, r"r1\.fastq\.gz.*cut off or damaged"):
            self.read_all(self.path("r1.fastq.gz"), r2_path)

    def test_a_missing_file_is_the_usual_os_error(self):
        r2_path = self.write_gz("r2.fastq.gz", [record(b"read1")])
        with self.assertRaises(FileNotFoundError):
            self.read_all(self.path("nope.fastq.gz"), r2_path)

    def test_a_failed_open_closes_the_file_that_was_already_open(self):
        reader = F.FastqPairReader(self.write_gz("r1.fastq.gz", [record(b"a")]), self.path("missing.fastq.gz"))
        with self.assertRaises(FileNotFoundError):
            reader.__enter__()
        self.assertTrue(reader._r1._handle.closed)

    def test_pairs_without_a_with_block_says_what_to_do(self):
        r1_path, r2_path = self.pair_files([record(b"a")], [record(b"a")])
        with self.assertRaisesRegex(ValueError, "with block"):
            list(F.FastqPairReader(r1_path, r2_path).pairs())

    def test_the_files_are_closed_when_the_reader_is_left(self):
        r1_path, r2_path = self.pair_files([record(b"a")], [record(b"a")])
        with F.FastqPairReader(r1_path, r2_path) as reader:
            pass
        with self.assertRaises(ValueError):
            list(reader.pairs())


class WriterTests(FastqTestCase):
    def outputs_for(self, *keys):
        return {key: (self.path(f"{key}_R1.tmp"), self.path(f"{key}_R2.tmp")) for key in keys}

    def pair(self, read_id=b"read1", sequence=b"ACGT"):
        return F.ReadPair(read_id, record(read_id, sequence), record(read_id, sequence[::-1]))

    def test_each_key_gets_an_r1_and_an_r2_file_holding_its_pairs(self):
        paths = self.outputs_for("97", "193")
        first, second, third = self.pair(b"a"), self.pair(b"b"), self.pair(b"c")
        with F.TileOutputs(paths) as outputs:
            outputs.write("97", first)
            outputs.write("193", second)
            outputs.write("97", third)
        self.assertEqual(self.gunzip(paths["97"][0]), first.r1 + third.r1)
        self.assertEqual(self.gunzip(paths["97"][1]), first.r2 + third.r2)
        self.assertEqual(self.gunzip(paths["193"][0]), second.r1)

    def test_the_counts_say_how_many_pairs_each_key_received(self):
        with F.TileOutputs(self.outputs_for("97", "193")) as outputs:
            outputs.write("97", self.pair(b"a"))
            outputs.write("97", self.pair(b"b"))
        self.assertEqual(outputs.counts, {"97": 2, "193": 0})

    def test_a_key_with_no_pairs_still_gets_valid_empty_files(self):
        paths = self.outputs_for("97")
        with F.TileOutputs(paths):
            pass
        self.assertEqual(self.gunzip(paths["97"][0]), b"")

    def test_one_pair_can_go_to_two_keys(self):
        paths = self.outputs_for("97", "193")
        both = self.pair()
        with F.TileOutputs(paths) as outputs:
            outputs.write("97", both)
            outputs.write("193", both)
        self.assertEqual(self.gunzip(paths["97"][0]), self.gunzip(paths["193"][0]))

    def test_writing_to_a_key_that_was_not_given_is_an_error(self):
        with F.TileOutputs(self.outputs_for("97")) as outputs:
            with self.assertRaises(KeyError):
                outputs.write("999", self.pair())

    def test_writing_before_the_outputs_are_opened_says_what_to_do(self):
        outputs = F.TileOutputs(self.outputs_for("97"))
        with self.assertRaisesRegex(ValueError, "with block"):
            outputs.write("97", self.pair())

    def test_the_files_are_valid_gzip_even_if_the_job_failed(self):
        paths = self.outputs_for("97")
        with self.assertRaises(RuntimeError):
            with F.TileOutputs(paths) as outputs:
                outputs.write("97", self.pair(b"a"))
                raise RuntimeError("the job failed")
        self.assertEqual(self.gunzip(paths["97"][0]), self.pair(b"a").r1)

    def test_paths_lists_every_file_it_writes(self):
        paths = self.outputs_for("97", "193")
        self.assertEqual(F.TileOutputs(paths).paths(), [paths["97"][0], paths["97"][1], paths["193"][0], paths["193"][1]])

    def test_discard_deletes_the_files_it_was_given_and_nothing_else(self):
        paths = self.outputs_for("97")
        for name in ("notes.tmp", "other_R1.tmp", "97_R1.tmp.keep"):
            open(self.path(name), "w").close()
        with F.TileOutputs(paths) as outputs:
            outputs.write("97", self.pair())
            outputs.discard()
        self.assertEqual(sorted(os.listdir(self.root)), ["97_R1.tmp.keep", "notes.tmp", "other_R1.tmp"])

    def test_discard_after_the_files_are_gone_is_fine(self):
        outputs = F.TileOutputs(self.outputs_for("97"))
        outputs.discard()

    def test_discard_never_deletes_a_link(self):
        target = self.path("precious.txt")
        open(target, "w").close()
        link = self.path("97_R1.tmp")
        os.symlink(target, link)
        F.TileOutputs({"97": (link, self.path("97_R2.tmp"))}).discard()
        self.assertTrue(os.path.exists(target) and os.path.islink(link))

    def test_a_file_that_cannot_be_created_closes_the_ones_already_open(self):
        paths = {"97": (self.path("97_R1.tmp"), self.path("97_R2.tmp")), "193": (self.path("no_such_folder/193_R1.tmp"), self.path("193_R2.tmp"))}
        with self.assertRaises(FileNotFoundError):
            with F.TileOutputs(paths):
                pass
        self.assertEqual(self.gunzip(paths["97"][0]), b"")           # opened first, and closed properly

    def test_reading_and_writing_again_gives_the_same_bytes(self):
        r1 = [record(b"read%d" % number, b"ACGT" * (number % 5 + 1)) for number in range(50)]
        r2 = [record(b"read%d" % number, b"TTGA" * (number % 3 + 1)) for number in range(50)]
        paths = self.outputs_for("97")
        with F.FastqPairReader(*self.pair_files(r1, r2)) as reader, F.TileOutputs(paths) as outputs:
            for pair in reader.pairs():
                outputs.write("97", pair)
        self.assertEqual(self.gunzip(paths["97"][0]), b"".join(r1))
        self.assertEqual(self.gunzip(paths["97"][1]), b"".join(r2))


class MemoryTests(FastqTestCase):
    """The pass must not hold the reads: ten times the pairs must not need ten times the memory."""

    def peak_bytes(self, pairs):
        records = [record(b"read%d" % number, b"ACGT" * 25) for number in range(pairs)]
        r1_path, r2_path = self.pair_files(records, records)
        paths = {"97": (self.path("o_R1.tmp"), self.path("o_R2.tmp"))}
        tracemalloc.start()
        try:
            with F.FastqPairReader(r1_path, r2_path) as reader, F.TileOutputs(paths) as outputs:
                for pair in reader.pairs():
                    outputs.write("97", pair)
            return tracemalloc.get_traced_memory()[1]
        finally:
            tracemalloc.stop()

    def test_ten_times_the_pairs_needs_about_the_same_memory(self):
        small, large = self.peak_bytes(4000), self.peak_bytes(40000)
        self.assertLess(large, small * 2 + 100000, f"peak memory grew from {small} to {large} bytes")


if __name__ == "__main__":
    unittest.main()