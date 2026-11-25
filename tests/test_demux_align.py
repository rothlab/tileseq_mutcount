# These tests were generated with Claude (Anthropic) and then reviewed, edited and run by the maintainers.
"""
Tests for reading bowtie2's SAM output (slice C of ticket 1.4): what a line says about an alignment, and how a stream of
lines is joined to the reads as they come, with unaligned reads missing.

Run with:  python -m pytest tests/test_demux_align.py -v     (from the repository root)
Without pytest: python -m unittest discover -s tests -p "test_demux_*.py"
"""
import unittest

from TileSeqDemux import align as A


def sam(name, position, cigar="50M", length=50, flag=0):
    """A SAM line as bowtie2 writes it with --no-head, as bytes."""
    fields = [name, str(flag), "gene", str(position), "42", cigar, "*", "0", "0", "A" * length, "I" * length, "AS:i:100"]
    return ("\t".join(fields) + "\n").encode()


class ParseTests(unittest.TestCase):
    def test_position_length_and_span_of_a_plain_match(self):
        self.assertEqual(A.parse_sam_line(sam("read1", 101, "50M", 50)), A.Alignment(101, 50, 50))

    def test_span_counts_the_template_bases_the_cigar_uses(self):
        cases = {"50M": 50, "5S45M": 45, "20M2D30M": 52, "20M2I28M": 48, "30M100N20M": 150, "10=1X39=": 50, "3H47M": 47,
                 "10M5P10M": 20}
        for cigar, span in cases.items():
            with self.subTest(cigar=cigar):
                self.assertEqual(A.parse_sam_line(sam("read1", 1, cigar)).span, span)

    def test_length_is_the_length_of_the_sequence_not_of_the_alignment(self):
        self.assertEqual(A.parse_sam_line(sam("read1", 10, "5S45M", 50)).length, 50)

    def test_a_line_with_no_tabs_is_named_when_the_stream_reaches_it(self):
        with self.assertRaisesRegex(A.AlignmentFormatError, "cannot read this SAM line"):
            A.AlignmentStream([b"just text\n"])

    def test_a_line_with_too_few_fields_is_named_when_it_is_taken(self):
        stream = A.AlignmentStream([b"read1\t0\tgene\n"])
        with self.assertRaisesRegex(A.AlignmentFormatError, "cannot read this SAM line"):
            stream.take(b"read1")

    def test_a_position_that_is_not_a_number_is_named_when_it_is_taken(self):
        stream = A.AlignmentStream([sam("read1", "x")])
        with self.assertRaisesRegex(A.AlignmentFormatError, "cannot read this SAM line"):
            stream.take(b"read1")

    def test_a_cigar_of_star_gives_a_span_of_zero(self):
        """bowtie2 writes * for a read with no alignment; a CIGAR is not otherwise checked, as bowtie2 writes valid ones."""
        self.assertEqual(A.parse_sam_line(sam("read1", 5, "*")).span, 0)


class StreamTests(unittest.TestCase):
    def test_a_read_that_aligned_gets_its_alignment(self):
        stream = A.AlignmentStream([sam("read1", 101), sam("read2", 205, "40M", 40)])
        self.assertEqual(stream.take(b"read1"), A.Alignment(101, 50, 50))
        self.assertEqual(stream.take(b"read2"), A.Alignment(205, 40, 40))

    def test_a_read_that_did_not_align_gets_none_and_does_not_use_up_the_next_record(self):
        stream = A.AlignmentStream([sam("read2", 205)])
        self.assertIsNone(stream.take(b"read1"))
        self.assertEqual(stream.take(b"read2").position, 205)

    def test_reads_in_order_with_gaps_where_reads_are_missing(self):
        stream = A.AlignmentStream([sam("read1", 1), sam("read3", 3), sam("read6", 6)])
        got = {name: stream.take(name) for name in (b"read1", b"read2", b"read3", b"read4", b"read5", b"read6")}
        self.assertEqual({name: alignment is not None for name, alignment in got.items()},
                         {b"read1": True, b"read2": False, b"read3": True, b"read4": False, b"read5": False, b"read6": True})
        self.assertIsNone(stream.leftover)

    def test_a_name_is_matched_whole_so_read1_is_not_read10(self):
        stream = A.AlignmentStream([sam("read10", 10)])
        self.assertIsNone(stream.take(b"read1"))
        self.assertIsNotNone(stream.take(b"read10"))

    def test_further_lines_for_the_same_read_are_skipped(self):
        stream = A.AlignmentStream([sam("read1", 10), sam("read1", 99, flag=256), sam("read2", 20)])
        self.assertEqual(stream.take(b"read1").position, 10)
        self.assertEqual(stream.take(b"read2").position, 20)
        self.assertIsNone(stream.leftover)

    def test_the_stream_counts_the_lines_it_has_seen(self):
        stream = A.AlignmentStream([sam("read1", 10), sam("read1", 99, flag=256), sam("read2", 20)])
        stream.take(b"read1")
        stream.take(b"read2")
        self.assertEqual(stream.records_read, 3)

    def test_an_empty_stream_has_nothing_for_anyone_and_nothing_left_over(self):
        stream = A.AlignmentStream([])
        self.assertIsNone(stream.take(b"read1"))
        self.assertIsNone(stream.leftover)
        self.assertEqual(stream.records_read, 0)

    def test_a_line_that_was_never_taken_is_left_over(self):
        """A record for a read nobody asked about means the output was not in input order, or the names do not match."""
        stream = A.AlignmentStream([sam("read1", 10), sam("stranger", 20)])
        stream.take(b"read1")
        self.assertEqual(stream.leftover.split(b"\t")[0], b"stranger")

    def test_it_works_on_a_generator_and_does_not_read_ahead_more_than_one_line(self):
        pulled = []

        def lines():
            for name in ("read1", "read2", "read3"):
                pulled.append(name)
                yield sam(name, 1)
        stream = A.AlignmentStream(lines())
        self.assertEqual(pulled, ["read1"])
        stream.take(b"read1")
        self.assertEqual(pulled, ["read1", "read2"])


if __name__ == "__main__":
    unittest.main()