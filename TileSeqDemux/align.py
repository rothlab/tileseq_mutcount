"""Reading what bowtie2 writes: one alignment per read, joined to the reads as they stream past (slice C of ticket 1.4).

    Alignment        where a mate aligned: its first position, the length of its sequence and the template it covers
    AlignmentStream  the SAM lines of one bowtie2 run, which answers "how did this read align?" read by read

bowtie2 is run with --reorder, so its records come in the order of the input reads, and with --no-unal, so reads that
did not align have no record. A read has therefore aligned exactly when the next record carries its name, and the join
needs no read IDs held in memory. The caller checks at the end that the stream was used up (`leftover`) and that its
count agrees with bowtie2's own summary.

Standard library only; nothing here imports another module of the package.
"""
import re
from typing import NamedTuple

REF_OPS = b"MDN=X"                                # CIGAR operations that use up template (reference) bases
CIGAR = re.compile(rb"(\d+)([MIDNSHP=X])")


class AlignmentFormatError(Exception):
    """A line of bowtie2's output cannot be read as an alignment."""


class Alignment(NamedTuple):
    position: int     # POS: the 1-based template position of the first aligned base
    length: int       # the length of SEQ, as in tileseq_mut's test of where a read ends
    span: int         # the template bases the alignment covers, from its CIGAR


def parse_sam_line(line):
    """The Alignment described by one SAM line (bytes)."""
    fields = line.split(b"\t", 10)
    try:
        span = sum(int(count) for count, operation in CIGAR.findall(fields[5]) if operation in REF_OPS)
        return Alignment(int(fields[3]), len(fields[9]), span)
    except (IndexError, ValueError):
        raise AlignmentFormatError(f"cannot read this SAM line: {line[:80]!r}") from None


class AlignmentStream:
    """The SAM lines of one bowtie2 run, taken one read at a time. `lines` is anything that yields lines as bytes,
    such as the stdout of the bowtie2 process or a list in a test."""

    def __init__(self, lines):
        self._lines = iter(lines)
        self._line = None
        self._name = None
        self.records_read = 0               # SAM lines seen so far, compared with bowtie2's summary at the end
        self._advance()

    def _advance(self):
        self._line = next(self._lines, None)
        if self._line is None:
            self._name = None
            return
        self.records_read += 1
        tab = self._line.find(b"\t")
        if tab < 0:
            raise AlignmentFormatError(f"cannot read this SAM line: {self._line[:80]!r}")
        self._name = self._line[:tab]

    def take(self, read_id):
        """The Alignment of the read called `read_id`, or None if it did not align. Call it for every read in order."""
        if self._name != read_id:
            return None
        alignment = parse_sam_line(self._line)
        while self._name == read_id:        # any further lines for the same read are skipped
            self._advance()
        return alignment

    @property
    def leftover(self):
        """The first line that was never matched to a read, or None when the stream was used up. A line left over means
        the output was not in input order or the read names do not match."""
        return self._line