"""FASTQ streams: reading R1 and R2 together, and writing the output files of a job (slice B of ticket 1.4).

    ReadPair          one read pair: its ID, and its R1 and R2 records exactly as they were read
    FastqPairReader   reads an R1 file and an R2 file in step, one pair at a time, and checks what it reads
    TileOutputs       the gzip files a job writes, an R1 and an R2 file for each key (an output sample)

A record is four lines. The reader keeps them as bytes and never decodes or rebuilds them, so what is written is
byte for byte what was read. Only one pair is in memory at a time. As in the prototype, the read ID is the first
token of the first line after the "@", and R1 and R2 must give the same one (names ending /1 and /2 are not handled).

The reader raises FastqFormatError, naming the file and the record, for a record that is cut off, a first line that
does not start with "@", a third line that does not start with "+", a sequence and a quality line of different length,
R1 and R2 with different IDs or different numbers of records, and a gzip file that is cut off or damaged.

Standard library only; nothing here imports another module of the package.
"""
import gzip
import itertools
import os
import zlib
from contextlib import ExitStack
from typing import NamedTuple


class FastqFormatError(Exception):
    """A FASTQ file is not what a read pair needs. The message names the file and the record."""


class ReadPair(NamedTuple):
    read_id: bytes    # the first token of the header, without the "@"
    r1: bytes         # the four lines of the R1 record, line ends included
    r2: bytes         # the four lines of the R2 record


class _MateStream:
    """One gzipped FASTQ file, read a record at a time."""

    def __init__(self, path):
        self.path = path
        self._handle = None
        self._readline = None

    def open(self):
        self._handle = gzip.open(self.path, "rb")
        self._readline = self._handle.readline

    def close(self):
        if self._handle is not None:
            self._handle.close()

    def read(self, number):
        """The next record as (read ID, the four lines joined), or None at the end of the file."""
        try:
            header = self._readline()
            if not header:
                return None
            sequence, plus, quality = self._readline(), self._readline(), self._readline()
        except (EOFError, OSError, zlib.error) as error:
            raise FastqFormatError(f"{self.path}: cannot be read at record {number} ({error}); "
                                   f"the file is cut off or damaged") from error
        problem = self._problem(header, sequence, plus, quality)
        if problem:
            raise FastqFormatError(f"{self.path}, record {number}: {problem}")
        return self._read_id(header, number), header + sequence + plus + quality

    @staticmethod
    def _problem(header, sequence, plus, quality):
        """Why the four lines are not a FASTQ record, or None when they are. The first check that fails is reported."""
        checks = ((not quality, "the record is cut off"),
                  (header[:1] != b"@", "the first line does not start with @"),
                  (plus[:1] != b"+", "the third line does not start with +"),
                  (len(sequence) != len(quality) and len(sequence.rstrip(b"\r\n")) != len(quality.rstrip(b"\r\n")),
                   "the sequence and the quality line have different lengths"))      # a last line may lack its line end
        for broken, message in checks:
            if broken:
                return message
        return None

    def _read_id(self, header, number):
        try:
            return header[1:].split(None, 1)[0]
        except IndexError:
            raise FastqFormatError(f"{self.path}, record {number}: the first line has no read ID") from None


class FastqPairReader:
    """Reads an R1 file and an R2 file together. Use it as a context manager, which opens and closes both files."""

    def __init__(self, r1_path, r2_path):
        self._r1 = _MateStream(r1_path)
        self._r2 = _MateStream(r2_path)

    def __enter__(self):
        try:
            self._r1.open()
            self._r2.open()
        except BaseException:
            self.__exit__()
            raise
        return self

    def __exit__(self, *error):
        self._r1.close()
        self._r2.close()

    def pairs(self):
        """Yield a ReadPair for every pair, in file order. Pairs are counted from 1 in the messages."""
        for number in itertools.count(1):
            first, second = self._r1.read(number), self._r2.read(number)
            if first is None and second is None:
                return
            yield self._pair(first, second, number)

    def _pair(self, first, second, number):
        if first is None or second is None:
            shorter = self._r1 if first is None else self._r2
            raise FastqFormatError(f"{shorter.path} has fewer records than its mate: it ends after record {number - 1}")
        if first[0] != second[0]:
            raise FastqFormatError(f"R1 and R2 are out of step at pair {number}: "
                                   f"the read IDs are {first[0]!r} and {second[0]!r}")
        return ReadPair(first[0], first[1], second[1])


class TileOutputs:
    """The gzip files a job writes: for each key an R1 file and an R2 file, at the paths it is given. The caller
    chooses the paths, which output.py makes .tmp names, so this class does not know what a finished name is.

    Use it as a context manager: opening creates every file, leaving closes them all, also after an error, so a file is
    always a valid gzip stream even if the job failed. `counts` says how many pairs each key has received."""

    def __init__(self, paths, compresslevel=6):
        self._paths = {key: tuple(pair) for key, pair in paths.items()}
        self._compresslevel = compresslevel
        self._handles = {}
        self._stack = ExitStack()
        self.counts = {key: 0 for key in self._paths}

    def __enter__(self):
        try:
            for key, (r1_path, r2_path) in self._paths.items():
                self._handles[key] = tuple(self._stack.enter_context(gzip.open(path, "wb", compresslevel=self._compresslevel))
                                           for path in (r1_path, r2_path))
        except BaseException:
            self._stack.close()
            raise
        return self

    def __exit__(self, *error):
        self._stack.close()

    def write(self, key, pair):
        """Write a pair's R1 record to the key's R1 file and its R2 record to the R2 file."""
        r1_file, r2_file = self._handles[key]
        r1_file.write(pair.r1)
        r2_file.write(pair.r2)
        self.counts[key] += 1

    def paths(self):
        """Every file this object writes, in key order."""
        return list(itertools.chain.from_iterable(self._paths.values()))

    def discard(self):
        """Close the files and delete exactly the files this object was given, and nothing else."""
        self._stack.close()
        for path in self.paths():
            if os.path.isfile(path) and not os.path.islink(path):
                os.remove(path)