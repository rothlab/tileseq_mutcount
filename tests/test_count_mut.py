# These tests were generated with Claude (Anthropic) and then reviewed, edited and run by the maintainers.
"""
Tests for how count_mut.py handles a read pair that raises DuplicatedMutations.

Background
    count_mut.process_wrapper catches DuplicatedMutations, so the read pair is skipped and
    recorded in skippedreads_<sample>.csv instead of aborting the whole count. How the parser
    produces that error is tested in test_locate_mut.py.

Run with:  python -m pytest tests/test_count_mut.py -v
"""


import logging
import re

import pandas as pd
import pytest

from TileSeqMut import count_mut, locate_mut
from TileSeqMut.locate_mut import DuplicatedMutations, SameRefandAlt


### test data ###

READ_ROW = {
  "pos_start_r1" : "491", 
  "cigar_r1" : "10M3I4D10M",
  "pos_start_r2" : "520",
  "cigar_r2" : "30M"
}

HEADER_LABELS = ["#Sample", "#Raw read depth", "#Number of read pairs without mutations",
                 "#Number of read pairs did not map to gene", "#Mutation positions outside of the tile",
                 "#Number of reads outside of the tile", "#Final read-depth",
                 "#Total read pairs with mutations", "#Comment", "#Comment",
                 "#Read len for R1", "#Read len for R2"]

### test sam file data ###
def sam_line(name, cigar, md, pos=1, rname="ref"):
    length = sum(int(n) for n, op in re.findall(r"(\d+)([MIDS])", cigar) if op in "MIS")
    return f"{name}\t0\t{rname}\t{pos}\t60\t{cigar}\t*\t0\t0\t{'A' * length}\t{'I' * length}\tMD:Z:{md}\n"


def sam_pair(name, r1_cigar, r1_md, r2_cigar="20M", r2_md="20"):
    """(R1 line, R2 line) of one read pair; R2 defaults to a clean 20M alignment."""
    return sam_line(name, r1_cigar, r1_md), sam_line(name, r2_cigar, r2_md)


SKIP_1 = sam_pair("p1", "10M3I4D10M", "10^ATTC10")      # reaches the parser, then skipped
SKIP_2 = sam_pair("p2", "20M9I9D20M", "20^CCATTCAAG20")  # reaches the parser, then skipped
NO_MUT = sam_pair("p3", "20M", "20")                     # no mutation: never reaches the parser
UNMAPPED = ("p4\t4\t*\t0\t0\t*\t*\t0\t0\tAAAA\tIIII\n",  # R1 unmapped: never reaches the parser
            sam_line("p4", "20M", "20"))


### fake classes to replace the parser and the process pool ###
class ListLogger:
    """Collects log messages so tests can look at them."""
    def __init__(self):
        self.messages = []

    def _add(self, msg, *args, **kwargs):
        self.messages.append(str(msg))

    info = warning = error = debug = exception = _add


class FakeResult:
    def __init__(self, value):
        self._value = value

    def get(self):
        return self._value


class FakePool:
    """Runs every job immediately in the calling process."""
    def __init__(self, *args, **kwargs):
        pass

    def apply_async(self, func, args):
        return FakeResult(func(*args))

    def close(self):
        pass


class FakeMutParser:
    """Stands in for MutParser: reads with an insertion and a deletion raise DuplicatedMutations."""
    def __init__(self, row, *args, **kwargs):
        self._row = row

    def _main(self):
        cigar = self._row["cigar_r1"]
        if "I" in cigar and "D" in cigar:
            raise DuplicatedMutations()
        raise AssertionError("only skipped pairs should reach the parser in these tests")


def fake_parser(behaviour):
    """A MutParser stand-in whose _main() returns (or raises) whatever `behaviour()` does."""
    class FakeParser:
        def __init__(self, *args, **kwargs):
            pass

        def _main(self):
            return behaviour()
    return FakeParser


### tests ###
class TestWriteSkippedReads:
    PAIRS = [
        ("491", "10M3I4D10M", "520", "30M"),
        ("466", "11S33M9I9D63M35S", "467", "104M47S"),
    ]

    def test_writes_header_and_one_row_per_pair(self, tmp_path):
        path = tmp_path / "skippedreads_S1.csv"
        assert count_mut.write_skipped_reads(str(path), self.PAIRS) == 2
        assert path.read_text().splitlines() == [
            "r1_start,r1_cigar,r2_start,r2_cigar",
            "491,10M3I4D10M,520,30M",
            "466,11S33M9I9D63M35S,467,104M47S",
        ]

    def test_no_pairs_writes_only_the_header(self, tmp_path):
        path = tmp_path / "skippedreads_S1.csv"
        assert count_mut.write_skipped_reads(str(path), []) == 0
        assert path.read_text().splitlines() == ["r1_start,r1_cigar,r2_start,r2_cigar"]


class TestProcessWrapper:
    @staticmethod
    def run():
        return count_mut.process_wrapper(
            READ_ROW, None, None, None, 1, 100, 0.9,
            logging.getLogger("test_count_mut"),
            0.0025, 33, False, False, [])

    def test_skips_pair_on_duplicated_mutations(self, monkeypatch):
        def raise_duplicated():
            raise DuplicatedMutations()

        monkeypatch.setattr(locate_mut, "MutParser", fake_parser(raise_duplicated))
        result = self.run()
        assert len(result) == 7
        hgvs, outside_mut, all_df, r1_clusters, r2_clusters, track_df, skipped_pair = result
        assert (hgvs, outside_mut, r1_clusters, r2_clusters) == ([], [], [], [])
        assert all_df.empty and track_df.empty
        assert skipped_pair == ("491", "10M3I4D10M", "520", "30M")

    def test_passes_normal_results_through(self, monkeypatch):
        all_df, track_df = pd.DataFrame({"a": [1]}), pd.DataFrame({"b": [2]})
        monkeypatch.setattr(locate_mut, "MutParser",
                            fake_parser(lambda: ("c.1A>G", [5], all_df, ["c1"], ["c2"], track_df)))
        hgvs, outside_mut, all_out, r1_clusters, r2_clusters, track_out, skipped_pair = self.run()
        assert (hgvs, outside_mut, r1_clusters, r2_clusters) == ("c.1A>G", [5], ["c1"], ["c2"])
        assert all_out is all_df and track_out is track_df
        assert skipped_pair is None

    def test_does_not_swallow_other_errors(self, monkeypatch):
        def raise_value_error():
            raise ValueError("unrelated failure")

        monkeypatch.setattr(locate_mut, "MutParser", fake_parser(raise_value_error))
        with pytest.raises(ValueError):
            self.run()

    def test_does_not_swallow_same_ref_and_alt(self, monkeypatch):
        def raise_same_ref_alt():
            raise SameRefandAlt(pd.DataFrame({"ref": ["A"], "alt": ["A"]}))

        monkeypatch.setattr(locate_mut, "MutParser", fake_parser(raise_same_ref_alt))
        with pytest.raises(SameRefandAlt):
            self.run()


class TestMultiCore:
    """readSam.multi_core with the process pool and the parser replaced by fakes.

    Only pairs that are skipped, unmapped or mutation-free are used, so the per-position
    coverage merging at the end of the job loop is not exercised here.
    """

    @staticmethod
    def make_reader(tmp_path, pairs):
        r1, r2 = tmp_path / "S1_R1.sam", tmp_path / "S1_R2.sam"
        r1.write_text("".join(p[0] for p in pairs))
        r2.write_text("".join(p[1] for p in pairs))
        reader = count_mut.readSam.__new__(count_mut.readSam)  # skip __init__ (needs a parameter sheet)
        reader._r1, reader._r2, reader._cores = str(r1), str(r2), 1
        reader._mut_log = ListLogger()
        reader._cds_start, reader._tile_begins, reader._tile_ends, reader._min_map_len = 1, 1, 100, 10
        reader._seq = reader._cds_seq = reader._seq_lookup = None
        reader._qual, reader._mutrate, reader._base = 0.9, 0.0025, 33
        reader._posteriorQC = reader._error_override = False
        reader._output_counts_dir, reader._sample_id = str(tmp_path), "S1"
        reader._sample_counts_f = str(tmp_path / "counts_sample_S1.csv")
        reader._sample_offreads_f = str(tmp_path / "offreads_S1.csv")
        (tmp_path / "counts_sample_S1.csv").write_text("#Sample:S1\n")  # what __init__ writes
        (tmp_path / "offreads_S1.csv").write_text("r1_start, r1_len, r2_start, r2_len\n")
        track = pd.DataFrame({"pos": range(1, 101), "m_r1": 0, "m_r2": 0, "m_both": 0,
                              "passed": 0, "wt_passed": 0})
        reader._track_reads = track.set_index("pos")
        return reader

    @classmethod
    def run(cls, tmp_path, monkeypatch, pairs):
        monkeypatch.setattr(count_mut.mp, "Pool", FakePool)
        monkeypatch.setattr(locate_mut, "MutParser", FakeMutParser)
        reader = cls.make_reader(tmp_path, pairs)
        reader.multi_core([])
        return reader

    def test_records_skipped_pairs_and_still_writes_outputs(self, tmp_path, monkeypatch):
        reader = self.run(tmp_path, monkeypatch, [SKIP_1, SKIP_2, NO_MUT, UNMAPPED])
        assert (tmp_path / "skippedreads_S1.csv").read_text().splitlines() == [
            "r1_start,r1_cigar,r2_start,r2_cigar",
            "1,10M3I4D10M,1,20M",
            "1,20M9I9D20M,1,20M",
        ]
        assert (tmp_path / "coverage_S1.csv").exists()
        counts = (tmp_path / "counts_sample_S1.csv").read_text().splitlines()
        assert counts[-1] == "HGVS,count"
        assert "#Raw read depth:4" in counts
        assert "#Final read-depth:3" in counts  # skipped pairs are not subtracted from the depth
        assert any("SKIPPED_READ_PAIRS count=2" in m for m in reader._mut_log.messages)

    def test_does_not_change_the_counts_header(self, tmp_path, monkeypatch):
        self.run(tmp_path, monkeypatch, [SKIP_1, SKIP_2, NO_MUT, UNMAPPED])
        counts = (tmp_path / "counts_sample_S1.csv").read_text().splitlines()
        assert [line.split(":")[0] for line in counts if line.startswith("#")] == HEADER_LABELS

    def test_survives_a_failure_to_write_the_skipped_file(self, tmp_path, monkeypatch):
        def fail(path, skipped):
            raise OSError("disk full")

        monkeypatch.setattr(count_mut, "write_skipped_reads", fail)
        reader = self.run(tmp_path, monkeypatch, [SKIP_1, NO_MUT])
        assert (tmp_path / "coverage_S1.csv").exists()
        assert not (tmp_path / "skippedreads_S1.csv").exists()
        assert any("Could not write" in m for m in reader._mut_log.messages)

    def test_writes_no_skipped_file_when_nothing_is_skipped(self, tmp_path, monkeypatch):
        reader = self.run(tmp_path, monkeypatch, [NO_MUT, UNMAPPED])
        assert not (tmp_path / "skippedreads_S1.csv").exists()
        assert not any("SKIPPED" in m for m in reader._mut_log.messages)
        assert (tmp_path / "coverage_S1.csv").exists()