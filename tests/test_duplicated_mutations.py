"""
Tests for the DuplicatedMutations handling.

Background: an insertion immediately followed by a deletion (CIGAR ...nI mD...) gets the same
position key in MutParser._parse_cigar_mdz, so MutParser._convert_mut raises DuplicatedMutations.
count_mut.process_wrapper catches it, and the read pair is skipped and recorded instead of
aborting the whole count.

Run with:  pytest -v tests/test_duplicated_mutations.py
"""

import logging
import re

import pandas as pd
import pytest

from TileSeqMut import count_mut, locate_mut
from TileSeqMut.locate_mut import DuplicatedMutations, MutParser, SameRefandAlt

INS_THEN_DEL = {
  "cigar": [("10", "M"), ("3", "I"), ("4", "D"), ("10", "M")],
  "mdz": "10^ATTC10",
  "read": "A" * 10 + "GAG" + "A" * 10,
  "pos": 491
}

DEL_THEN_INS = {
  "cigar": [("10", "M"), ("4", "D"), ("3", "I"), ("10", "M")],
  "mdz": "10^ATTC10",
  "read": "A" * 10 + "GAG" + "A" * 10,
  "pos": 491
}

NINE_BASE = {
  "cigar": [("20", "M"), ("9", "I"), ("9", "D"), ("20", "M")],
  "mdz": "20^CCATTCAAG20",
  "read": "A" * 20 + "TTACCTGGT" + "A" * 20,
  "pos": 480
}

TWO_INS_THEN_DEL = {
  "cigar": [("13", "M"), ("1", "I"), ("8", "M"), ("3", "I"), ("4", "D"), ("10", "M")],
  "mdz": "21^ATTC10",
  "read": "A" * 13 + "G" + "A" * 8 + "GAG" + "A" * 10,
  "pos": 480
}

REPORTED_READ = {
  "cigar": [("1", "S"), ("26", "M"), ("4", "I"), ("5", "D"), ("120", "M")],
  "mdz": "6G3G8G6^ACAAG29G10A21T57",
  "read": "N" + "A" * 26 + "TTTT" + "A" * 120,
  "pos": 100
}

INS_BEFORE_MISMATCH = {
  "cigar": [("10", "M"), ("1", "I"), ("8", "M")],
  "mdz": "10A7",
  "read": "C" * 10 + "G" + "C" * 8,
  "pos": 491
}

READ_ROW = {
  "pos_start_r1" : "491", 
  "cigar_r1" : "10M3I4D10M",
  "pos_start_r2" : "520",
  "cigar_r2" : "30M"
}


def make_parser():
  # skips __init__ 
  # test only use parsing helpers, only needs logger
  parser = MutParser.__new__(MutParser)
  parser._logging = logging.getLogger("test_duplicated_mutations")
  return parser


def parse(cigar, mdz, read, pos):
  qual = "I" * len(read)
  # no ref, setting that to None, since we don't need it for this test
  snp, delins, _pos_map = make_parser()._parse_cigar_mdz(cigar, mdz, None, read, pos, qual)
  return snp, delins


def positions(records):
  # indel records
  out = []
  for record in records:
    fields = record.split("|")
    pos = int(fields[0])
    mut_type = fields[2]
    out.append((pos, mut_type))
  return out


### tests global test var conventions ###
def test_ins_then_del_share_a_position_key():
  assert positions(parse(**INS_THEN_DEL)[1]) == [(501, "ins"), (501, "del")]


def test_del_then_ins_do_not_collide():
    assert positions(parse(**DEL_THEN_INS)[1]) == [(505, "ins"), (501, "del")]


def test_nine_base_ins_then_del_share_a_position_key():
    assert positions(parse(**NINE_BASE)[1]) == [(500, "ins"), (500, "del")]


def test_two_insertions_then_deletion():
    assert positions(parse(**TWO_INS_THEN_DEL)[1]) == [(493, "ins"), (501, "ins"), (501, "del")]


def test_reported_read_shares_a_position_key():
    assert positions(parse(**REPORTED_READ)[1]) == [(126, "ins"), (126, "del")]


### tests functions ###
def test_isolated_insertion_converts_cleanly():
   _, delins = parse([("10", "M"), ("2", "I"), ("10", "M")], "20", "A" * 10 + "GG" + "A" * 10, 491)
   assert len(make_parser()._convert_mut(delins, "r1")) == 1


def test_isolated_deletion_is_flattened_per_base():
    _, delins = parse([("10", "M"), ("3", "D"), ("10", "M")], "10^ATT10", "A" * 20, 491)
    assert list(make_parser()._convert_mut(delins, "r1")["pos"]) == ["501", "502", "503"]


def test_snp_is_parsed():
   snp, delins = parse([("20", "M")], "10A9", "A" * 10 + "C" + "A" * 9, 491)
   assert snp == ["501|A|C|I"]
   assert delins == []


### tests yet to pass, will fail if fixed ###
@pytest.mark.xfail(strict=True, raises=DuplicatedMutations,
                   reason="known: an insertion directly followed by a deletion shares a position key")
def test_convert_mut_accepts_ins_then_del():
    _, delins = parse(**INS_THEN_DEL)
    assert len(make_parser()._convert_mut(delins, "r1")) == 5  # 1 insertion + 4 deleted bases


@pytest.mark.xfail(strict=True, raises=DuplicatedMutations,
                   reason="known: an insertion directly followed by a deletion shares a position key")
def test_convert_mut_accepts_two_ins_then_del():
    _, delins = parse(**TWO_INS_THEN_DEL)
    assert len(make_parser()._convert_mut(delins, "r1")) == 6  # 2 insertions + 4 deleted bases


### test count_mut.write_skipped_reads ###
SKIPPED_PAIRS = [
   ("491", "10M3I4D10M", "520", "30M"),
   ("466", "11S33M9I9D63M35S", "467", "104M47S")
]


def test_write_skipped_reads(tmp_path):
    path = tmp_path / "skippedreads_S1.csv"
    assert count_mut.write_skipped_reads(str(path), SKIPPED_PAIRS) == 2
    assert path.read_text().splitlines() == [
        "r1_start,r1_cigar,r2_start,r2_cigar",
        "491,10M3I4D10M,520,30M",
        "466,11S33M9I9D63M35S,467,104M47S",
    ]


def test_write_skipped_reads_with_no_pairs_writes_only_the_header(tmp_path):
    path = tmp_path / "skippedreads_S1.csv"
    assert count_mut.write_skipped_reads(str(path), []) == 0
    assert path.read_text().splitlines() == ["r1_start,r1_cigar,r2_start,r2_cigar"]


### test count_mut.process_wrapper ###
def fake_parser(behaviour):
    class FakeParser:
        def __init__(self, *args, **kwargs):
            pass
        
        def _main(self):
            return behaviour()
    return FakeParser


def run_wrapper():
  return count_mut.process_wrapper(
     READ_ROW, None, None, None, 1, 100, 0.9,
     logging.getLogger("test_duplicated_mutations"),
     0.0025, 33, False, False, []
  )


def test_process_wrapper_skips_pair_on_duplicated_mutations(monkeypatch):
    def raise_duplicated():
        raise DuplicatedMutations()

    monkeypatch.setattr(locate_mut, "MutParser", fake_parser(raise_duplicated))
    result = run_wrapper()
    assert len(result) == 7
    hgvs, outside_mut, all_df, r1_clusters, r2_clusters, track_df, skipped_pair = result
    assert (hgvs, outside_mut, r1_clusters, r2_clusters) == ([], [], [], [])
    assert all_df.empty and track_df.empty
    assert skipped_pair == ("491", "10M3I4D10M", "520", "30M")


def test_process_wrapper_passes_normal_results_through(monkeypatch):
    all_df, track_df = pd.DataFrame({"a": [1]}), pd.DataFrame({"b": [2]})
    monkeypatch.setattr(locate_mut, "MutParser",
                        fake_parser(lambda: ("c.1A>G", [5], all_df, ["c1"], ["c2"], track_df)))
    hgvs, outside_mut, all_out, r1_clusters, r2_clusters, track_out, skipped_pair = run_wrapper()
    assert (hgvs, outside_mut, r1_clusters, r2_clusters) == ("c.1A>G", [5], ["c1"], ["c2"])
    assert all_out is all_df and track_out is track_df
    assert skipped_pair is None


def test_process_wrapper_does_not_swallow_other_errors(monkeypatch):
    def raise_value_error():
        raise ValueError("unrelated failure")

    monkeypatch.setattr(locate_mut, "MutParser", fake_parser(raise_value_error))
    with pytest.raises(ValueError):
        run_wrapper()


def test_process_wrapper_does_not_swallow_same_ref_and_alt(monkeypatch):
    def raise_same_ref_alt():
        raise SameRefandAlt(pd.DataFrame({"ref": ["A"], "alt": ["A"]}))

    monkeypatch.setattr(locate_mut, "MutParser", fake_parser(raise_same_ref_alt))
    with pytest.raises(SameRefandAlt):
        run_wrapper()


### test count_mut.readSam.multi_core ###
# Only pairs that are skipped, unmapped or mutation-free are used here, so the
# per-position coverage merging at the end of the job loop is not exercised.
HEADER_LABELS = ["#Sample", "#Raw read depth", "#Number of read pairs without mutations",
                 "#Number of read pairs did not map to gene", "#Mutation positions outside of the tile",
                 "#Number of reads outside of the tile", "#Final read-depth",
                 "#Total read pairs with mutations", "#Comment", "#Comment",
                 "#Read len for R1", "#Read len for R2"]


class ListLogger:
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


def sam_line(name, cigar, md, pos=1, rname="ref"):
    length = sum(int(n) for n, op in re.findall(r"(\d+)([MIDS])", cigar) if op in "MIS")
    return f"{name}\t0\t{rname}\t{pos}\t60\t{cigar}\t*\t0\t0\t{'A' * length}\t{'I' * length}\tMD:Z:{md}\n"


def pair(name, r1_cigar, r1_md, r2_cigar="20M", r2_md="20"):
    return sam_line(name, r1_cigar, r1_md), sam_line(name, r2_cigar, r2_md)


SKIP_1 = pair("p1", "10M3I4D10M", "10^ATTC10")
SKIP_2 = pair("p2", "20M9I9D20M", "20^CCATTCAAG20")
NO_MUT = pair("p3", "20M", "20")
UNMAPPED = ("p4\t4\t*\t0\t0\t*\t*\t0\t0\tAAAA\tIIII\n", sam_line("p4", "20M", "20"))


def make_reader(tmp_path, pairs):
    r1, r2 = tmp_path / "S1_R1.sam", tmp_path / "S1_R2.sam"
    r1.write_text("".join(p[0] for p in pairs))
    r2.write_text("".join(p[1] for p in pairs))
    reader = count_mut.readSam.__new__(count_mut.readSam)
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


def run_multi_core(tmp_path, monkeypatch, pairs):
    monkeypatch.setattr(count_mut.mp, "Pool", FakePool)
    monkeypatch.setattr(locate_mut, "MutParser", FakeMutParser)
    reader = make_reader(tmp_path, pairs)
    reader.multi_core([])
    return reader


def test_multi_core_records_skipped_pairs_and_still_writes_outputs(tmp_path, monkeypatch):
    reader = run_multi_core(tmp_path, monkeypatch, [SKIP_1, SKIP_2, NO_MUT, UNMAPPED])
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


def test_multi_core_does_not_change_the_counts_header(tmp_path, monkeypatch):
    run_multi_core(tmp_path, monkeypatch, [SKIP_1, SKIP_2, NO_MUT, UNMAPPED])
    counts = (tmp_path / "counts_sample_S1.csv").read_text().splitlines()
    assert [line.split(":")[0] for line in counts if line.startswith("#")] == HEADER_LABELS


def test_multi_core_survives_a_failure_to_write_the_skipped_file(tmp_path, monkeypatch):
    def fail(path, skipped):
        raise OSError("disk full")

    monkeypatch.setattr(count_mut, "write_skipped_reads", fail)
    reader = run_multi_core(tmp_path, monkeypatch, [SKIP_1, NO_MUT])
    assert (tmp_path / "coverage_S1.csv").exists()
    assert not (tmp_path / "skippedreads_S1.csv").exists()
    assert any("Could not write" in m for m in reader._mut_log.messages)


def test_multi_core_writes_no_skipped_file_when_nothing_is_skipped(tmp_path, monkeypatch):
    reader = run_multi_core(tmp_path, monkeypatch, [NO_MUT, UNMAPPED])
    assert not (tmp_path / "skippedreads_S1.csv").exists()
    assert not any("SKIPPED" in m for m in reader._mut_log.messages)
    assert (tmp_path / "coverage_S1.csv").exists()