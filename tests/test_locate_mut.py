# These tests were generated with Claude (Anthropic) and then reviewed, edited and run by the maintainers.
"""
Tests for how MutParser (locate_mut.py) keys insertions and deletions.

Background
    An insertion immediately followed by a deletion (CIGAR ...nI mD...) gets the same position
    key in MutParser._parse_cigar_mdz, so MutParser._convert_mut raises DuplicatedMutations.
    count_mut.process_wrapper handles that error; its tests are in test_count_mut.py.

Run with:  python -m pytest tests/test_locate_mut.py -v
"""

import logging
import pytest
from TileSeqMut.locate_mut import DuplicatedMutations, MutParser

### test data ###

# 3-base insertion directly followed by a 4-base deletion (the failing pattern)
INS_THEN_DEL = {
    "cigar": [("10", "M"), ("3", "I"), ("4", "D"), ("10", "M")],
    "mdz": "10^ATTC10",
    "read": "A" * 10 + "GAG" + "A" * 10,
    "pos": 491,
}

# the same events in the opposite order (does not collide)
DEL_THEN_INS = {
    "cigar": [("10", "M"), ("4", "D"), ("3", "I"), ("10", "M")],
    "mdz": "10^ATTC10",
    "read": "A" * 10 + "GAG" + "A" * 10,
    "pos": 491,
}

# 9-base insertion directly followed by a 9-base deletion
NINE_BASE = {
    "cigar": [("20", "M"), ("9", "I"), ("9", "D"), ("20", "M")],
    "mdz": "20^CCATTCAAG20",
    "read": "A" * 20 + "TTACCTGGT" + "A" * 20,
    "pos": 480,
}

# two insertions and a deletion in one read; only the second insertion collides
TWO_INS_THEN_DEL = {
    "cigar": [("13", "M"), ("1", "I"), ("8", "M"), ("3", "I"), ("4", "D"), ("10", "M")],
    "mdz": "21^ATTC10",
    "read": "A" * 13 + "G" + "A" * 8 + "GAG" + "A" * 10,
    "pos": 480,
}

# a read from the original bug report (starts with a soft clip)
REPORTED_READ = {
    "cigar": [("1", "S"), ("26", "M"), ("4", "I"), ("5", "D"), ("120", "M")],
    "mdz": "6G3G8G6^ACAAG29G10A21T57",
    "read": "N" + "A" * 26 + "TTTT" + "A" * 120,
    "pos": 100,
}

# reads with a single kind of event, which must keep working
ISOLATED_INSERTION = {
    "cigar": [("10", "M"), ("2", "I"), ("10", "M")],
    "mdz": "20",
    "read": "A" * 10 + "GG" + "A" * 10,
    "pos": 491,
}

ISOLATED_DELETION = {
    "cigar": [("10", "M"), ("3", "D"), ("10", "M")],
    "mdz": "10^ATT10",
    "read": "A" * 20,
    "pos": 491,
}

SNP_ONLY = {
    "cigar": [("20", "M")],
    "mdz": "10A9",
    "read": "A" * 10 + "C" + "A" * 9,
    "pos": 491,
}


### helper functions ###

def make_parser():
    # Skip __init__ (it needs a full reference/tile setup): these tests only use the
    # parsing helpers, which need nothing but a logger.
    parser = MutParser.__new__(MutParser)
    parser._logging = logging.getLogger("test_locate_mut")
    return parser


def parse(read):
    """(snp_list, delins_list) that MutParser._parse_cigar_mdz returns for an example read."""
    qual = "I" * len(read["read"])
    # ref is not used by the parser for these tests, so None
    snp, delins, _pos_map = make_parser()._parse_cigar_mdz(
        read["cigar"], read["mdz"], None, read["read"], read["pos"], qual)
    return snp, delins


def convert(delins):
    return make_parser()._convert_mut(delins, "r1")


def positions(records):
    """[(position, type), ...] for insertion/deletion records such as '501|GAG|ins|I,I,I'."""
    out = []
    for record in records:
        fields = record.split("|")
        out.append((int(fields[0]), fields[2]))
    return out



### tests ###

class TestPositionKeys:
    """Documents the key (position) the parser gives each insertion and deletion.

    If the convention is changed on purpose, update the expected values here.
    """

    @pytest.mark.parametrize("read, expected", [
        pytest.param(INS_THEN_DEL, [(501, "ins"), (501, "del")], id="ins_then_del_collide"),
        pytest.param(DEL_THEN_INS, [(505, "ins"), (501, "del")], id="del_then_ins_do_not_collide"),
        pytest.param(NINE_BASE, [(500, "ins"), (500, "del")], id="nine_base_ins_then_del_collide"),
        pytest.param(TWO_INS_THEN_DEL, [(493, "ins"), (501, "ins"), (501, "del")],
                     id="only_the_second_insertion_collides"),
        pytest.param(REPORTED_READ, [(126, "ins"), (126, "del")], id="reported_read_collides"),
    ])
    def test_position_keys(self, read, expected):
        assert positions(parse(read)[1]) == expected


class TestParsingStillWorks:
    """Reads with a single kind of event must keep working."""

    def test_isolated_insertion_converts_cleanly(self):
        _, delins = parse(ISOLATED_INSERTION)
        assert len(convert(delins)) == 1

    def test_isolated_deletion_is_flattened_per_base(self):
        _, delins = parse(ISOLATED_DELETION)
        assert list(convert(delins)["pos"]) == ["501", "502", "503"]

    def test_snp_is_parsed(self):
        snp, delins = parse(SNP_ONLY)
        assert snp == ["501|A|C|I"]
        assert delins == []


class TestKnownLimitation:
    """Behavior we want but don't have yet.

    These are strict xfails: once the key problem is fixed they will start passing,
    pytest reports that as a failure, and the xfail marker should then be removed.
    """

    @pytest.mark.xfail(strict=True, raises=DuplicatedMutations,
                       reason="known: an insertion directly followed by a deletion shares a position key")
    @pytest.mark.parametrize("read, expected_rows", [
        pytest.param(INS_THEN_DEL, 5, id="1_insertion_and_4_deleted_bases"),
        pytest.param(TWO_INS_THEN_DEL, 6, id="2_insertions_and_4_deleted_bases"),
    ])
    def test_convert_mut_accepts_insertion_then_deletion(self, read, expected_rows):
        _, delins = parse(read)
        assert len(convert(delins)) == expected_rows