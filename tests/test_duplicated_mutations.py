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