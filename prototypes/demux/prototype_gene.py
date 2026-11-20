#!/usr/bin/env python3
"""Prototype for ticket 1.2, whole-gene version. Not production code.

Each mate is aligned to the whole template with tileseq_mut's bowtie2 flags (alignment.py). Each read pair is then
assigned to a tile with tileseq_mut's tile test (count_mut.py) plus an overlap check, in one streaming pass over the
FASTQs: the two bowtie2 outputs (kept in input order with --reorder) are joined to the FASTQ records as they stream
past, so no read IDs are held in memory."""
import argparse
import csv
import glob
import gzip
import math
import os
import re
import resource
import subprocess
import sys
import time
from collections import Counter
from contextlib import ExitStack
from itertools import takewhile, zip_longest

BOWTIE = ["bowtie2", "--no-head", "--no-sq", "--no-unal", "--reorder", "--local", "--rdg", "12,1", "--rfg", "12,1"]
STRAND = {1: ["--norc"], 2: ["--nofw"]}      # alignment.py: R1 forward strand only, R2 reverse only (--rc removes this)
REF_OPS = b"MDN=X"                           # CIGAR operations that consume the reference
CIGAR = re.compile(rb"(\d+)([MIDNSHP=X])")
UNMAPPED = {(False, False): None, (True, False): "r1_unmapped", (False, True): "r2_unmapped", (True, True): "both_unmapped"}
AMBIGUOUS = {"drop": lambda hits: [], "both": lambda hits: hits}


# ---------------------------------------------------------------- sheet and reference
def read_sheet(path):
    """Template sequence, CDS start, {tile: (start_aa, end_aa)} and the minimum coverage from the parameter sheet CSV."""
    with open(path, newline="") as fh:
        rows = [[c.strip().strip('"') for c in r] for r in csv.reader(fh)]

    def find(label):
        return [r for r in rows if r and r[0] == label]

    def one(label):
        hits = find(label)
        if len(hits) != 1:
            sys.exit("expected exactly one %r row in %s, found %d" % (label, path, len(hits)))
        return hits[0]

    table = takewhile(lambda r: r and r[0], rows[rows.index(one("Tile Number")) + 1:])
    tiles = {r[0]: (int(r[1]), int(r[2])) for r in table}
    cover = find("Minimum coverage:")
    min_cover = float(cover[0][1]) if cover and len(cover[0]) > 1 and cover[0][1] else 0.6   # csv2json's default
    return one("Sequence:")[1].upper(), int(one("CDS start:")[1]), tiles, min_cover


def tile_geometry(start_aa, end_aa, min_cover):
    """tileseq_mut's arithmetic (count_mut.py): CDS positions of the tile and the minimum mapped length."""
    begins, ends = start_aa * 3 - 2, end_aa * 3
    return {"begins": begins, "ends": ends, "min_map": math.ceil((ends - begins) * min_cover)}


def build_index(name, seq, ref_dir):
    fasta = os.path.join(ref_dir, name + ".fasta")
    with open(fasta, "w") as fh:
        fh.write(">%s\n%s\n" % (name, seq))
    prefix = os.path.join(ref_dir, name)
    subprocess.run(["bowtie2-build", "--quiet", "-f", fasta, prefix], check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return prefix


# ---------------------------------------------------------------- alignment streams
def start_bowtie(index, fastq, mate, threads, rc, log_path):
    cmd = BOWTIE + ([] if rc else STRAND[mate]) + ["-p", str(threads), "-x", index, "-U", os.path.abspath(fastq)]
    log = open(log_path, "w")                 # stderr goes to a file, never into the SAM stream
    return subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=log, bufsize=1 << 20), log


def parse_sam(line):
    """(POS, length of SEQ, reference bases covered by the alignment) from one SAM line."""
    f = line.split(b"\t", 10)
    span = sum(int(n) for n, op in CIGAR.findall(f[5]) if op in REF_OPS)
    return int(f[3]), len(f[9]), span


class SamStream:
    """One bowtie2 output, read with one line of lookahead. Records come in input order (--reorder) and unaligned
    reads are absent (--no-unal), so a read is aligned exactly when the next record carries its name."""

    def __init__(self, proc):
        self.lines = iter(proc.stdout)
        self.count = 0
        self.line = self.name = None
        self._advance()

    def _advance(self):
        self.line = next(self.lines, None)
        self.name = self.line[:self.line.index(b"\t")] if self.line is not None else None
        if self.line is not None:
            self.count += 1

    def take(self, read_id):
        if self.name != read_id:
            return None
        aln = parse_sam(self.line)
        while self.name == read_id:          # skip any further lines for the same read
            self._advance()
        return aln


# ---------------------------------------------------------------- assignment rules
def mate_ok(aln, cds_start, tile):
    """One mate against one tile: tileseq_mut's end test plus the overlap check (the addition for pooled tiles)."""
    pos, seq_len, span = aln
    reaches = pos + seq_len - cds_start >= tile["begins"] + tile["min_map"]           # count_mut.py
    first = pos - cds_start + 1                                                       # CDS position of the first aligned base
    overlap = min(first + span - 1, tile["ends"]) - max(first, tile["begins"]) + 1
    return reaches and overlap >= tile["min_map"]


def classify(a1, a2, geom, cds_start, ambiguous):
    """(category, tiles to write the pair to)."""
    kind = UNMAPPED[(a1 is None, a2 is None)]
    if kind:
        return kind, []
    hits = [t for t, g in geom.items() if mate_ok(a1, cds_start, g) and mate_ok(a2, cds_start, g)]
    if len(hits) == 1:
        return "assigned:" + hits[0], hits
    if not hits:
        return "off_tile", []
    return "ambiguous", AMBIGUOUS[ambiguous](hits)


# ---------------------------------------------------------------- FASTQ
def records(path):
    """Yield raw FASTQ records as [header, sequence, plus, quality] byte lines."""
    with gzip.open(path, "rb") as fh:
        while True:
            rec = [fh.readline() for _ in range(4)]
            if not rec[0]:
                return
            if not rec[3]:
                raise ValueError("truncated FASTQ record in %s" % path)
            yield rec


def read_id(rec):
    return rec[0][1:].split(None, 1)[0]


def find_fastq(folder, physical, mate):
    hits = glob.glob(os.path.join(folder, "%s_*_R%d_*.fastq.gz" % (physical, mate)))
    if len(hits) != 1:
        sys.exit("expected one R%d file for %s in %s, found %s" % (mate, physical, folder, hits))
    return hits[0]


def open_outputs(stack, out_dir, plan):
    def path(logical, mate):
        return os.path.join(out_dir, "%s_S%s_R%d_001.fastq.gz" % (logical, logical, mate))
    return {tile: [stack.enter_context(gzip.open(path(logical, m), "wb", compresslevel=6)) for m in (1, 2)]
            for tile, logical in plan.items()}


def load_baseline(folder, plan):
    """{read ID: tile} from the baseline ID lists (<logical>.ids, lines start with '@')."""
    ids = {}
    for tile, logical in plan.items():
        with open(os.path.join(folder, logical + ".ids"), "rb") as fh:
            for line in fh:
                ids[line.rstrip(b"\r\n")[1:]] = tile
    return ids


# ---------------------------------------------------------------- the pass
def run_pass(fastq, streams, geom, cds_start, ambiguous, outs, baseline):
    s1, s2 = streams
    cats, written, xtab = Counter(), Counter(), Counter()
    n = 0
    for a, b in zip_longest(records(fastq[1]), records(fastq[2])):
        if a is None or b is None:
            raise ValueError("R1 and R2 contain different numbers of reads")
        rid = read_id(a)
        if rid != read_id(b):
            raise ValueError("R1/R2 read IDs differ at read %d: %r vs %r" % (n + 1, rid, read_id(b)))
        n += 1
        cat, tiles = classify(s1.take(rid), s2.take(rid), geom, cds_start, ambiguous)
        cats[cat] += 1
        for t in tiles:
            outs[t][0].write(b"".join(a))
            outs[t][1].write(b"".join(b))
            written[t] += 1
        if baseline is not None:
            xtab[(baseline.get(rid, "-"), cat)] += 1
    return n, cats, written, xtab


def check_stream(proc, log, stream, log_path, mate, n):
    """The join is only valid if bowtie2's output was fully consumed and agrees with its own summary."""
    if stream.line is not None:
        proc.kill()
        raise ValueError("R%d: bowtie2 output not in input order or read names do not match (left over: %r)"
                         % (mate, stream.line[:80]))
    code = proc.wait()
    log.close()
    if code:
        raise RuntimeError("bowtie2 for R%d exited %d, see %s" % (mate, code, log_path))
    text = open(log_path).read()
    total = int(re.search(r"(\d+) reads; of these:", text).group(1))
    aligned = sum(int(x) for x in re.findall(r"(\d+) \([\d.]+%\) aligned (?:exactly 1 time|>1 times)", text))
    if total != n or aligned != stream.count:
        raise ValueError("R%d: bowtie2 reports %d reads and %d aligned, but the pass saw %d reads and %d SAM lines"
                         % (mate, total, aligned, n, stream.count))
    return aligned


def print_crosstab(xtab, cats):
    cols = sorted(cats)
    rows = sorted({b for b, _ in xtab})
    print("baseline tile -> category of the same pair in this run (\"-\" = not in the baseline)")
    print("  %-8s" % "tile" + "".join(" %14s" % c[-14:] for c in cols))
    for r in rows:
        print("  %-8s" % r + "".join(" %14d" % xtab.get((r, c), 0) for c in cols))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sheet", required=True)
    ap.add_argument("--fastq-dir", required=True)
    ap.add_argument("--physical", required=True)
    ap.add_argument("--tile", action="append", required=True, metavar="TILE:LOGICAL")
    ap.add_argument("--baseline-dir")
    ap.add_argument("--out", required=True)
    ap.add_argument("--threads", type=int, default=1, help="threads per bowtie2 process (there are two)")
    ap.add_argument("--rc", action="store_true", help="align both strands, like tileseq_mut -rc")
    ap.add_argument("--ambiguous", choices=sorted(AMBIGUOUS), default="drop", help="pairs that pass two tiles")
    ap.add_argument("--min-cover", type=float, help="override the sheet's Minimum coverage")
    a = ap.parse_args()
    start = time.time()

    plan = dict(t.split(":") for t in a.tile)                       # {tile: logical sample ID}
    seq, cds_start, sheet_tiles, sheet_cover = read_sheet(a.sheet)
    unknown = sorted(set(plan) - set(sheet_tiles))
    if unknown:
        sys.exit("tiles not in sheet: %s (sheet has %s)" % (unknown, sorted(sheet_tiles, key=int)))
    min_cover = a.min_cover if a.min_cover is not None else sheet_cover
    geom = {t: tile_geometry(*sheet_tiles[t], min_cover) for t in plan}
    fastq = {m: find_fastq(a.fastq_dir, a.physical, m) for m in (1, 2)}

    ref_dir, log_dir, out_dir = (os.path.join(a.out, d) for d in ("ref", "logs", "tile_fastqs"))
    for d in (ref_dir, log_dir, out_dir):
        os.makedirs(d, exist_ok=True)
    index = build_index("gene", seq, ref_dir)
    baseline = load_baseline(a.baseline_dir, plan) if a.baseline_dir else None
    print("minimum coverage %.2f, CDS start %d; tiles: %s" % (min_cover, cds_start, {t: g for t, g in geom.items()}))

    logs = {m: os.path.join(log_dir, "gene_R%d.log" % m) for m in (1, 2)}
    procs = {m: start_bowtie(index, fastq[m], m, a.threads, a.rc, logs[m]) for m in (1, 2)}
    with ExitStack() as stack:
        outs = open_outputs(stack, out_dir, plan)
        streams = (SamStream(procs[1][0]), SamStream(procs[2][0]))
        n, cats, written, xtab = run_pass(fastq, streams, geom, cds_start, a.ambiguous, outs, baseline)
    aligned = [check_stream(procs[m][0], procs[m][1], streams[m - 1], logs[m], m, n) for m in (1, 2)]
    if sum(cats.values()) != n:
        raise ValueError("categories add up to %d, not %d" % (sum(cats.values()), n))

    print("reads in %d; aligned R1 %d, R2 %d; wall %.0f s" % (n, aligned[0], aligned[1], time.time() - start))
    with open(os.path.join(a.out, "summary.tsv"), "w") as fh:
        for cat in sorted(cats):
            print("  %-16s %10d  (%.2f%%)" % (cat, cats[cat], 100.0 * cats[cat] / n))
            fh.write("%s\t%d\n" % (cat, cats[cat]))
    print("written per tile:", dict(written))
    if baseline is not None:
        print_crosstab(xtab, cats)
    own, kids = resource.getrusage(resource.RUSAGE_SELF), resource.getrusage(resource.RUSAGE_CHILDREN)
    print("peak RSS: python %.0f MB, largest child (bowtie2) %.0f MB" % (own.ru_maxrss / 1024.0, kids.ru_maxrss / 1024.0))
    print("CPU seconds: python %.0f, bowtie2 %.0f; wall %.0f"
          % (own.ru_utime + own.ru_stime, kids.ru_utime + kids.ru_stime, time.time() - start))


if __name__ == "__main__":
    main()
