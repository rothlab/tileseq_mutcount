"""The FASTQ file naming contract that tileseq_mut expects, stated once.

tileseq_mut finds its input like this:
    - only files that end in .fastq.gz and sit directly inside the -f folder are seen
    - a sample's ID is the first "_" token of the file name
    - R1 is recognised by "_R1_" in the name, and its mate by "_R2_"

Everything in tileseq_demux that decides which file belongs to which sample, or names a file it writes, goes through
this module, so the rules live in one place. Standard library only.
"""
SUFFIX = ".fastq.gz"


def sample_id(file_name):
    """The sample ID tileseq_mut reads from a file name: its first "_" token. None when the name has no "_"."""
    return file_name.split("_")[0] if "_" in file_name else None


def is_mate_file(file_name, sample, mate):
    """True when `file_name` is the R1 (mate=1) or R2 (mate=2) file of `sample`."""
    after_id = file_name[len(sample) + 1:]
    return file_name.startswith(sample + "_") and file_name.endswith(SUFFIX) and f"_R{mate}_" in after_id


def tile_file_name(sample, mate):
    """The name tileseq_demux gives the R1 or R2 file it writes for `sample`: {ID}_S{ID}_R{mate}_001.fastq.gz."""
    return f"{sample}_S{sample}_R{mate}_001{SUFFIX}"


def group_by_sample_id(file_names):
    """File names grouped by sample ID, in the order given. Names with no "_" have no sample ID and are left out."""
    groups = {}
    for name in file_names:
        sample = sample_id(name)
        if sample is not None:
            groups.setdefault(sample, []).append(name)
    return groups