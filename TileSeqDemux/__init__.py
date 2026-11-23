"""tileseq_demux: a preprocessor that splits FASTQ pairs in which two tiles share one sequencing index.

It runs before tileseq_mut and shares no code with it; the two programs meet only at the parameter sheet and the
output folder. Nothing in this package imports TileSeqMut.
"""