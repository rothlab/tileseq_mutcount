"""Labels and defaults of the parameter sheet, kept in one place.

The labels of the Index Samples section are provisional until confirmed.
Everything else matches what `tsm csv2json` reads, so one sheet serves both tools.
"""

LABEL_SEQUENCE = "Sequence:"
LABEL_CDS_START = "CDS start:"
LABEL_CDS_END = "CDS end:"
LABEL_MIN_COVER = "Minimum coverage:"
LABEL_TILES = "Tile Number"          # first cell of the tile table's header row
LABEL_SAMPLES = "Sample ID"          # first cell of the sample table's header row
LABEL_INDEX = "Index Samples"        # title row of the new section (provisional)

HEADER_INDEX = ("Shared Sample ID", "Unique Sample ID")     # header row under the title (provisional)
SECOND_COLUMN_SAMPLES = "Tile ID"

DEFAULT_MIN_COVER = 0.6              # what tsm csv2json uses when the sheet has no Minimum coverage row
ALLOWED_BASES = frozenset("ACGTN")