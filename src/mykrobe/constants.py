import os

DB_PREFIX = "atlas"
STOP = "*"
K = os.environ.get("KMER_SIZE", 21)
ONT_E_RATE = 0.08
ONT_PLOIDY = "haploid"
ILLUMINA_E_RATE = 0.05
NON_METADATA_KEYS = {"INFO", "FILTER", "ALT", "FORMAT", "contig"}
# Spec is a bit weak on which metadata lines are singular, like fileformat
# and which can have repeats, like contig
SINGULAR_METADATA = {'fileformat', 'fileDate', 'reference'}
MCCORTEX_BINARY_ENV_VAR = "MYKROBE_MCCORTEX"
