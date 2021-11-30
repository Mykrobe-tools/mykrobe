from __future__ import print_function

import argparse
import os

from mykrobe import K, ONT_E_RATE, ONT_PLOIDY, ILLUMINA_E_RATE

sequence_or_graph_parser_mixin = argparse.ArgumentParser(add_help=False)
sequence_or_graph_parser_mixin.add_argument(
    "-s", "--sample", required=True, type=str, help="Sample identifier [REQUIRED]"
)
sequence_or_graph_parser_mixin.add_argument(
    "-k",
    "--kmer",
    metavar="kmer",
    type=int,
    help="K-mer length (default: %(default)d)",
    default=K,
)
sequence_or_graph_parser_mixin.add_argument(
    "--tmp", help="Directory to write temporary files to", default=None
)
sequence_or_graph_parser_mixin.add_argument(
    "--keep_tmp", help="Don't remove temporary files", action="store_true"
)
sequence_or_graph_parser_mixin.add_argument(
    "--skeleton_dir",
    help="Directory for skeleton binaries",
    default="mykrobe/data/skeletons/",
)
sequence_or_graph_parser_mixin.add_argument(
    "-t", "--threads", type=int, help="Number of threads to use", default=1
)
sequence_or_graph_parser_mixin.add_argument(
    "-m",
    "--memory",
    type=str,
    help="Memory to allocate for graph constuction (default: %(default)s)",
    default="1GB",
)
sequence_or_graph_parser_mixin.add_argument(
    "--expected_depth", type=int, help="Expected depth", default=None
)


SEQUENCE_FILES_HELP_STRING = "Sequence files (fasta,fastq,bam)"

sequence_or_binary_parser_mixin = argparse.ArgumentParser(
    parents=[sequence_or_graph_parser_mixin], add_help=False
)
sequence_or_binary_parser_mixin.add_argument(
    "-1",
    "-i",
    "--seq",
    metavar="seq",
    type=str,
    nargs="+",
    help=SEQUENCE_FILES_HELP_STRING,
)
sequence_or_binary_parser_mixin.add_argument(
    "-c", "--ctx", metavar="ctx", type=str, help="Cortex graph binary"
)

probe_set_mixin = argparse.ArgumentParser(add_help=False)
probe_set_mixin.add_argument(
    "probe_set", metavar="probe_set", type=str, help="probe_set"
)

force_mixin = argparse.ArgumentParser(add_help=False)
force_mixin.add_argument(
    "-f", "--force", action="store_true", help="Force override any skeleton files"
)

genotyping_mixin = argparse.ArgumentParser(add_help=False)
genotyping_mixin.add_argument(
    "--ont",
    action="store_true",
    help=f"Set defaults for ONT data. Sets `-e {ONT_E_RATE} --ploidy {ONT_PLOIDY}`",
)
genotyping_mixin.add_argument(
    "--guess_sequence_method",
    action="store_true",
    help="Guess if ONT or Illumia based on error rate. If error rate is > 10%%, ploidy is set to haploid and a confidence threshold is used ",
)
genotyping_mixin.add_argument(
    "--ignore_minor_calls",
    action="store_true",
    help="Ignore minor calls when running resistance prediction",
)
genotyping_mixin.add_argument(
    "--ignore_filtered", help="Don't include filtered genotypes", default=False
)
genotyping_mixin.add_argument(
    "--model",
    metavar="model",
    choices=["median_depth", "kmer_count"],
    type=str,
    help="Genotype model used. Options kmer_count, median_depth (default: %(default)s)",
    default="kmer_count",
)
genotyping_mixin.add_argument(
    "--ploidy",
    metavar="ploidy",
    choices=["diploid", "haploid"],
    type=str,
    help="Use a diploid (includes 0/1 calls) or haploid genotyping model (default: %(default)s)",
    default="diploid",
)
genotyping_mixin.add_argument(
    "--filters",
    help="Don't include specific filtered genotypes (default: %(default)s)",
    nargs="+",
    default=["MISSING_WT", "LOW_PERCENT_COVERAGE", "LOW_GT_CONF", "LOW_TOTAL_DEPTH"],
    required=False,
)
genotyping_mixin.add_argument(
    "-A",
    "--report_all_calls",
    help="Report all calls",
    action="store_true",
    default=False,
)
genotyping_mixin.add_argument(
    "-e",
    "--expected_error_rate",
    help="Expected sequencing error rate (default: %(default).3f)",
    default=ILLUMINA_E_RATE,
    type=float,
)
genotyping_mixin.add_argument(
    "--min_variant_conf",
    help="Minimum genotype confidence for variant genotyping (default: %(default)d)",
    default=150,
    type=int,
)
genotyping_mixin.add_argument(
    "--min_gene_conf",
    help="Minimum genotype confidence for gene genotyping (default: %(default)d)",
    default=1,
    type=int,
)
genotyping_mixin.add_argument(
    "-D",
    "--min_proportion_expected_depth",
    help="Minimum depth required on the sum of both alleles (default: %(default).2f)",
    default=0.3,
    type=float,
)
genotyping_mixin.add_argument(
    "--min_gene_percent_covg_threshold",
    help="All genes alleles found above this percent coverage will be reported (default: %(default)d (only best alleles reported))",
    default=100,
    type=int,
)
genotyping_mixin.add_argument(
    "-o",
    "--output",
    type=str,
    help="File path to save output file as. Default is to stdout",
    default="",
)


panels_mixin = argparse.ArgumentParser(add_help=False)
panels_mixin.add_argument(
    "--panels_dir",
    metavar="DIRNAME",
    help="Name of directory that contains panel data (default: %(default)s)",
    default=os.path.abspath(os.path.join(os.path.dirname(__file__), "data")),
)
