#! /usr/bin/env python
from __future__ import print_function

import argparse
import logging

from mykrobe import K
from mykrobe.base import force_mixin
from mykrobe.base import genotyping_mixin
from mykrobe.base import panels_mixin
from mykrobe.base import sequence_or_binary_parser_mixin
from mykrobe.version import __version__

logging.basicConfig(
    level=logging.INFO,
    format="[mykrobe %(asctime)s %(levelname)s] %(message)s",
    datefmt="%Y-%m-%dT%H:%M:%S",
)
logger = logging.getLogger()


def run_subtool(parser, args):
    logger.debug(f"command: {args.command}")
    logger.debug(f"args: {args}")
    if args.command == "variants":
        if args.sub_command == "add":
            from mykrobe.cmds.variants.add import run
        elif args.sub_command == "dump-probes":
            from mykrobe.cmds.dump import run
        elif args.sub_command == "make-probes":
            from mykrobe.cmds.makeprobes import run
    elif args.command == "predict":
        from mykrobe.cmds.amr import run
    elif args.command == "panels":
        if args.sub_command == "describe":
            from mykrobe.cmds.panels import describe as run
        elif args.sub_command == "update_metadata":
            from mykrobe.cmds.panels import update_metadata as run
        elif args.sub_command == "update_species":
            from mykrobe.cmds.panels import update_species as run
    elif args.command == "genotype":
        raise NotImplementedError(
            "The 'genotype' option is no longer available - please use 'predict' with "
            "the species as 'custom' and the option --custom_probe_set_path "
            "probes.fasta (and also see the options --custom_variant_to_resistance_json, "
            "--custom_lineage_json) "
        )
    else:
        raise NotImplementedError(f"{args.command} is not a recognised mykrobe command")

    # run the chosen submodule.
    logger.debug(f"Start run {run.__module__}")
    run(parser, args)
    logger.debug(f"End run {run.__module__}")


class LogLevelWarn(argparse.Action):
    def __init__(self, nargs=0, **kw):
        if nargs != 0:
            raise ValueError("nargs for LogLevelWarn must be 0, it's a flag")
        super().__init__(nargs=nargs, **kw)

    def __call__(self, parser, namespace, values, option_string=None):
        logger.setLevel(logging.WARN)


class LogLevelDebug(argparse.Action):
    def __init__(self, nargs=0, **kw):
        if nargs != 0:
            raise ValueError("nargs for LogLevelDebug must be 0, it's a flag")
        super().__init__(nargs=nargs, **kw)

    def __call__(self, parser, namespace, values, option_string=None):
        logger.setLevel(logging.DEBUG)


class ArgumentParserWithDefaults(argparse.ArgumentParser):
    def __init__(self, *args, **kwargs):
        super(ArgumentParserWithDefaults, self).__init__(*args, **kwargs)
        self.add_argument(
            "-q",
            "--quiet",
            help="Only output warnings/errors to stderr",
            action=LogLevelWarn,
            nargs=0,
        )
        self.add_argument(
            "-d",
            "--debug",
            help="Output debugging information to stderr",
            action=LogLevelDebug,
            nargs=0,
        )


#########################################
# create the top-level parser
#########################################
parser = argparse.ArgumentParser(
    prog="mykrobe", formatter_class=argparse.ArgumentDefaultsHelpFormatter
)
parser.add_argument(
    "--version",
    help="mykrobe version",
    action="version",
    version="mykrobe " + str(__version__),
)
subparsers = parser.add_subparsers(
    title="[sub-commands]", dest="command", parser_class=ArgumentParserWithDefaults
)

db_parser_mixin = argparse.ArgumentParser(add_help=False)
db_parser_mixin.add_argument(
    "--db_name", metavar="db_name", type=str, help="db_name", default="default"
)
db_parser_mixin.add_argument(
    "--db_uri", metavar="db_uri", type=str, help="e.g.: mongodb://admin:qwerty@localhost/db_name", default=None
)

# ##########
# # AMR predict
# ##########

parser_amr = subparsers.add_parser(
    "predict",
    parents=[
        sequence_or_binary_parser_mixin,
        force_mixin,
        genotyping_mixin,
        panels_mixin,
    ],
    help="Predict the sample's drug susceptibility",
)
parser_amr.add_argument(
    "-S",
    "--species",
    metavar="species",
    type=str,
    required=True,
    help="Species name, or 'custom' to use custom data, in which case --custom_probe_set_path is required. Run `mykrobe panels describe` to see list of options [REQUIRED]",
)
parser_amr.add_argument(
    "--panel",
    metavar="panel",
    type=str,
    help="Name of panel to use. Ignored if species is 'custom'. Run `mykrobe panels describe` to see list of options",
)
parser_amr.add_argument(
    "-P",
    "--custom_probe_set_path",
    metavar="FILENAME",
    type=str,
    help="Required if species is 'custom'. Ignored otherwise. File path to fasta file from `mykrobe make-probes`.",
    default=None,
)
parser_amr.add_argument(
    "-R",
    "--custom_variant_to_resistance_json",
    metavar="FILENAME",
    type=str,
    help="For use with `--panel custom`. Ignored otherwise. File path to JSON with key,value pairs of variant names and induced drug resistance.",
    default=None,
)
parser_amr.add_argument(
    "-L",
    "--custom_lineage_json",
    metavar="FILENAME",
    type=str,
    help="For use with `--panel custom`. Ignored otherwise. File path to JSON made by --lineage option of make-probes",
    default=None,
)
parser_amr.add_argument(
    "--min_depth",
    metavar="min_depth",
    type=int,
    help="Minimum depth (default: %(default)d)",
    default=1,
)
parser_amr.add_argument(
    "--conf_percent_cutoff",
    metavar="conf_percent_cutoff",
    type=float,
    help=(
        "Number between 0 and 100. Determines --min_variant_conf, by simulating "
        "variants and choosing the cutoff that would keep x%% of the variants "
        "(default: %(default)d)"
    ),
    default=100,
)
parser_amr.add_argument(
    "-O",
    "--format",
    dest="output_format",
    type=str,
    help="Choose output format (default: %(default)s)",
    choices=["json", "csv", "json_and_csv"],
    default="csv",
)
parser_amr.set_defaults(func=run_subtool)


# ##################
# ###  Panels   ####
# ##################
parser_panels = subparsers.add_parser("panels", help="Add, update, or remove panels",)
panels_subparsers = parser_panels.add_subparsers(
    title="[sub-commands]", dest="sub_command", help="help"
)

# -------- panels describe -------------#
parser_describe = panels_subparsers.add_parser(
    "describe", help="Describe all known panels", parents=[panels_mixin],
)
parser_describe.set_defaults(func=run_subtool)

# -------- panels update_metadata ------#
parser_update_metadata = panels_subparsers.add_parser(
    "update_metadata",
    help="Update metadata about available species and their panels",
    parents=[panels_mixin],
)
parser_update_metadata.add_argument(
    "--filename", help=argparse.SUPPRESS,  # for testing. Keep hidden from the user.
)
parser_update_metadata.set_defaults(func=run_subtool)


# -------- panels update_species -------#
parser_update_species = panels_subparsers.add_parser(
    "update_species", help="Update species panel(s)", parents=[panels_mixin],
)
parser_update_species.add_argument(
    "species", help="Name of species to update, or 'all' to update all species",
)
parser_update_species.add_argument(
    "--remove", help="Remove species instead of updating", action="store_true",
)
parser_update_species.set_defaults(func=run_subtool)


# ##################
# ### Variants ##
# ##################

parser_variants = subparsers.add_parser(
    "variants", help="build variant probes", aliases=["vars"]
)

variant_subparsers = parser_variants.add_subparsers(
    title="[sub-commands]", dest="sub_command", help="help"
)

##########
# Add
##########
parser_add = variant_subparsers.add_parser(
    "add",
    help="adds a set of variants to the database",
    parents=[db_parser_mixin, force_mixin],
)
parser_add.add_argument("vcf", type=str, help="a vcf file")
parser_add.add_argument("reference_set", type=str, help="reference set")
parser_add.add_argument(
    "-m",
    "--method",
    type=str,
    help="variant caller method (e.g. CORTEX)",
    default="NotSpecified",
)
parser_add.set_defaults(func=run_subtool)

# # ##########
# # # Dump panel
# # ##########
parser_dump = variant_subparsers.add_parser(
    "dump-probes",
    help="dump a probe set of variant alleles from VCFs stored in database",
    parents=[db_parser_mixin, force_mixin],
)
parser_dump.add_argument(
    "reference_filepath",
    metavar="reference_filepath",
    type=str,
    help="reference_filepath",
)
parser_dump.add_argument(
    "--kmer", metavar="kmer", type=int, help="kmer length", default=K
)
parser_dump.add_argument("-v", "--verbose", default=False, action="store_true")
parser_dump.set_defaults(func=run_subtool)


# ##################
# ### Make Probes ##
# ##################

parser_make_probes = variant_subparsers.add_parser(
    "make-probes", help="make probes from a list of variants", parents=[db_parser_mixin]
)
parser_make_probes.add_argument(
    "reference_filepath",
    metavar="reference_filepath",
    type=str,
    help="reference_filepath",
)
parser_make_probes.add_argument(
    "-f", "--vcf", type=str, help="Use variants defined in a VCF file", default=[]
)
parser_make_probes.add_argument(
    "-v",
    "--variants",
    type=str,
    action="append",
    help="Variant in DNA positions e.g. A1234T",
    default=[],
)
parser_make_probes.add_argument(
    "-t", "--text_file", type=str, help="Text file containing variants as rows A1234T"
)
parser_make_probes.add_argument(
    "-g", "--genbank", type=str, help="Genbank file containing genes as features"
)
parser_make_probes.add_argument(
    "-k", "--kmer", type=int, help="kmer length", default=K
)
parser_make_probes.add_argument(
    "--no-backgrounds",
    help="Build probe set against reference only ignoring nearby variants",
    default=False,
    action="store_true",
)
parser_make_probes.add_argument(
    "--lineage", type=str, help="Write lineages to output file", metavar="FILENAME",
)
parser_make_probes.set_defaults(func=run_subtool)
