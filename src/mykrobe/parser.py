#! /usr/bin/env python
from __future__ import print_function
import os
import sys
import logging
from mykrobe.version import __version__
from mykrobe.base import DEFAULT_DB_NAME
from mykrobe.base import sequence_parser_mixin
from mykrobe.base import sequence_or_binary_parser_mixin
from mykrobe.base import probe_set_mixin
from mykrobe.base import force_mixin
from mykrobe.base import genotyping_mixin

import argparse
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def run_subtool(parser, args):
    if args.command == "variants":
        if args.sub_command == 'add':
            from mykrobe.cmds.variants.add import run
    # elif args.command == "add-gt":
    #     from mykrobe.cmds.atlasadd import run
        elif args.sub_command == "dump-probes":
            from mykrobe.cmds.dump import run
        elif args.sub_command == "make-probes":
            from mykrobe.cmds.makeprobes import run
    elif args.command == "predict":
        from mykrobe.cmds.amr import run
    elif args.command == "panels":
        from mykrobe.cmds.amr import describe_panels as run        
    elif args.command == "genotype":
        from mykrobe.cmds.genotype import run

    elif args.command == "place":
        from mykrobe.cmds.place import run
    elif args.command == "diff":
        from mykrobe.cmds.diff import run
    # run the chosen submodule.
    run(parser, args)

class ArgumentParserWithDefaults(argparse.ArgumentParser):

    def __init__(self, *args, **kwargs):
        super(ArgumentParserWithDefaults, self).__init__(*args, **kwargs)
        self.add_argument(
            "-q",
            "--quiet",
            help="do not output warnings to stderr",
            action="store_true",
            dest="quiet")
#########################################
# create the top-level parser
#########################################
parser = argparse.ArgumentParser(
    prog='mykrobe-atlas',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--version", help="mykrobe-atlas version",
                    action="version",
                    version="%%(prog)s " + str(__version__))
subparsers = parser.add_subparsers(
    title='[sub-commands]',
    dest='command',
    parser_class=ArgumentParserWithDefaults
)

db_parser_mixin = argparse.ArgumentParser(add_help=False)
db_parser_mixin.add_argument(
    '--db_name',
    metavar='db_name',
    type=str,
    help='db_name',
    default="default")

# ##########
# # AMR predict
# ##########

parser_amr = subparsers.add_parser(
    'predict',
    parents=[sequence_or_binary_parser_mixin,
             force_mixin, genotyping_mixin],
    help="predict the sample's drug susceptibility")
parser_amr.add_argument(
    'species',
    metavar='species',
    choices=['staph', 'tb'],
    type=str,
    help='species')
parser_amr.add_argument(
    '--panel',
    metavar='panel',
    type=str,
    help='variant panel (default:201901). custom requires custom_probe_set_path and custom_variant_to_resistance_json to be set',
    choices=['bradley-2015', 'walker-2015', '201901', 'atlas', 'custom'],
    default='201901')
parser_amr.add_argument(
    '--custom_probe_set_path',
    metavar='custom_probe_set_path',
    type=str,
    help='For use with `--panel custom`. File path to fasta file from `mykrobe make-probes`.',
    default=None)
parser_amr.add_argument(
    '--custom_variant_to_resistance_json',
    metavar='custom_variant_to_resistance_json',
    type=str,
    help='For use with `--panel custom`. File path to JSON with key,value pairs of variant names and induced drug resistance.',
    default=None)
parser_amr.add_argument(
    '--min_depth',
    metavar='min_depth',
    type=int,
    help='min_depth',
    default=1)
parser_amr.add_argument(
    '--conf_percent_cutoff',
    metavar='conf_percent_cutoff',
    type=float,
    help='Number between 0 and 100. Determines --min_variant_conf, by simulating variants and choosing the cutoff that would keep x%% of the variants. Default is 90 if --ont, otherwise --min_variant_conf is used as the cutoff',
    default=-1)
parser_amr.add_argument(
    '--format',
    dest='output_format',
    type=str,
    help='Choose output format. Default: csv.',
    choices=['json', 'csv'],
    default='csv')
parser_amr.set_defaults(func=run_subtool)

parser_panels = subparsers.add_parser(
    'panels',
    help="A description of the AMR panels available within Mykrobe predict")
parser_panels.set_defaults(func=run_subtool)

# ##################
# ### Variants ##
# ##################

parser_variants = subparsers.add_parser(
    'variants',
    help="build variant probes", aliases=['vars'])

variant_subparsers = parser_variants.add_subparsers(
    title='[sub-commands]',
    dest='sub_command',
    help="help"
)

##########
# Add
##########
parser_add = variant_subparsers.add_parser(
    'add',
    help='adds a set of variants to the database',
    parents=[db_parser_mixin, force_mixin])
parser_add.add_argument('vcf', type=str, help='a vcf file')
parser_add.add_argument('reference_set', type=str, help='reference set')
parser_add.add_argument(
    '-m',
    '--method',
    type=str,
    help='variant caller method (e.g. CORTEX)',
    default="NotSpecified")
parser_add.set_defaults(func=run_subtool)

# # ##########
# # # Dump panel
# # ##########
parser_dump = variant_subparsers.add_parser(
    'dump-probes',
    help='dump a probe set of variant alleles from VCFs stored in database',
    parents=[db_parser_mixin, force_mixin])
parser_dump.add_argument(
    'reference_filepath',
    metavar='reference_filepath',
    type=str,
    help='reference_filepath')
parser_dump.add_argument(
    '--kmer',
    metavar='kmer',
    type=int,
    help='kmer length',
    default=31)
parser_dump.add_argument(
    '-v',
    '--verbose',
    default=False,
    action="store_true")
parser_dump.set_defaults(func=run_subtool)

# parser_add_gt = subparsers.add_parser(
#     'add-gt',
#     help='adds a set of atlas genotype calls to the atlas',
#     parents=[db_parser_mixin, force_mixin])
# parser_add_gt.add_argument('jsons', type=str, nargs='+',
#                            help='json output from `atlas genotype`')
# parser_add_gt.add_argument(
#     '-m',
#     '--method',
#     type=str,
#     help='variant caller method (e.g. CORTEX)',
#     default="atlas")
# parser_add_gt.set_defaults(func=run_subtool)

# ##################
# ### Make Probes ##
# ##################

parser_make_probes = variant_subparsers.add_parser(
    'make-probes', help='make probes from a list of variants',
    parents=[db_parser_mixin])
parser_make_probes.add_argument(
    'reference_filepath',
    metavar='reference_filepath',
    type=str,
    help='reference_filepath')
parser_make_probes.add_argument(
    '-f',
    '--vcf',
    type=str,
    help='Use variants defined in a VCF file',
    default=[])
parser_make_probes.add_argument(
    '-v',
    '--variant',
    type=str,
    action='append',
    help='Variant in DNA positions e.g. A1234T',
    default=[])
parser_make_probes.add_argument(
    '-t',
    '--text_file',
    type=str,
    help='Text file containing variants as rows A1234T')
parser_make_probes.add_argument(
    '-g',
    '--genbank',
    type=str,
    help='Genbank file containing genes as features')
parser_make_probes.add_argument(
    '-k',
    '--kmer',
    type=int,
    help='kmer length',
    default=31)
parser_make_probes.add_argument(
    '--no-backgrounds',
    help='Build probe set against reference only ignoring nearby variants',
    default=False,
    action="store_true")
parser_make_probes.set_defaults(func=run_subtool)

# ##########
# # Genotype
# ##########
parser_geno = subparsers.add_parser(
    'genotype',
    parents=[
        sequence_or_binary_parser_mixin,
        probe_set_mixin,
        force_mixin,
        genotyping_mixin],
    help='genotype a sample using a probe set')
parser_geno.set_defaults(func=run_subtool)
