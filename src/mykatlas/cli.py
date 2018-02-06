#! /usr/bin/env python
from __future__ import print_function
import os
import sys
import logging
import argparse


sys.path.append(
    os.path.realpath(
        os.path.join(
            os.path.dirname(__file__),
            "..")))


from mykatlas.version import __version__
from mykatlas.base import ArgumentParserWithDefaults
from mykatlas.base import DEFAULT_DB_NAME
from mykatlas.base import sequence_parser_mixin
from mykatlas.base import sequence_or_binary_parser_mixin
from mykatlas.base import probe_set_mixin
from mykatlas.base import force_mixin
from mykatlas.base import genotyping_mixin


logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
logger = logging.getLogger(__name__)

DEFAULT_KMER_SIZE = os.environ.get("KMER_SIZE", 31)


def run_subtool(parser, args):
    # if args.command == 'add':
    #     from mykatlas.cmds.add import run
    # elif args.command == "add-gt":
    #     from mykatlas.cmds.atlasadd import run
    # elif args.command == "dump-probes":
    #     from mykatlas.cmds.dump import run

    if args.command == "genotype":
        from mykatlas.cmds.genotype import run
    elif args.command == "make-probes":
        from mykatlas.cmds.makeprobes import run
    elif args.command == "walk":
        from mykatlas.cmds.walk import run
    elif args.command == "place":
        from mykatlas.cmds.place import run
    elif args.command == "diff":
        from mykatlas.cmds.diff import run
    # run the chosen submodule.
    run(parser, args)


def main():
    #########################################
    # create the top-level parser
    #########################################
    parser = argparse.ArgumentParser(
        prog='atlas',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--version", help="atlas version",
                        action="version",
                        version="%(prog)s " + str(__version__))
    subparsers = parser.add_subparsers(
        title='[sub-commands]',
        dest='command',
        parser_class=ArgumentParserWithDefaults)

    db_parser_mixin = argparse.ArgumentParser(add_help=False)
    db_parser_mixin.add_argument(
        '--db_name',
        metavar='db_name',
        type=str,
        help='db_name',
        default="default")

    # ##########
    # # Add
    # ##########
    # parser_add = subparsers.add_parser(
    #     'add',
    #     help='adds a set of variants to the atlas',
    #     parents=[db_parser_mixin, force_mixin])
    # parser_add.add_argument('vcf', type=str, help='a vcf file')
    # parser_add.add_argument('reference_set', type=str, help='reference set')
    # parser_add.add_argument(
    #     '-m',
    #     '--method',
    #     type=str,
    #     help='variant caller method (e.g. CORTEX)',
    #     default="NotSpecified")
    # parser_add.set_defaults(func=run_subtool)

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

    # # ##########
    # # # Dump panel
    # # ##########
    # parser_dump = subparsers.add_parser(
    #     'dump-probes',
    #     help='dump a probe set of variant alleles',
    #     parents=[db_parser_mixin, force_mixin])
    # parser_dump.add_argument(
    #     'reference_filepath',
    #     metavar='reference_filepath',
    #     type=str,
    #     help='reference_filepath')
    # parser_dump.add_argument(
    #     '--kmer',
    #     metavar='kmer',
    #     type=int,
    #     help='kmer length',
    #     default=31)
    # parser_dump.add_argument(
    #     '-v',
    #     '--verbose',
    #     default=False,
    #     action="store_true")
    # parser_dump.set_defaults(func=run_subtool)

    # ##################
    # ### Make Probes ##
    # ##################

    parser_make_probes = subparsers.add_parser(
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

    ##############
    ## Walk ##
    #############
    parser_walk = subparsers.add_parser(
        'walk',
        parents=[
            sequence_or_binary_parser_mixin,
            probe_set_mixin,
            force_mixin,
            genotyping_mixin],
        help='Walk through a graph using an existing sequence probe set as seeds. default walking algorithm is a depth first search')
    parser_walk.add_argument(
        '--also-genotype',
        default=False,
        action="store_true")
    parser_walk.add_argument('--show-all-paths', action="store_true")
    parser_walk.set_defaults(func=run_subtool)

    ##############
    ## Place ##
    #############
    parser_place = subparsers.add_parser(
        'place', help='Place a sample on a prebuilt tree', parents=[
            db_parser_mixin, force_mixin])
    parser_place.add_argument(
        'sample',
        metavar='sample',
        type=str,
        help='sample id')
    parser_place.add_argument('--tree', metavar='tree', type=str, help='tree')
    parser_place.add_argument(
        '--searchable_samples',
        metavar='searchable_samples',
        type=str,
        help='list of samples (file)')
    parser_place.add_argument(
        '--no-cache',
        default=False,
        action="store_true")
    parser_place.set_defaults(func=run_subtool)

    ##

    ##############
    ## Place ##
    #############
    parser_diff = subparsers.add_parser(
        'diff',
        help='Outputs novel sequence by calculating the difference between the sequence and combined graph',
        parents=[sequence_or_binary_parser_mixin])
    parser_diff.add_argument(
        'graph',
        metavar='graph',
        type=str,
        help='The graph to compare new sample against')
    parser_diff.add_argument(
        '--add',
        default=False,
        action="store_true",
        help="after comparing, add the new sample to the graph")
    parser_diff.set_defaults(func=run_subtool)
    ##

    args = parser.parse_args()
    args.func(parser, args)


if __name__ == "__main__":
    main()
