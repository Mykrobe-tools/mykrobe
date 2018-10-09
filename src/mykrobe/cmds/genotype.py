# Read the kmer counts into a hash
from mykrobe.utils import check_args
from mykrobe.typing import Genotyper
from mykrobe.typing import CoverageParser
from mykrobe.version import __version__

from pprint import pprint
import json
import logging
logger = logging.getLogger(__name__)


def run_main(parser, args):
    args = parser.parse_args()
    verbose = True
    if args.ont:
        args.expected_error_rate = 0.15
        args.filters = ["LOW_GT_CONF"]
        args.model = "kmer_count"
        logger.debug("Setting expected error rate to %s (--ont)" %
                     args.expected_error_rate)
        logger.debug(
            "Removing LOW_PERCENT_COVERAGE filter (increases sensitivity - in particular for ONT data)")

    if args.min_variant_conf is None:
        args.min_variant_conf = 100
    cp = CoverageParser(
        sample=args.sample,
        panel_file_paths=[args.probe_set],
        seq=args.seq,
        ctx=args.ctx,
        kmer=args.kmer,
        force=args.force,
        verbose=verbose,
        tmp_dir=args.tmp,
        skeleton_dir=args.skeleton_dir,
        threads=args.threads,
        memory=args.memory,
        mccortex31_path=args.mccortex31_path)
    cp.run()
    if args.expected_depth is None:
        args.expected_depth = cp.estimate_depth()

    base_json = {args.sample: {}}
    base_json[args.sample]["probe_set"] = args.probe_set
    if args.seq:
        base_json[args.sample]["files"] = args.seq
    else:
        base_json[args.sample]["files"] = args.ctx
    base_json[args.sample]["kmer"] = args.kmer
    base_json[args.sample]["version"] = __version__
    gt = Genotyper(
        sample=args.sample,
        expected_error_rate=args.expected_error_rate,
        expected_depths=[
            args.expected_depth],
        variant_covgs=cp.variant_covgs,
        gene_presence_covgs=cp.covgs["presence"],
        base_json=base_json,
        contamination_depths=[],
        ignore_filtered=args.ignore_filtered,
        filters=args.filters,
        model=args.model,
        report_all_calls=args.report_all_calls,
        variant_confidence_threshold=args.min_variant_conf,
        sequence_confidence_threshold=args.min_gene_conf,
        min_gene_percent_covg_threshold=args.min_gene_percent_covg_threshold,
        kmer_size=args.kmer,
        min_proportion_expected_depth=args.min_proportion_expected_depth)
    gt.run()
    if args.output:
        with open(args.output, 'w') as outfile:
            json.dump(gt.out_json, outfile, indent=4)
                
    if not args.keep_tmp:
        cp.remove_temporary_files()
    return gt.out_json


def run(parser, args):
    print(json.dumps(run_main(parser, args), indent=1))
