# Read the kmer counts into a hash
import json
import logging
import tempfile

from mykrobe.typing import CoverageParser
from mykrobe.typing import Genotyper
from mykrobe.utils import load_json
from mykrobe.version import __version__
from mykrobe import ONT_E_RATE

logger = logging.getLogger(__name__)


def run_main(parser, args):
    args = parser.parse_args()
    verbose = True
    if args.ont:
        args.expected_error_rate = ONT_E_RATE
        logger.debug(
            "Setting expected error rate to %s (--ont)" % args.expected_error_rate
        )

    if args.min_variant_conf is None:
        args.min_variant_conf = 100

    if args.tmp is None:
        args.tmp = tempfile.mkdtemp() + "/"

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
    )
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
    if args.lineage is None:
        lineage_dict = None
    else:
        lineage_dict = load_json(args.lineage)
    gt = Genotyper(
        sample=args.sample,
        expected_error_rate=args.expected_error_rate,
        expected_depths=[args.expected_depth],
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
        min_proportion_expected_depth=args.min_proportion_expected_depth,
        ploidy=args.ploidy,
        lineage_variants=lineage_dict,
    )
    gt.run()
    if args.output:
        with open(args.output, "w") as outfile:
            json.dump(gt.out_json, outfile, indent=4)

    if not args.keep_tmp:
        cp.remove_temporary_files()
    return gt.out_json


def run(parser, args):
    print(json.dumps(run_main(parser, args), indent=1))
