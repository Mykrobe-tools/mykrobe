from __future__ import print_function
import logging
logger = logging.getLogger(__name__)

from pprint import pprint
import json
import os
from mykrobe.utils import check_args
from mykrobe.typing import CoverageParser
from mykrobe.typing import Genotyper
from mykrobe.predict import TBPredictor
from mykrobe.predict import StaphPredictor
from mykrobe.predict import GramNegPredictor
from mykrobe.predict import MykrobePredictorSusceptibilityResult
from mykrobe.metagenomics import AMRSpeciesPredictor
from mykrobe.metagenomics import MykrobePredictorPhylogeneticsResult
from mykrobe.version import __version__ as predictor_version
from mykrobe.version import __version__ as atlas_version

from mongoengine import EmbeddedDocumentField
from mongoengine import IntField
from mongoengine import DictField
from mongoengine import StringField

STAPH_PANELS = ["data/panels/staph-species-160227.fasta.gz",
                "data/panels/staph-amr-bradley_2015-feb-17-2017.fasta.gz"]

GN_PANELS = [
    "data/panels/gn-amr-genes",
    "data/panels/Escherichia_coli",
    "data/panels/Klebsiella_pneumoniae",
    "data/panels/gn-amr-genes-extended"]


class MykrobePredictorResult(object):

    def __init__(self, susceptibility, phylogenetics,
                 variant_calls, sequence_calls,
                 kmer, probe_sets, files, version, model):
        self.susceptibility = susceptibility
        self.phylogenetics = phylogenetics
        self.variant_calls = variant_calls
        self.sequence_calls = sequence_calls
        self.kmer = kmer
        self.probe_sets = probe_sets
        self.files = files
        self.version = version
        self.model = model

    def to_dict(self):
        return {"susceptibility": list(self.susceptibility.to_dict().values())[0],
                "phylogenetics": list(self.phylogenetics.to_dict().values())[0],
                "variant_calls": self.variant_calls,
                "sequence_calls": self.sequence_calls,
                "kmer": self.kmer,
                "probe_sets": self.probe_sets,
                "files": self.files,
                "version": self.version,
                "genotype_model": self.model
                }
    # For database document
    # susceptibility = EmbeddedDocumentField("MykrobePredictorSusceptibilityResult")
    # phylogenetics = EmbeddedDocumentField("MykrobePredictorPhylogeneticsResult")
    # kmer = IntField()
    # probe_sets = StringField()
    # variant_calls = DictField()
    # sequence_calls = DictField()


def run(parser, args):
    base_json = {args.sample: {}}
    args = parser.parse_args()
    hierarchy_json_file = None
    if args.panel is not None:
        if args.panel == "bradley-2015":
            TB_PANELS = [
                "data/panels/tb-species-170421.fasta.gz",
                "data/panels/tb-bradley-probe-set-feb-09-2017.fasta.gz"]
        elif args.panel == "walker-2015":
            TB_PANELS = [
                "data/panels/tb-species-170421.fasta.gz",
                "data/panels/tb-walker-probe-set-feb-09-2017.fasta.gz"]
    Predictor = None
    if not args.species:
        panels = TB_PANELS + GN_PANELS + STAPH_PANELS
        panel_name = "tb-gn-staph-amr"
    elif args.species == "staph":
        panels = STAPH_PANELS
        panel_name = "staph-amr"
        Predictor = StaphPredictor
        args.kmer = 15  # Forced
    elif args.species == "tb":
        panels = TB_PANELS
        panel_name = "tb-amr"
        hierarchy_json_file = "data/phylo/mtbc_hierarchy.json"
        Predictor = TBPredictor
    elif args.species == "gn":
        panels = GN_PANELS
        panel_name = "gn-amr"
        Predictor = GramNegPredictor
    logger.info("Running AMR prediction with panels %s" % ", ".join(panels))
    version = {}
    version["mykrobe-predictor"] = predictor_version
    version["mykrobe-atlas"] = atlas_version
    # Get real paths for panels
    panels = [
        os.path.realpath(
            os.path.join(
                os.path.dirname(__file__),
                "..",
                f)) for f in panels]
    if hierarchy_json_file is not None:
        hierarchy_json_file = os.path.realpath(
            os.path.join(
                os.path.dirname(__file__),
                "..",
                hierarchy_json_file))
    if args.ont:
        args.expected_error_rate = 0.15
        logger.debug("Setting expected error rate to %s (--ont)" %
                     args.expected_error_rate)
        args.filters = ["LOW_GT_CONF"]
        args.model = "kmer_count"
    # Run Cortex
    cp = CoverageParser(
        sample=args.sample,
        panel_file_paths=panels,
        seq=args.seq,
        kmer=args.kmer,
        force=args.force,
        threads=1,
        verbose=False,
        tmp_dir=args.tmp,
        skeleton_dir=args.skeleton_dir,
        mccortex31_path=args.mccortex31_path)
    cp.run()
    logger.debug('CoverageParser complete')

    # Detect species
    species_predictor = AMRSpeciesPredictor(
        phylo_group_covgs=cp.covgs.get(
            "complex",
            cp.covgs.get(
                "phylo_group",
                {})),
        sub_complex_covgs=cp.covgs.get(
            "sub-complex",
            {}),
        species_covgs=cp.covgs["species"],
        lineage_covgs=cp.covgs.get(
            "sub-species",
            {}),
        hierarchy_json_file=hierarchy_json_file)
    phylogenetics = species_predictor.run()

    # ## AMR prediction

    depths = []
    if species_predictor.is_saureus_present():
        depths = [species_predictor.out_json["phylogenetics"]
                  ["phylo_group"]["Staphaureus"]["median_depth"]]
    elif species_predictor.is_mtbc_present():
        depths = [species_predictor.out_json["phylogenetics"]["phylo_group"][
            "Mycobacterium_tuberculosis_complex"]["median_depth"]]
    elif species_predictor.is_gram_neg_present():
        Predictor = GramNegPredictor
        try:
            depths = [species_predictor.out_json["phylogenetics"][
                "species"]["Klebsiella_pneumoniae"]["median_depth"]]
        except KeyError:
            depths = [species_predictor.out_json["phylogenetics"]
                      ["species"]["Escherichia_coli"]["median_depth"]]
    # pprint (species_predictor.out_json["phylogenetics"]["species"])
    # Genotype
    q = args.quiet
    args.quiet = True
    variant_calls_dict = {}
    sequence_calls_dict = {}
    if args.force and not depths:
        depths = [1]
    gt = None
    if depths or args.force:
        gt = Genotyper(sample=args.sample,
                       expected_depths=depths,
                       expected_error_rate=args.expected_error_rate,
                       variant_covgs=cp.variant_covgs,
                       gene_presence_covgs=cp.covgs["presence"],
                       base_json=base_json,
                       contamination_depths=[],
                       report_all_calls=True,
                       ignore_filtered=True,
                       filters=args.filters,
                       variant_confidence_threshold=args.min_variant_conf,
                       sequence_confidence_threshold=args.min_gene_conf,
                       model=args.model
                       )
        gt.run()
        variant_calls_dict = gt.variant_calls_dict
        sequence_calls_dict = gt.sequence_calls_dict
    else:
        depths = [cp.estimate_depth()]
    args.quiet = q
    mykrobe_predictor_susceptibility_result = MykrobePredictorSusceptibilityResult()
    if gt is not None and (max(depths) > args.min_depth or args.force):
        predictor = Predictor(variant_calls=gt.variant_calls,
                              called_genes=gt.sequence_calls_dict,
                              base_json=base_json[args.sample],
                              depth_threshold=args.min_depth,
                              ignore_filtered=True,
                              ignore_minor_calls=args.ont)
        mykrobe_predictor_susceptibility_result = predictor.run()
    base_json[
        args.sample] = MykrobePredictorResult(
        susceptibility=mykrobe_predictor_susceptibility_result,
        phylogenetics=phylogenetics,
        variant_calls=variant_calls_dict,
        sequence_calls=sequence_calls_dict,
        probe_sets=panels,
        files=args.seq,
        kmer=args.kmer,
        version=version,
        model=args.model).to_dict()
    if not args.keep_tmp:
        cp.remove_temporary_files()

    # write to file is specified by user, otherwise send to stdout
    if args.output:
        with open(args.output, 'w') as outfile:
            json.dump(obj=base_json, fp=outfile, indent=4)
    else:
        print(json.dumps(base_json, indent=4))
