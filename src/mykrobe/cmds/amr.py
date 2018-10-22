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


import numpy as np
from mykrobe.typing.models.base import ProbeCoverage
from mykrobe.typing.models.variant import VariantProbeCoverage
from mykrobe.typing.typer.variant import VariantTyper
import random
random.seed(42)

class ConfThresholder:
    def __init__(self, error_rate, mean_depth, kmer_length, iterations=10000):
        self.error_rate = error_rate
        self.mean_depth = mean_depth
        self.iterations = iterations
        self.kmer_length = kmer_length

        # For now, store the coverage as well, in case we need to debug.
        # In future, could just store the confidences, as that's all we
        # need to decide on the cutoff
        self.log_conf_and_covg = []


    @classmethod
    def _simulate_percent_coverage(cls, kmers_to_sample, kmers_in_allele):
        '''Simulates sampling kmers, returning the percent coverage.
        This means the % of kmers that are found by the simulation'''
        found_kmers = set()

        for _ in range(kmers_to_sample):
            j = random.randint(1, kmers_in_allele)
            found_kmers.add(j)
            if len(found_kmers) == kmers_in_allele:
                return 100.0

        return 100 * len(found_kmers) / kmers_in_allele


    def _simulate_snps(self):
        ref_covg = np.random.poisson(lam=self.mean_depth, size=self.iterations)
        alt_covg = np.random.binomial(self.mean_depth, self.error_rate, size=self.iterations)
        f = open('test.covs', 'w')
        print('Ref_cov', 'Alt_cov', 'Cov', 'Conf', sep='\t', file=f)
        probe_coverage_list = []
        vtyper = VariantTyper([self.mean_depth], error_rate=self.error_rate, kmer_size=self.kmer_length)

        for i in range(self.iterations):
            c1 = ref_covg[i]
            c2 = alt_covg[i]
            if c1 + c2 == 0:
                continue

            min_depth = 1 # not used?
            # Check what allele_length means in depth_to_expected_kmer_coun()! Probably need to change next two lines...
            ref_k_count = ((2 + self.kmer_length) * ref_covg[i]) + 0.01 # see KmerCountGenotypeModel.depth_to_expected_kmer_count()
            alt_k_count = ((2 + self.kmer_length) * alt_covg[i]) + 0.01 # see KmerCountGenotypeModel.depth_to_expected_kmer_count()
            logging.debug('ref_k_count ' + str(ref_k_count) + '. alt_k_count ' + str(alt_k_count))

            ref_percent_coverage = ConfThresholder._simulate_percent_coverage(int(ref_k_count), 2 + self.kmer_length)
            alt_percent_coverage = ConfThresholder._simulate_percent_coverage(int(alt_k_count), 2 + self.kmer_length)
            logging.debug('ref_percent_coverage ' + str(ref_percent_coverage) + '.  alt_percent_coverage ' + str(alt_percent_coverage))

            ref_probe_coverage = ProbeCoverage(ref_percent_coverage, self.mean_depth, min_depth, ref_k_count, self.kmer_length + 2)
            alt_probe_coverage = ProbeCoverage(alt_percent_coverage, self.mean_depth, min_depth, alt_k_count, self.kmer_length + 2)
            vpc = VariantProbeCoverage([ref_probe_coverage], [alt_probe_coverage])
            call = vtyper.type(vpc)

            cov = np.log10(c1 + c2)
            conf = call['info']['conf']
            self.log_conf_and_covg.append((conf, cov))
            print(c1, c2, cov, conf, sep='\t', file=f)

        f.close()
        self.log_conf_and_covg.sort(reverse=True)
        with open('test.log_conf_and_covg.tsv', 'w') as f:
            print('Conf\tCov', file=f)
            for t in self.log_conf_and_covg:
                print(*t, sep='\t', file=f)


    def get_conf_threshold(self, percent_to_keep=95):
        '''percent_to_keep determines the confidence cutoff in terms of
        simulatead confience scores. Should be in (0,100]. eg default of 95
        means that we choose a cutoff that would keep 95% of the data'''
        if len(self.log_conf_and_covg) == 0:
            self._simulate_snps()

        conf_cutoff_index = int(0.01 * percent_to_keep * len(self.log_conf_and_covg))
        return self.log_conf_and_covg[conf_cutoff_index][0]



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
        elif args.panel == "atlas":
            TB_PANELS = [
                "data/panels/tb-species-170421.fasta.gz",
                "data/panels/tb-walker-probe-set-feb-09-2017.fasta.gz",
                "data/panels/tb-k21-probe-set-feb-09-2017.fasta.gz"]
        elif args.panel == "custom":
            if not args.custom_probe_set_path:
                raise ValueError("Custom panel requires custom_probe_set_path")
            TB_PANELS = [
                args.custom_probe_set_path,
                "data/panels/tb-species-170421.fasta.gz"
            ]
    Predictor = None
    if not args.species:
        panels = TB_PANELS + GN_PANELS + STAPH_PANELS
    elif args.species == "staph":
        panels = STAPH_PANELS
        Predictor = StaphPredictor
        args.kmer = 15  # Forced
    elif args.species == "tb":
        panels = TB_PANELS
        hierarchy_json_file = "data/phylo/mtbc_hierarchy.json"
        Predictor = TBPredictor
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
        args.ploidy = "haploid"
        args.ignore_minor_calls = True
        logger.warning("Setting ploidy to haploid")
        logger.warning("Setting ignore_minor_calls to True")
        logger.warning("Setting expected error rate to %s (--ont)" %
                     args.expected_error_rate)
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
                       model=args.model,
                       kmer_size=args.kmer,
                       min_proportion_expected_depth=args.min_proportion_expected_depth,
                       ploidy=args.ploidy
                       )
        gt.run()
        variant_calls_dict = gt.variant_calls_dict
        sequence_calls_dict = gt.sequence_calls_dict
        kmer_count_error_rate = gt.estimate_kmer_count_error_rate()
        logger.info("Estimated error rate for kmer count model: " + str(round(100 * kmer_count_error_rate, 2)) + "%")
        logger.info("Expected depth: " + str(depths[0]))
        conf_thresholder = ConfThresholder(kmer_count_error_rate, depths[0], args.kmer)
        import time
        time_start = time.time()
        conf_threshold = conf_thresholder.get_conf_threshold()
        time_end = time.time()
        time_to_sim = time_end - time_start
        logger.info('Simulation time: ' + str(time_to_sim))
        logger.info("Confidence cutoff: " + str(conf_threshold))
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
                              ignore_minor_calls=args.ignore_minor_calls,
                              variant_to_resistance_json_fp=args.custom_variant_to_resistance_json)
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
            json.dump(base_json, outfile, indent=4)
    else:
        print(json.dumps(base_json, indent=4))
