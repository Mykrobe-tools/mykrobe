from __future__ import print_function

import copy
import json
import logging
import random
import sys
import tempfile
import time

import numpy as np

from mykrobe import ONT_E_RATE, ONT_PLOIDY
from mykrobe.metagenomics import AMRSpeciesPredictor
from mykrobe.mformat import json_to_csv
from mykrobe.predict import BasePredictor
from mykrobe.predict import MykrobePredictorSusceptibilityResult
from mykrobe.species_data import DataDir
from mykrobe.typing import CoverageParser
from mykrobe.typing import Genotyper
from mykrobe.typing.models.base import ProbeCoverage
from mykrobe.typing.models.variant import VariantProbeCoverage
from mykrobe.typing.typer.variant import VariantTyper
from mykrobe.utils import fix_amino_acid_X_variants_keys
from mykrobe.utils import load_json
from mykrobe.version import __version__ as atlas_version
from mykrobe.version import __version__ as predictor_version

logger = logging.getLogger(__name__)
random.seed(42)


class ConfThresholder:
    def __init__(
        self,
        error_rate,
        mean_depth,
        kmer_length,
        incorrect_kmer_to_pc_cov,
        iterations=10000,
    ):
        self.error_rate = error_rate
        self.mean_depth = mean_depth
        self.iterations = iterations
        self.incorrect_kmer_to_pc_cov = incorrect_kmer_to_pc_cov
        self.kmer_length = kmer_length

        # For now, store the coverage as well, in case we need to debug.
        # In future, could just store the confidences, as that's all we
        # need to decide on the cutoff
        self.log_conf_and_covg = []

    @classmethod
    def _simulate_percent_coverage(cls, kmers_to_sample, kmers_in_allele):
        """Simulates sampling kmers, returning the percent coverage.
        This means the % of kmers that are found by the simulation"""
        found_kmers = set()

        for _ in range(kmers_to_sample):
            j = random.randint(1, kmers_in_allele)
            found_kmers.add(j)
            if len(found_kmers) == kmers_in_allele:
                return 100.0

        return 100 * len(found_kmers) / kmers_in_allele

    @classmethod
    def _get_incorrect_kmer_percent_cov(cls, k_count, incorrect_kmer_to_pc_cov):
        i = k_count
        while i >= 0:
            if i in incorrect_kmer_to_pc_cov:
                return incorrect_kmer_to_pc_cov[i]
            i -= 1

        return 0

    def _simulate_snps(self):
        correct_covg = np.random.poisson(lam=self.mean_depth, size=self.iterations)
        incorrect_covg = np.random.binomial(
            self.mean_depth, self.error_rate, size=self.iterations
        )
        probe_coverage_list = []
        vtyper = VariantTyper(
            [self.mean_depth], error_rate=self.error_rate, kmer_size=self.kmer_length
        )

        for i in range(self.iterations):
            if correct_covg[i] + incorrect_covg[i] == 0:
                continue

            min_depth = 1  # not used?
            # Check what allele_length means in depth_to_expected_kmer_coun()! Probably need to change next two lines...
            correct_k_count = (
                self.kmer_length * correct_covg[i]
            ) + 0.01  #  see KmerCountGenotypeModel.depth_to_expected_kmer_count()
            incorrect_k_count = (
                self.kmer_length * incorrect_covg[i]
            ) + 0.01  #  see KmerCountGenotypeModel.depth_to_expected_kmer_count()

            correct_percent_coverage = 100
            incorrect_percent_coverage = (
                ConfThresholder._get_incorrect_kmer_percent_cov(
                    int(incorrect_k_count), self.incorrect_kmer_to_pc_cov
                )
            )
            correct_probe_coverage = ProbeCoverage(
                correct_percent_coverage,
                self.mean_depth,
                min_depth,
                correct_k_count,
                self.kmer_length,
            )
            incorrect_probe_coverage = ProbeCoverage(
                incorrect_percent_coverage,
                self.mean_depth,
                min_depth,
                incorrect_k_count,
                self.kmer_length,
            )
            vpc = VariantProbeCoverage(
                [correct_probe_coverage], [incorrect_probe_coverage]
            )
            call = vtyper.type(vpc)

            cov = np.log10(correct_covg[i] + incorrect_covg[i])
            conf = call["info"]["conf"]
            self.log_conf_and_covg.append((conf, cov))

        self.log_conf_and_covg.sort(reverse=True)

    def get_conf_threshold(self, percent_to_keep=95):
        """percent_to_keep determines the confidence cutoff in terms of
        simulatead confience scores. Should be in (0,100]. eg default of 95
        means that we choose a cutoff that would keep 95% of the data"""
        if len(self.log_conf_and_covg) == 0:
            self._simulate_snps()

        conf_cutoff_index = min(
            int(0.01 * percent_to_keep * len(self.log_conf_and_covg)),
            len(self.log_conf_and_covg) - 1,
        )
        return self.log_conf_and_covg[conf_cutoff_index][0]


def ref_data_from_args(args):
    if args.species == "custom":
        if args.custom_probe_set_path is None:
            raise ValueError(
                "Must use --custom_probe_set_path option if the species is 'custom'"
            )
        ref_data = {
            "fasta_files": [args.custom_probe_set_path],
            "var_to_res_json": args.custom_variant_to_resistance_json,
            "hierarchy_json": None,
            "lineage_json": args.custom_lineage_json,
            "kmer": args.kmer,
            "version": "custom",
            "species_phylo_group": None,
        }
    else:
        data_dir = DataDir(args.panels_dir)
        species_dir = data_dir.get_species_dir(args.species)
        if args.panel is not None:
            species_dir.set_panel(args.panel)
        ref_data = {
            "fasta_files": species_dir.fasta_files(),
            "var_to_res_json": species_dir.json_file("amr"),
            "hierarchy_json": species_dir.json_file("hierarchy"),
            "lineage_json": species_dir.json_file("lineage"),
            "kmer": species_dir.kmer(),
            "version": species_dir.version(),
            "species_phylo_group": species_dir.species_phylo_group(),
        }

    if ref_data["lineage_json"] is None:
        ref_data["lineage_dict"] = None
    else:
        ref_data["lineage_dict"] = load_json(ref_data["lineage_json"])

    return ref_data


def detect_species_and_get_depths(cov_parser, hierarchy_json, wanted_phylo_group):
    depths = []
    if wanted_phylo_group is None:
        return {}, depths

    species_predictor = AMRSpeciesPredictor(
        phylo_group_covgs=cov_parser.covgs.get(
            "complex", cov_parser.covgs.get("phylo_group", {})
        ),
        sub_complex_covgs=cov_parser.covgs.get("sub-complex", {}),
        species_covgs=cov_parser.covgs["species"],
        lineage_covgs=cov_parser.covgs.get("sub-species", {}),
        hierarchy_json_file=hierarchy_json,
    )
    phylogenetics = species_predictor.run()

    if wanted_phylo_group in species_predictor.out_json["phylogenetics"]["phylo_group"]:
        depths = [
            species_predictor.out_json["phylogenetics"]["phylo_group"][
                wanted_phylo_group
            ]["median_depth"]
        ]
    return phylogenetics, depths


def write_outputs(args, base_json):
    outputs = {}

    if args.output_format in ["csv", "json_and_csv"]:
        outputs["csv"] = json_to_csv(base_json)
    if args.output_format in ["json", "json_and_csv"]:
        # Verbose json output requires --report_all_calls
        if not args.report_all_calls:
            del base_json[args.sample]["variant_calls"]
            del base_json[args.sample]["sequence_calls"]
            del base_json[args.sample]["lineage_calls"]
        outputs["json"] = json.dumps(base_json, indent=4)

    if len(outputs) == 0:
        raise ValueError(
            (
                f"Output format must be one of: csv,json,json_and_csv. Got "
                f"'{args.output_format}'"
            )
        )

    for output_type, output in outputs.items():
        # write to file is specified by user, otherwise send to stdout
        if args.output:
            if args.output_format == "json_and_csv":
                outfile = args.output + "." + output_type
            else:
                outfile = args.output
            with open(outfile, "w") as f:
                f.write(output)
        else:
            print(output)


def fix_X_amino_acid_variants(sample_json):
    if "susceptibility" in sample_json:
        for drug_dict in sample_json["susceptibility"].values():
            if "called_by" in drug_dict:
                fix_amino_acid_X_variants_keys(drug_dict["called_by"])

    if "variant_calls" in sample_json:
        fix_amino_acid_X_variants_keys(sample_json["variant_calls"])


def run(parser, args):
    logger.info(f"Start runnning mykrobe predict. Command line: {' '.join(sys.argv)}")
    base_json = {args.sample: {}}
    ref_data = ref_data_from_args(args)
    if (
        args.species == "custom"
        and ref_data["var_to_res_json"] is None
        and ref_data["lineage_json"] is None
    ):
        logger.info(
            "Forcing --report_all_calls because species is 'custom' and options --custom_variant_to_resistance_json,--custom_lineage_json were not used"
        )
        args.report_all_calls = True
    logger.info(
        f"Running mykrobe predict using species {args.species}, and panel version {ref_data['version']}"
    )

    if args.tmp is None:
        args.tmp = tempfile.mkdtemp() + "/"

    if args.ont:
        args.expected_error_rate = ONT_E_RATE
        logger.info(
            f"Set expected error rate to {args.expected_error_rate} because --ont flag was used"
        )
        args.ploidy = ONT_PLOIDY
        logger.info(f"Set ploidy to {args.ploidy} because --ont flag used")

    # Run Cortex
    cp = CoverageParser(
        sample=args.sample,
        panel_file_paths=ref_data["fasta_files"],
        seq=args.seq,
        kmer=ref_data["kmer"],
        force=args.force,
        threads=args.threads,
        verbose=False,
        tmp_dir=args.tmp,
        skeleton_dir=args.skeleton_dir,
        memory=args.memory,
    )
    cp.run()
    logger.debug("CoverageParser complete")

    if ref_data["species_phylo_group"] is None:
        phylogenetics = {}
        depths = [cp.estimate_depth()]
    else:
        phylogenetics, depths = detect_species_and_get_depths(
            cp, ref_data["hierarchy_json"], ref_data["species_phylo_group"]
        )

    # Genotype
    variant_calls_dict = {}
    sequence_calls_dict = {}
    lineage_calls_dict = {}
    lineage_predict_dict = {}
    if args.force and len(depths) == 0:
        depths = [cp.estimate_depth()]
    gt = None

    if len(depths) > 0 or args.force:
        # Running the genotyper changes the contents of cp.covgs["presence"].
        # The changes can cause later second run (if it happens) of the
        # genotytper to crash.
        # Store the original cp.covgs["presence"] so we can use it again later.
        original_covgs_presence = copy.deepcopy(cp.covgs["presence"])
        gt = Genotyper(
            sample=args.sample,
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
            kmer_size=ref_data["kmer"],
            min_proportion_expected_depth=args.min_proportion_expected_depth,
            ploidy=args.ploidy,
            lineage_variants=ref_data["lineage_dict"],
        )
        gt.run()
        (
            kmer_count_error_rate,
            incorrect_kmer_to_pc_cov,
        ) = gt.estimate_kmer_count_error_rate_and_incorrect_kmer_to_percent_cov()
        logger.debug(
            "Estimated error rate for kmer count model: "
            + str(round(100 * kmer_count_error_rate, 2))
            + "%"
        )
        if args.guess_sequence_method and kmer_count_error_rate > 0.001:
            logger.warning("Guess sequence method is on, and we've guessed ONT")
            # this is duplicated from the if args.ont statement above
            args.expected_error_rate = ONT_E_RATE
            logger.info(f"Set expected error rate to {args.expected_error_rate}")
            args.ploidy = ONT_PLOIDY
            logger.info(f"Set ploidy to {args.ploidy}")

        # conf_percent_cutoff == 100 means that we want to keep all variant calls,
        # in which case there is no need to run the simulations
        if args.conf_percent_cutoff < 100:
            logger.debug("Expected depth: " + str(depths[0]))
            conf_thresholder = ConfThresholder(
                kmer_count_error_rate,
                depths[0],
                ref_data["kmer"],
                incorrect_kmer_to_pc_cov,
            )
            time_start = time.time()
            conf_threshold = conf_thresholder.get_conf_threshold(
                percent_to_keep=args.conf_percent_cutoff
            )
            time_end = time.time()
            time_to_sim = time_end - time_start
            logger.debug("Simulation time: " + str(time_to_sim))
            logger.debug(
                "Confidence cutoff (using percent cutoff "
                + str(args.conf_percent_cutoff)
                + "%): "
                + str(conf_threshold)
            )
            cp.covgs["presence"] = copy.deepcopy(original_covgs_presence)
            gt = Genotyper(
                sample=args.sample,
                expected_depths=depths,
                expected_error_rate=kmer_count_error_rate,
                variant_covgs=cp.variant_covgs,
                gene_presence_covgs=cp.covgs["presence"],
                base_json=base_json,
                contamination_depths=[],
                report_all_calls=True,
                ignore_filtered=True,
                filters=args.filters,
                variant_confidence_threshold=conf_threshold,
                sequence_confidence_threshold=args.min_gene_conf,
                model=args.model,
                kmer_size=ref_data["kmer"],
                min_proportion_expected_depth=args.min_proportion_expected_depth,
                ploidy=args.ploidy,
                lineage_variants=ref_data["lineage_dict"],
            )
            gt.run()

        variant_calls_dict = gt.variant_calls_dict
        sequence_calls_dict = gt.sequence_calls_dict
        lineage_predict_dict, lineage_calls_dict = gt.predict_lineage()
    else:
        depths = [cp.estimate_depth()]

    mykrobe_predictor_susceptibility_result = MykrobePredictorSusceptibilityResult()
    if (
        gt is not None
        and (max(depths) > args.min_depth or args.force)
        and ref_data["var_to_res_json"] is not None
    ):
        predictor = BasePredictor(
            variant_calls=gt.variant_calls,
            called_genes=gt.sequence_calls_dict,
            base_json=base_json[args.sample],
            depth_threshold=args.min_depth,
            ignore_filtered=True,
            ignore_minor_calls=args.ignore_minor_calls,
            variant_to_resistance_json_fp=ref_data["var_to_res_json"],
        )
        mykrobe_predictor_susceptibility_result = predictor.run()
        logger.info("Progress: finished making AMR predictions")

    base_json[args.sample] = {
        "susceptibility": list(
            mykrobe_predictor_susceptibility_result.to_dict().values()
        )[0],
        "phylogenetics": {}
        if phylogenetics == {}
        else list(phylogenetics.to_dict().values())[0],
        "variant_calls": variant_calls_dict,
        "sequence_calls": sequence_calls_dict,
        "lineage_calls": lineage_calls_dict,
        "kmer": ref_data["kmer"],
        "probe_sets": ref_data["fasta_files"],
        "files": args.seq,
        "version": {
            "mykrobe-predictor": predictor_version,
            "mykrobe-atlas": atlas_version,
            "panel": ref_data["version"],
        },
        "genotype_model": args.model,
    }
    if len(lineage_predict_dict) > 0:
        base_json[args.sample]["phylogenetics"]["lineage"] = lineage_predict_dict

    if not args.keep_tmp:
        cp.remove_temporary_files()

    logger.info("Progress: writing output")
    fix_X_amino_acid_variants(base_json[args.sample])
    write_outputs(args, base_json)
    logger.info("Progress: finished")
