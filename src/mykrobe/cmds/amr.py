from __future__ import print_function

import logging

logger = logging.getLogger(__name__)

import json

import numpy as np
import os
import random
import time
from enum import Enum
from typing import NamedTuple, List, Union, Optional
from mykrobe.mformat import json_to_csv
from mykrobe.typing import CoverageParser
from mykrobe.typing import Genotyper
from mykrobe.typing.models.base import ProbeCoverage
from mykrobe.typing.models.variant import VariantProbeCoverage
from mykrobe.typing.typer.variant import VariantTyper
from mykrobe.predict import TBPredictor
from mykrobe.predict import StaphPredictor
from mykrobe.predict import MykrobePredictorSusceptibilityResult
from mykrobe.metagenomics import AMRSpeciesPredictor
from mykrobe.version import __version__ as predictor_version
from mykrobe.version import __version__ as atlas_version


PathLike = Union[str, os.PathLike]

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
        # f = open('test.simulate_snps_data.tsv', 'w')
        # print('Correct_cov', 'Incorrect_cov', 'correct_k_count', 'incorrect_k_count', 'correct_percent_coverage', 'incorrect_percent_coverage', 'Cov', 'Conf', sep='\t', file=f)
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

            # correct_percent_coverage = ConfThresholder._simulate_percent_coverage(int(correct_k_count), 2 + self.kmer_length)
            # incorrect_percent_coverage = ConfThresholder._simulate_percent_coverage(int(incorrect_k_count), 2 + self.kmer_length)
            correct_percent_coverage = 100
            incorrect_percent_coverage = ConfThresholder._get_incorrect_kmer_percent_cov(
                int(incorrect_k_count), self.incorrect_kmer_to_pc_cov
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
            # print(correct_covg[i], incorrect_covg[i], correct_k_count, incorrect_k_count, correct_percent_coverage, incorrect_percent_coverage, cov, conf, sep='\t', file=f)

        # f.close()
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


class MykrobePredictorResult(object):
    def __init__(
        self,
        susceptibility,
        phylogenetics,
        variant_calls,
        sequence_calls,
        kmer,
        probe_sets,
        files,
        version,
        model,
    ):
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
        return {
            "susceptibility": list(self.susceptibility.to_dict().values())[0],
            "phylogenetics": list(self.phylogenetics.to_dict().values())[0],
            "variant_calls": self.variant_calls,
            "sequence_calls": self.sequence_calls,
            "kmer": self.kmer,
            "probe_sets": self.probe_sets,
            "files": self.files,
            "version": self.version,
            "genotype_model": self.model,
        }

    # For database document
    # susceptibility = EmbeddedDocumentField("MykrobePredictorSusceptibilityResult")
    # phylogenetics = EmbeddedDocumentField("MykrobePredictorPhylogeneticsResult")
    # kmer = IntField()
    # probe_sets = StringField()
    # variant_calls = DictField()
    # sequence_calls = DictField()


class Species(Enum):
    TB = "tb"
    STAPH = "staph"
    GN = "gn"


class TbPanel(Enum):
    BRADLEY = "bradley-2015"
    WALKER = "walker-2015"
    NEJM_WALKER = "201901"
    ATLAS = "atlas"
    CUSTOM = "custom"


class StaphPanel(Enum):
    DEFAULT = "default"
    CUSTOM = "custom"


class GnPanel(Enum):
    DEFAULT = "default"


PanelName = Union[StaphPanel, TbPanel, GnPanel]


class Panel(NamedTuple):
    paths: List[str]
    name: PanelName

    @staticmethod
    def from_species_and_name(species: Species, name: str) -> "Panel":
        if species is Species.STAPH:
            panel_name = StaphPanel(name)
        elif species is Species.TB:
            panel_name = TbPanel(name)
        elif species is Species.GN:
            panel_name = GnPanel(name)
        else:
            raise NameError(f"{species} is not a known species.")

        paths: List[str] = PANELS[species][panel_name]
        return Panel(paths, panel_name)

    def add_path(self, *args):
        for path in args:
            self.paths.append(path)


PANELS = {
    Species.STAPH: {
        StaphPanel.DEFAULT: [
            "data/panels/staph-species-160227.fasta.gz",
            "data/panels/staph-amr-bradley_2015-feb-17-2017.fasta.gz",
        ],
        StaphPanel.CUSTOM: ["data/panels/staph-species-160227.fasta.gz"],
    },
    Species.GN: {
        GnPanel.DEFAULT: [
            "data/panels/gn-amr-genes",
            "data/panels/Escherichia_coli",
            "data/panels/Klebsiella_pneumoniae",
            "data/panels/gn-amr-genes-extended",
        ]
    },
    Species.TB: {
        TbPanel.ATLAS: [
            "data/panels/tb-species-170421.fasta.gz",
            "data/panels/tb-walker-probe-set-jan-2019.fasta.gz",
            "data/panels/tb-k21-probe-set-feb-09-2017.fasta.gz",
        ],
        TbPanel.BRADLEY: [
            "data/panels/tb-species-170421.fasta.gz",
            "data/panels/tb-bradley-probe-set-jan-2019.fasta.gz",
        ],
        TbPanel.WALKER: [
            "data/panels/tb-species-170421.fasta.gz",
            "data/panels/tb-walker-probe-set-jan-2019.fasta.gz",
        ],
        TbPanel.NEJM_WALKER: [
            "data/panels/tb-species-170421.fasta.gz",
            "data/panels/tb-hunt-probe-set-jan-03-2019.fasta.gz",
        ],
        TbPanel.CUSTOM: ["data/panels/tb-species-170421.fasta.gz"],
    },
}


def describe_panels(parser, args):
    all_panels = {
        "tb": {
            "bradley-2015": {
                "description": "AMR panel described in Bradley, P et al. Rapid antibiotic-resistance predictions from genome sequence data for Staphylococcus aureus and Mycobacterium tuberculosis. Nat. Commun. 6:10063 doi: 10.1038/ncomms10063 (2015).",
                "reference": "NC_000962.3",
            },
            "walker-2015": {
                "description": "AMR panel described in Walker, Timothy M et al. Whole-genome sequencing for prediction of Mycobacterium tuberculosis drug susceptibility and resistance: a retrospective cohort study. The Lancet Infectious Diseases , Volume 15 , Issue 10 , 1193 - 1202",
                "reference": "NC_000962.3",
            },
            "201901": {
                "description": "AMR panel based on first line drugs from NEJM-2018 variants (DOI 10.1056/NEJMoa1800474), and second line drugs from Walker 2015 panel.",
                "reference": "NC_000962.3",
            },
        },
        "staph": {
            "default": {
                "description": "AMR panel described in Bradley, P et al. Rapid antibiotic-resistance predictions from genome sequence data for Staphylococcus aureus and Mycobacterium tuberculosis. Nat. Commun. 6:10063 doi: 10.1038/ncomms10063 (2015).",
                "reference": "BX571856.1",
            }
        },
    }

    print("\t".join(["species", "panel-name", "description", "reference"]))
    for species, panels in all_panels.items():
        for panel_name, data in panels.items():
            print(
                "\t".join([species, panel_name, data["description"], data["reference"]])
            )


def run(parser, args):
    base_json = {args.sample: {}}
    args = parser.parse_args()
    hierarchy_json_file = None
    variant_to_resistance_json_fp: Optional[PathLike] = None
    species = Species(args.species)
    if species is not Species.TB and args.panel != "custom":
        args.panel = "default"
    panels = Panel.from_species_and_name(species, args.panel)

    if species is Species.TB and panels.name is TbPanel.NEJM_WALKER:
        data_dir = os.path.abspath(
            os.path.join(os.path.dirname(__file__), "../data/predict/tb/")
        )
        variant_to_resistance_json_fp = os.path.join(
            data_dir, "variant_to_resistance_drug-jan-03-2019.json"
        )
    if panels.name in (TbPanel.CUSTOM, StaphPanel.CUSTOM):
        if not args.custom_probe_set_path:
            raise ValueError("Custom panel requires custom_probe_set_path")

        if not os.path.exists(args.custom_probe_set_path):
            raise FileNotFoundError(
                f"Custom probe path {args.custom_probe_set_path} does not exist!"
            )
        panels.add_path(args.custom_probe_set_path)

        if not os.path.exists(args.custom_variant_to_resistance_json):
            raise FileNotFoundError(
                (
                    "Custom variant to resistance json "
                    f"{args.custom_variant_to_resistance_json} does not exist!"
                )
            )
        variant_to_resistance_json_fp = args.custom_variant_to_resistance_json

    # todo: the following two lines are flagged for deletion. Based on the current CLI
    # implementation, we cannot get to here without there being a species
    # if not species:
    #     panels = TB_PANELS + GN_PANELS + STAPH_PANELS
    if species is Species.STAPH:
        Predictor = StaphPredictor
        args.kmer = 15  # Forced
    elif species is Species.TB:
        hierarchy_json_file = "data/phylo/mtbc_hierarchy.json"
        Predictor = TBPredictor
    elif species is Species.GN:
        raise NotImplementedError("Predictor not implement for gram negatives.")
    else:
        raise ValueError(f"Unrecognised species {species}")

    logger.info("Running AMR prediction with panels %s" % ", ".join(panels.paths))
    version = dict()
    version["mykrobe-predictor"] = predictor_version
    version["mykrobe-atlas"] = atlas_version
    # Get real paths for panels
    panels = [
        os.path.realpath(os.path.join(os.path.dirname(__file__), "..", f))
        for f in panels.paths
    ]
    if hierarchy_json_file is not None:
        hierarchy_json_file = os.path.realpath(
            os.path.join(os.path.dirname(__file__), "..", hierarchy_json_file)
        )
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
    )
    cp.run()
    logger.debug("CoverageParser complete")

    # Detect species
    species_predictor = AMRSpeciesPredictor(
        phylo_group_covgs=cp.covgs.get("complex", cp.covgs.get("phylo_group", {})),
        sub_complex_covgs=cp.covgs.get("sub-complex", {}),
        species_covgs=cp.covgs["species"],
        lineage_covgs=cp.covgs.get("sub-species", {}),
        hierarchy_json_file=hierarchy_json_file,
    )
    phylogenetics = species_predictor.run()

    # ## AMR prediction

    depths = []
    if species_predictor.is_saureus_present():
        depths = [
            species_predictor.out_json["phylogenetics"]["phylo_group"]["Staphaureus"][
                "median_depth"
            ]
        ]
    elif species_predictor.is_mtbc_present():
        depths = [
            species_predictor.out_json["phylogenetics"]["phylo_group"][
                "Mycobacterium_tuberculosis_complex"
            ]["median_depth"]
        ]
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
            kmer_size=args.kmer,
            min_proportion_expected_depth=args.min_proportion_expected_depth,
            ploidy=args.ploidy,
        )
        gt.run()
        (
            kmer_count_error_rate,
            incorrect_kmer_to_pc_cov,
        ) = gt.estimate_kmer_count_error_rate_and_incorrect_kmer_to_percent_cov()
        logger.info(
            "Estimated error rate for kmer count model: "
            + str(round(100 * kmer_count_error_rate, 2))
            + "%"
        )
        if args.guess_sequence_method and kmer_count_error_rate > 0.001:
            logger.warning("Guess sequence method is on, and we've guessed ONT")
            args.ont = True

        if args.ont:
            args.expected_error_rate = 0.15
            args.ploidy = "haploid"
            args.ignore_minor_calls = True
            logger.warning("Setting ploidy to haploid")
            logger.warning("Setting ignore_minor_calls to True")
            logger.warning(
                "Setting expected error rate to %s (--ont)" % args.expected_error_rate
            )
            args.model = "kmer_count"

        # If the user didn't specify the conf_percent_cutoff, then set it
        # depending on whether or not the --ont flag was used
        if args.conf_percent_cutoff == -1:
            args.conf_percent_cutoff = 90 if args.ont else 100

        # conf_percent_cutoff == 100 means that we want to keep all variant calls,
        # in which case there is no need to run the simulations
        if args.conf_percent_cutoff < 100:
            logger.info("Expected depth: " + str(depths[0]))
            conf_thresholder = ConfThresholder(
                kmer_count_error_rate, depths[0], args.kmer, incorrect_kmer_to_pc_cov
            )
            time_start = time.time()
            conf_threshold = conf_thresholder.get_conf_threshold(
                percent_to_keep=args.conf_percent_cutoff
            )
            time_end = time.time()
            time_to_sim = time_end - time_start
            logger.info("Simulation time: " + str(time_to_sim))
            logger.info(
                "Confidence cutoff (using percent cutoff "
                + str(args.conf_percent_cutoff)
                + "%): "
                + str(conf_threshold)
            )
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
                kmer_size=args.kmer,
                min_proportion_expected_depth=args.min_proportion_expected_depth,
                ploidy=args.ploidy,
            )
            gt.run()

        variant_calls_dict = gt.variant_calls_dict
        sequence_calls_dict = gt.sequence_calls_dict
    else:
        depths = [cp.estimate_depth()]
    args.quiet = q
    mykrobe_predictor_susceptibility_result = MykrobePredictorSusceptibilityResult()
    if gt is not None and (max(depths) > args.min_depth or args.force):
        predictor = Predictor(
            variant_calls=gt.variant_calls,
            called_genes=gt.sequence_calls_dict,
            base_json=base_json[args.sample],
            depth_threshold=args.min_depth,
            ignore_filtered=True,
            ignore_minor_calls=args.ignore_minor_calls,
            variant_to_resistance_json_fp=variant_to_resistance_json_fp,
        )
        mykrobe_predictor_susceptibility_result = predictor.run()
    base_json[args.sample] = MykrobePredictorResult(
        susceptibility=mykrobe_predictor_susceptibility_result,
        phylogenetics=phylogenetics,
        variant_calls=variant_calls_dict,
        sequence_calls=sequence_calls_dict,
        probe_sets=panels,
        files=args.seq,
        kmer=args.kmer,
        version=version,
        model=args.model,
    ).to_dict()
    if not args.keep_tmp:
        cp.remove_temporary_files()

    outputs = {}

    if args.output_format in ["csv", "json_and_csv"]:
        outputs["csv"] = json_to_csv(base_json)
    if args.output_format in ["json", "json_and_csv"]:
        # Verbose json output requires --report_all_calls
        if not args.report_all_calls:
            del base_json[args.sample]["variant_calls"]
            del base_json[args.sample]["sequence_calls"]
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
