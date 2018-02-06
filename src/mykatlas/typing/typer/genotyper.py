from __future__ import print_function
import os
import json
import csv
import glob
import logging
import subprocess
from copy import copy

from mykatlas.typing import SequenceProbeCoverage
from mykatlas.typing import VariantProbeCoverage
from mykatlas.typing import ProbeCoverage
from mykatlas.typing import Panel

from mykatlas.typing.typer.presence import GeneCollectionTyper
from mykatlas.typing.typer.variant import VariantTyper
from mykatlas.typing.typer.base import DEFAULT_MINOR_FREQ
from mykatlas.typing.typer.base import DEFAULT_ERROR_RATE

from ga4ghmongo.schema import VariantCallSet
from ga4ghmongo.schema import Variant

from mykatlas.cortex import McCortexGenoRunner

from mykatlas.utils import get_params
from mykatlas.utils import split_var_name
from mykatlas.utils import median

logger = logging.getLogger(__name__)


class CoverageParser(object):

    def __init__(
            self,
            sample,
            panel_file_paths,
            kmer,
            force,
            seq=None,
            ctx=None,
            threads=2,
            memory="1GB",
            panels=None,
            verbose=True,
            tmp_dir='tmp/',
            skeleton_dir='atlas/data/skeletons/',
            mccortex31_path="mccortex31"):
        self.sample = sample
        self.seq = seq
        self.ctx = ctx
        self.kmer = kmer
        self.force = force
        self.covgs = {"variant": {}, "presence": {}}
        self.variant_covgs = self.covgs["variant"]
        self.gene_presence_covgs = self.covgs["presence"]
        self.mc_cortex_runner = None
        self.verbose = verbose
        self.skeleton_dir = skeleton_dir
        self.tmp_dir = tmp_dir
        self.panel_file_paths = panel_file_paths
        self.panels = []
        self.mccortex31_path = mccortex31_path
        self.threads = threads
        self.memory = memory
        for panel_file_path in self.panel_file_paths:
            panel = Panel(panel_file_path)
            self.panels.append(panel)
        if self.seq and self.ctx:
            raise ValueError("Can't have both -1 and -c")

    def run(self):
        self._run_cortex()
        self._parse_covgs()

    def _run_cortex(self):
        self.mc_cortex_runner = McCortexGenoRunner(
            sample=self.sample,
            panels=self.panels,
            seq=self.seq,
            ctx=self.ctx,
            kmer=self.kmer,
            force=self.force,
            threads=self.threads,
            memory=self.memory,
            panel_name=self.panel_name,
            tmp_dir=self.tmp_dir,
            skeleton_dir=self.skeleton_dir,
            mccortex31_path=self.mccortex31_path)
        self.mc_cortex_runner.run()

    def estimate_depth(self):
        depth = []
        for variant_coverages in self.variant_covgs.values():
            for variant_covg in variant_coverages:
                if variant_covg.reference_coverage.median_depth > 0:
                    depth.append(variant_covg.reference_coverage.median_depth)
        for spcs in self.gene_presence_covgs.values():
            __median_depth = median(
                [spc.median_depth for spc in spcs.values()])
            if __median_depth > 0:
                depth.append(__median_depth)
        _median = median(depth)
        if _median < 1:
            return 1
        else:
            return _median

    def remove_temporary_files(self):
        self.mc_cortex_runner.remove_temporary_files()

    @property
    def panel_name(self):
        return "-".join([panel.name for panel in self.panels])

    def _parse_summary_covgs_row(self, row):
        try:
            return row[0], int(row[2]), int(row[3]), 100 * float(row[4]), int(row[5])
        except ValueError:
            logger.warning("Failed to parse %s" % str(row))
            return row[0], 0, 0, 0.0, 0

    def _parse_covgs(self):
        with open(self.mc_cortex_runner.covg_tmp_file_path, 'r') as infile:
            self.reader = csv.reader(infile, delimiter="\t")
            for row in self.reader:
                allele, median_depth, min_depth, percent_coverage, k_count = self._parse_summary_covgs_row(
                    row)
                allele_name = allele.split('?')[0]
                if self._is_variant_panel(allele_name):
                    self._parse_variant_panel(row)
                else:
                    self._parse_seq_panel(row)

    def _is_variant_panel(self, allele_name):
        try:
            alt_or_ref = allele_name.split('-')[0]
            return alt_or_ref in ["ref", "alt"]
        except ValueError:
            return False

    def _parse_seq_panel(self, row):
        allele, median_depth, min_depth, percent_coverage, k_count = self._parse_summary_covgs_row(
            row)
        probe_coverage = ProbeCoverage(
            percent_coverage=percent_coverage,
            median_depth=median_depth,
            min_depth=min_depth,
            k_count=k_count)

        allele_name = allele.split('?')[0]
        params = get_params(allele)
        panel_type = params.get("panel_type", "presence")
        name = params.get('name')
        version = params.get(
            'version',
            '1')
        if panel_type in ["variant", "presence"]:
            sequence_probe_coverage = SequenceProbeCoverage(
                name=name,
                probe_coverage=probe_coverage,
                version=version,
                length=params.get("length"))
            try:
                self.covgs[panel_type][name][version] = sequence_probe_coverage
            except KeyError:
                self.covgs[panel_type][name] = {}
                self.covgs[panel_type][name][version] = sequence_probe_coverage

        else:
            # Species panels are treated differently
            l = int(params.get("length", -1))
            try:
                self.covgs[panel_type][name]["total_bases"] += l
                if percent_coverage > 75 and median_depth > 0:
                    self.covgs[panel_type][name][
                        "percent_coverage"].append(percent_coverage)
                    self.covgs[panel_type][name]["length"].append(l)
                    self.covgs[panel_type][name]["median"].append(median_depth)
            except KeyError:
                if panel_type not in self.covgs:
                    self.covgs[panel_type] = {}
                self.covgs[panel_type][name] = {}
                self.covgs[panel_type][name]["total_bases"] = l
                if percent_coverage > 75 and median_depth > 0:
                    self.covgs[panel_type][name][
                        "percent_coverage"] = [percent_coverage]
                    self.covgs[panel_type][name]["length"] = [l]
                    self.covgs[panel_type][name]["median"] = [median_depth]
                else:
                    self.covgs[panel_type][name]["percent_coverage"] = []
                    self.covgs[panel_type][name]["length"] = []
                    self.covgs[panel_type][name]["median"] = []

    def _parse_variant_panel(self, row):
        allele, median_depth, min_depth, percent_coverage, k_count = self._parse_summary_covgs_row(
            row)
        params = get_params(allele)
        if 'var_name' in params:
            var_name = params.get('var_name')
        else:
            var_name = allele.split('?')[0].split('-')[1]

        num_alts = int(params.get("num_alts", 0))
        reference_coverages = [ProbeCoverage(
            percent_coverage=percent_coverage,
            median_depth=median_depth,
            min_depth=min_depth,
            k_count=k_count)]
        alt_or_ref = 'ref'
        alternate_coverages = []
        for i in range(num_alts-1):
            row = next(self.reader)
            ref_allele, median_depth, min_depth, percent_coverage, k_count = self._parse_summary_covgs_row(
                row)
            if ref_allele.split('-')[0] != 'ref':
                logger.warning(
                    "Fewer ref alleles than alt alleles for %s" % ref_allele)
                alternate_coverages.append(ProbeCoverage(
                    min_depth=min_depth,
                    k_count=k_count,
                    percent_coverage=percent_coverage,
                    median_depth=median_depth))
                num_alts -= 1
                break

            assert ref_allele.split('-')[0] == 'ref'
            reference_coverages.append(ProbeCoverage(
                percent_coverage=percent_coverage,
                median_depth=median_depth,
                min_depth=min_depth,
                k_count=k_count))
        for i in range(num_alts):
            row = next(self.reader)
            alt_allele, median_depth, min_depth, percent_coverage, k_count = self._parse_summary_covgs_row(
                row)
            assert alt_allele.split('-')[0] == 'alt'
            alternate_coverages.append(
                ProbeCoverage(
                    min_depth=min_depth,
                    k_count=k_count,
                    percent_coverage=percent_coverage,
                    median_depth=median_depth))
        variant_probe_coverage = VariantProbeCoverage(
            reference_coverages=reference_coverages,
            alternate_coverages=alternate_coverages,
            var_name=var_name,
            params=params)
        try:
            self.variant_covgs[allele].append(variant_probe_coverage)
        except KeyError:
            self.variant_covgs[allele] = [variant_probe_coverage]


class Genotyper(object):

    """Takes output of mccortex coverages and types"""

    def __init__(
            self,
            sample,
            expected_depths,
            variant_covgs,
            gene_presence_covgs,
            contamination_depths=[],
            base_json={},
            report_all_calls=False,
            ignore_filtered=False,
            filters=[],
            expected_error_rate=DEFAULT_ERROR_RATE,
            minor_freq=DEFAULT_MINOR_FREQ,
            variant_confidence_threshold=1,
            sequence_confidence_threshold=0,
            min_gene_percent_covg_threshold=100,
            model="median_depth"):
        self.sample = sample
        self.variant_covgs = variant_covgs
        self.gene_presence_covgs = gene_presence_covgs
        self.out_json = base_json
        if expected_depths[0] < 1:
            self.expected_depths = [1]
        else:
            self.expected_depths = expected_depths
        self.contamination_depths = contamination_depths
        self.variant_calls = {}
        self.sequence_calls = {}
        self.variant_calls_dict = {}
        self.sequence_calls_dict = {}
        self.expected_error_rate = expected_error_rate
        self.report_all_calls = report_all_calls
        self.filters = filters
        self.model = model
        self.ignore_filtered = ignore_filtered
        self.minor_freq = minor_freq
        self.variant_confidence_threshold = variant_confidence_threshold
        self.sequence_confidence_threshold = sequence_confidence_threshold
        self.min_gene_percent_covg_threshold = min_gene_percent_covg_threshold

    def run(self):
        self._type()

    def _type(self):
        self._type_genes()
        self._type_variants()

    def _type_genes(self):
        gt = GeneCollectionTyper(
            expected_depths=self.expected_depths,
            contamination_depths=self.contamination_depths,
            confidence_threshold=self.sequence_confidence_threshold)
        for gene_name, gene_collection in self.gene_presence_covgs.items():
            self.gene_presence_covgs[gene_name] = gt.type(
                gene_collection,
                min_gene_percent_covg_threshold=self.min_gene_percent_covg_threshold)
            self.sequence_calls_dict[gene_name] = [
                gpc.to_mongo().to_dict() for gpc in self.gene_presence_covgs[gene_name]]
        self.out_json[self.sample][
            "sequence_calls"] = self.sequence_calls_dict

    def _type_variants(self):
        self.out_json[self.sample]["variant_calls"] = {}
        gt = VariantTyper(
            expected_depths=self.expected_depths,
            error_rate=self.expected_error_rate,
            contamination_depths=self.contamination_depths,
            ignore_filtered=self.ignore_filtered,
            minor_freq=self.minor_freq,
            confidence_threshold=self.variant_confidence_threshold,
            filters=self.filters,
            model=self.model
        )
        genotypes = []
        filters = []
        for probe_name, probe_coverages in self.variant_covgs.items():
            probe_id = self._name_to_id(probe_name)
            variant = None

            call = gt.type(probe_coverages, variant=probe_name)

            genotypes.append(sum(call["genotype"]))
            filters.append(int(call["info"]["filter"] == "PASS"))
            if sum(call["genotype"]) > 0 or not call[
                    "genotype"] or self.report_all_calls:
                self.variant_calls[probe_name] = call
                self.variant_calls_dict[
                    probe_id] = call
        self.out_json[self.sample]["genotypes"] = genotypes
        self.out_json[self.sample]["filtered"] = filters
        self.out_json[self.sample]["variant_calls"] = self.variant_calls_dict

    def _name_to_id(self, probe_name):
        names = []
        params = get_params(probe_name)
        if params.get("mut"):
            names.append("_".join([params.get("gene"), params.get("mut")]))
            var_name = params.get("var_name")
        else:
            var_name = probe_name.split('?')[0].split('-')[1]
        names.append(var_name)
        return "-".join(names)

    def _create_variant(self, probe_name):
        names = []
        params = get_params(probe_name)
        if params.get("mut"):
            names.append("_".join([params.get("gene"), params.get("mut")]))
        var_name = probe_name.split('?')[0].split('-')[1]
        names.append(var_name)
        try:
            # If it's a variant panel we can create a variant
            ref, start, alt = split_var_name(var_name)
            return Variant.create(
                start=start,
                reference_bases=ref,
                alternate_bases=[alt],
                names=names,
                info=params)
        except AttributeError:
            return None
