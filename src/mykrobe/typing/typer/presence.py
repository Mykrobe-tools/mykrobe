from mykrobe.typing.typer.base import Typer
from mykrobe.typing.typer.base import MIN_LLK
from mykrobe.variants.schema.models import SequenceCall
from mykrobe.stats import log_lik_depth
from mykrobe.stats import percent_coverage_from_expected_coverage
from math import log


class PresenceTyper(Typer):

    "Initiated with expected depths and contamination depths"

    def __init__(self, expected_depths, contamination_depths=[],
                 confidence_threshold=1):
        super(
            PresenceTyper,
            self).__init__(
            expected_depths,
            contamination_depths,
            confidence_threshold=confidence_threshold)
        if len(expected_depths) > 1:
            raise NotImplementedError("Mixed samples not supported")

    def type(self, sequence_probe_coverage):
        "Takes a single SequenceCoverage object (or child) and returns genotype"
        call = self._type(sequence_probe_coverage)
        return call

    def _type(self, sequence_probe_coverage):
        hom_alt_likelihoods = []
        het_likelihoods = []
        hom_ref_likelihoods = []
        for expected_depth in self.expected_depths:
            hom_alt_likelihoods.append(
                self._hom_alt_likeihood(
                    median_depth=sequence_probe_coverage.median_depth,
                    expected_depth=expected_depth))
            if not self.has_contamination():
                het_likelihoods.append(
                    self._het_likelihood(
                        median_depth=sequence_probe_coverage.median_depth,
                        expected_depth=expected_depth))
            else:
                het_likelihoods.append(MIN_LLK)

            hom_ref_likelihoods.append(
                self._hom_ref_likelihood(
                    median_depth=sequence_probe_coverage.median_depth,
                    expected_depth=expected_depth))

            for contamination_depth in self.contamination_depths:
                hom_alt_likelihoods.append(
                    self._hom_alt_likeihood(
                        median_depth=sequence_probe_coverage.median_depth,
                        expected_depth=expected_depth +
                        contamination_depth))
                # NOTE : _HOM_ALT_LIKEIHOOD is not a typo
                hom_ref_likelihoods.append(
                    self._hom_alt_likeihood(
                        median_depth=sequence_probe_coverage.median_depth,
                        expected_depth=contamination_depth))
            # Posterior
        hom_ref_likelihood = self._log_post_hom_ref(max(hom_ref_likelihoods))
        hom_alt_likelihood = self._log_post_het_or_alt(
            max(hom_alt_likelihoods),
            expected_depth * 0.75,
            sequence_probe_coverage)
        het_likelihood = self._log_post_het_or_alt(
            max(het_likelihoods),
            expected_depth *
            self.minimum_detectable_frequency,
            sequence_probe_coverage)
        likelihoods = [hom_ref_likelihood, het_likelihood, hom_alt_likelihood]
        gt = self.likelihoods_to_genotype(likelihoods)
        info = {
            "copy_number": float(
                sequence_probe_coverage.median_depth) / expected_depth,
            "coverage": sequence_probe_coverage.coverage_dict,
            "expected_depths": self.expected_depths,
            "contamination_depths": self.contamination_depths
        }
        if sum([int(i) for i in gt.split("/")]) > 0:
            info["version"] = sequence_probe_coverage.version
        if sequence_probe_coverage.length is not None:
            info["length"] = sequence_probe_coverage.length

        return SequenceCall.create(
            sequence=None,
            call_set=None,
            genotype=gt,
            genotype_likelihoods=likelihoods,
            info=info)

    def _hom_alt_likeihood(self, median_depth, expected_depth):
        return log_lik_depth(median_depth, expected_depth * 0.75)

    def _het_likelihood(self, median_depth, expected_depth):
        return log_lik_depth(
            median_depth,
            expected_depth *
            self.minimum_detectable_frequency)

    def _hom_ref_likelihood(self, median_depth, expected_depth):
        return log_lik_depth(median_depth, expected_depth * 0.001)

    @property
    def minimum_detectable_frequency(self):
        if self.error_rate < 0.1:
            return 0.05
        else:
            return 0.25

    def _log_post_hom_ref(self, llk):
        return log(1) + llk

    def _log_post_het_or_alt(self, llk, expected_depth, sequence_coverage):
        expected_percentage_coverage = percent_coverage_from_expected_coverage(
            expected_depth)
        minimum_percentage_coverage_required = expected_percentage_coverage * \
            sequence_coverage.percent_coverage_threshold
        if sequence_coverage.percent_coverage > minimum_percentage_coverage_required:
            return self._log_post_hom_ref(llk)
        else:
            return MIN_LLK


class GeneCollectionTyper(Typer):

    """Types a collection of genes returning only the most likely version
        in the collection"""

    def __init__(
            self,
            expected_depths,
            contamination_depths=[],
            confidence_threshold=1):
        super(
            GeneCollectionTyper,
            self).__init__(
            expected_depths,
            contamination_depths,
            confidence_threshold=confidence_threshold)
        self.presence_typer = PresenceTyper(
            expected_depths, contamination_depths)

    def type(self, sequence_coverage_collection,
             min_gene_percent_covg_threshold=99):
        """Types a collection of genes returning the most likely gene version
            in the collection with it's genotype"""
        best_versions = self.get_best_version(
            sequence_coverage_collection.values(),
            min_gene_percent_covg_threshold)
        return [self.presence_typer.type(best_version)
                for best_version in best_versions]

    def get_best_version(
            self,
            sequence_coverages,
            min_gene_percent_covg_threshold):
        sequence_coverages=list(sequence_coverages)
        sequence_coverages.sort(key=lambda x: x.percent_coverage, reverse=True)
        current_best_gene = sequence_coverages[0]
        current_best_genes = [current_best_gene]
        for gene in sequence_coverages[1:]:
            # Report if above threshold regardless if it is
            if gene.percent_coverage >= min_gene_percent_covg_threshold:
                current_best_genes.append(gene)
            elif gene.percent_coverage < min_gene_percent_covg_threshold:
                return current_best_genes
        return current_best_genes
