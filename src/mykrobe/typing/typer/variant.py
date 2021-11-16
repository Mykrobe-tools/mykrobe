from mykrobe import K
from mykrobe.typing.typer.base import Typer
from mykrobe.stats import log_lik_R_S_coverage
from mykrobe.stats import log_lik_R_S_kmer_count
from mykrobe.typing.typer.base import MIN_LLK
from mykrobe.typing.typer.base import DEFAULT_MINOR_FREQ
from mykrobe.typing.typer.base import DEFAULT_ERROR_RATE


from mykrobe.stats import percent_coverage_from_expected_coverage
from mykrobe.stats import log_lik_probability_of_N_gaps
import logging
logger = logging.getLogger(__name__)


def likelihoods_to_confidence(l):
    if not len(l) > 1:
        raise ValueError(
            "Must have at least 2 likelihoods to calculate confidence")
    l_sorted = sorted(l, reverse=True)
    if l_sorted[2] == l[0]:
        # 0/1 and 1/1 are most likely - since haploid - conf compare with 0/0
        return int(round(l_sorted[0] - l[0]))
    else:
        # Otherwise, compare with 2 most likely
        return int(round(l_sorted[0] - l_sorted[1]))


class VariantTyper(Typer):

    def __init__(self, expected_depths, contamination_depths=None,
                 error_rate=DEFAULT_ERROR_RATE,
                 minor_freq=DEFAULT_MINOR_FREQ,
                 ignore_filtered=False,
                 filters=None,
                 confidence_threshold=3,
                 model="kmer_count",
                 kmer_size=K,
                 min_proportion_expected_depth=0.3,
                 ploidy="diploid"):
        if filters is None:
            filters = []
        if contamination_depths is None:
            contamination_depths = []
        super(
            VariantTyper,
            self).__init__(
            expected_depths,
            contamination_depths,
            error_rate,
            ignore_filtered=ignore_filtered,
            filters=filters,
            confidence_threshold=confidence_threshold)
        self.method = "MAP"
        self.error_rate = error_rate
        self.minor_freq = minor_freq
        self.kmer_size = kmer_size
        self.min_proportion_expected_depth = min_proportion_expected_depth
        self.ploidy = ploidy

        if model == "median_depth":
            self.model = DepthCoverageGenotypeModel(
                self.expected_depths, self.contamination_depths, self.error_rate, self.minor_freq, self.kmer_size)
        elif model == "kmer_count":
            # logger.debug("Genotyping using kc model")
            self.model = KmerCountGenotypeModel(
                self.expected_depths, self.contamination_depths, self.error_rate, self.minor_freq, self.kmer_size)
        self.ignore_filtered = ignore_filtered
        self.filters = filters

        if len(expected_depths) > 1:
            raise NotImplementedError("Mixed samples not handled yet")

    def type(self, variant_probe_coverages, variant=None):
        """
            Takes a list of VariantProbeCoverages and returns a Call for the Variant.
            Note, in the simplest case the list will be of length one. However, we may be typing the
            Variant on multiple backgrouds leading to multiple VariantProbes for a single Variant.

        """
        if not isinstance(variant_probe_coverages, list):
            variant_probe_coverages = [variant_probe_coverages]
        calls = []
        for variant_probe_coverage in variant_probe_coverages:
            calls.append(
                self._type_variant_probe_coverages(
                    variant_probe_coverage, variant_probe_coverage.var_name))
        hom_alt_calls = [c for c in calls if sum(c["genotype"]) > 1]
        het_calls = [c for c in calls if sum(c["genotype"]) == 1]
        if hom_alt_calls:
            hom_alt_calls.sort(key=lambda x: x["info"]["conf"], reverse=True)
            return hom_alt_calls[0]
        elif het_calls:
            het_calls.sort(key=lambda x: x["info"]["conf"], reverse=True)
            return het_calls[0]
        else:
            calls.sort(key=lambda x: x["info"]["conf"], reverse=True)
            return calls[0]

    @property
    def diploid_model(self):
        return self.ploidy=="diploid"        

    def _type_variant_probe_coverages(
            self, variant_probe_coverage, variant=None):
        hom_ref_likelihood = self.model.hom_ref_lik(variant_probe_coverage)
        hom_alt_likelihood = self.model.hom_alt_lik(variant_probe_coverage)
        if not self.has_contamination() and self.diploid_model:
            het_likelihood = self.model.het_lik(variant_probe_coverage)
        else:
            het_likelihood = MIN_LLK
        likelihoods = [hom_ref_likelihood, het_likelihood, hom_alt_likelihood]
        confidence = likelihoods_to_confidence(likelihoods)
        gt = self.likelihoods_to_genotype(
            likelihoods
        )
        info = {"coverage": variant_probe_coverage.coverage_dict,
                "expected_depths": self.expected_depths,
                "contamination_depths": self.contamination_depths,
                "filter": [],
                "conf": confidence}
        if gt == "-/-" and not self.ignore_filtered:
            if variant_probe_coverage.alternate_percent_coverage > variant_probe_coverage.reference_percent_coverage:
                gt = "1/1"
            else:
                gt = "0/0"
            if "MISSING_WT" in self.filters:
                info["filter"].append("MISSING_WT")
        elif "LOW_PERCENT_COVERAGE" in self.filters and variant_probe_coverage.alternate_percent_coverage < 100 and variant_probe_coverage.reference_percent_coverage < 100:
            info["filter"].append("LOW_PERCENT_COVERAGE")
            if self.ignore_filtered:
                gt = "0/0"
        if "LOW_GT_CONF" in self.filters and (confidence < self.confidence_threshold):
            info["filter"].append("LOW_GT_CONF")

        if "LOW_TOTAL_DEPTH" in self.filters:
            total_depth = variant_probe_coverage.reference_median_depth + variant_probe_coverage.alternate_median_depth
            expected_depth = self.expected_depths[0]
            if total_depth < self.min_proportion_expected_depth * expected_depth:
                logger.debug("%s" % variant_probe_coverage.var_name)
                info["filter"].append("LOW_TOTAL_DEPTH") 

        if not self.diploid_model:
            likelihoods=[likelihoods[0],likelihoods[2]]
        return {
            "variant": variant,
            "genotype": [
                int(i) for i in gt.split("/")],
            "genotype_likelihoods": likelihoods,
            "info": info,
            "_cls": "Call.VariantCall"}


class GenotypeModel(object):

    def __init__(self, expected_depths, contamination_depths, error_rate, minor_freq, kmer_size):
        self.expected_depths = expected_depths
        self.contamination_depths = contamination_depths
        self.error_rate = error_rate
        self.minor_freq = minor_freq
        self.kmer_size = kmer_size

    def hom_ref_lik(self, variant_probe_coverage):
        raise NotImplementedError

    def hom_alt_lik(self, variant_probe_coverage):
        raise NotImplementedError

    def het_lik(self, variant_probe_coverage):
        raise NotImplementedError


class KmerCountGenotypeModel(GenotypeModel):

    def __init__(self, expected_depths, contamination_depths, error_rate, minor_freq, kmer_size=K):
        super(KmerCountGenotypeModel, self).__init__(
            expected_depths, contamination_depths, error_rate, minor_freq,kmer_size)

    def depth_to_expected_kmer_count(self, depth, allele_length):
        # logger.debug("%f %f"% (allele_length,depth))
        return (allele_length*depth)+0.01

    def klen_to_dna_len(self, klen):
        return klen+self.kmer_size-1

    def hom_ref_lik(self, variant_probe_coverage):
        hom_ref_likes = []
        # Either alt+cov or alt_covg + contam_covg
        for expected_depth in self.expected_depths:
            # logger.debug("hom ref %s" % variant_probe_coverage.var_name)
            kmer_count_likelihood = log_lik_R_S_kmer_count(
                    variant_probe_coverage.reference_kmer_count,
                    variant_probe_coverage.alternate_kmer_count,
                    self.depth_to_expected_kmer_count(expected_depth, variant_probe_coverage.reference_klen),
                    self.depth_to_expected_kmer_count(expected_depth * self.error_rate/3, variant_probe_coverage.alternate_klen)
                    )
            ngaps_likelihood=log_lik_probability_of_N_gaps(expected_depth,
                                              variant_probe_coverage.reference_percent_coverage,
                                              self.klen_to_dna_len(variant_probe_coverage.reference_klen))+\
                             log_lik_probability_of_N_gaps(expected_depth * self.error_rate/3,
                                              variant_probe_coverage.alternate_percent_coverage,
                                              self.klen_to_dna_len(variant_probe_coverage.alternate_klen))
            logger.debug("hom ref kmer_count_likelihood: %f %s" % (kmer_count_likelihood,variant_probe_coverage.var_name))
            logger.debug("hom ref ngaps_likelihood: %f %s " % (ngaps_likelihood,variant_probe_coverage.var_name))

            hom_ref_likes.append(kmer_count_likelihood + ngaps_likelihood)
            # for contamination in self.contamination_depths:
            #     hom_ref_likes.append(
            #         log_lik_R_S_kmer_count(
            #             variant_probe_coverage.reference_kmer_count,
            #             variant_probe_coverage.alternate_kmer_count,
            #             self.depth_to_expected_kmer_count(expected_depth + contamination,variant_probe_coverage.reference_klen),
            #             self.depth_to_expected_kmer_count((expected_depth + contamination) * self.error_rate / 3),variant_probe_coverage.alternate_klen))
        return max(hom_ref_likes)

    def hom_alt_lik(self, variant_probe_coverage):
        hom_alt_liks = []
        # Either alt+cov or alt_covg + contam_covg
        # logger.debug("hom alt %s" % variant_probe_coverage.var_name)
        for expected_depth in self.expected_depths:
            kmer_count_likelihood = log_lik_R_S_kmer_count(
                    variant_probe_coverage.alternate_kmer_count,
                    variant_probe_coverage.reference_kmer_count,
                    self.depth_to_expected_kmer_count(expected_depth,variant_probe_coverage.reference_klen),
                    self.depth_to_expected_kmer_count(expected_depth * self.error_rate / 3,variant_probe_coverage.alternate_klen)
                    )
            ngaps_likelihood=log_lik_probability_of_N_gaps(expected_depth * self.error_rate / 3,
                                              variant_probe_coverage.reference_percent_coverage,
                                              self.klen_to_dna_len(variant_probe_coverage.reference_klen))+\
                             log_lik_probability_of_N_gaps(expected_depth,
                                              variant_probe_coverage.alternate_percent_coverage,
                                              self.klen_to_dna_len(variant_probe_coverage.alternate_klen))
            logger.debug("hom alt kmer_count_likelihood: %f %s" % (kmer_count_likelihood,variant_probe_coverage.var_name))
            logger.debug("hom alt ngaps_likelihood: %f %s " % (ngaps_likelihood,variant_probe_coverage.var_name))
                
            hom_alt_liks.append(kmer_count_likelihood + ngaps_likelihood)
            # for contamination in self.contamination_depths:
            #     hom_alt_liks.append(
            #         log_lik_R_S_kmer_count(
            #             variant_probe_coverage.alternate_kmer_count,
            #             variant_probe_coverage.reference_kmer_count,
            #             self.depth_to_expected_kmer_count(expected_depth + contamination,variant_probe_coverage.alternate_klen),
            #             self.depth_to_expected_kmer_count((expected_depth + contamination) * self.error_rate / 3,variant_probe_coverage.reference_klen)))
        return max(hom_alt_liks)

    def het_lik(self, variant_probe_coverage):
        # logger.debug("het %s" % variant_probe_coverage.var_name)

        if (variant_probe_coverage.alternate_kmer_count+variant_probe_coverage.reference_kmer_count) == 0:
            return MIN_LLK
        elif variant_probe_coverage.alternate_percent_coverage < 100 or variant_probe_coverage.reference_percent_coverage < 100:
            return MIN_LLK
        else:
            het_liks = []
            for expected_depth in self.expected_depths:
                kmer_count_likelihood = log_lik_R_S_kmer_count(
                        variant_probe_coverage.alternate_kmer_count,
                        variant_probe_coverage.reference_kmer_count,
                        self.depth_to_expected_kmer_count(expected_depth/2 + (expected_depth/2 * self.error_rate/3),variant_probe_coverage.alternate_klen),
                        self.depth_to_expected_kmer_count(expected_depth/2 + (expected_depth/2 * self.error_rate/3),variant_probe_coverage.reference_klen)
                        )
                ngaps_likelihood_ref=log_lik_probability_of_N_gaps(expected_depth/2 + (expected_depth/2 * self.error_rate/3),
                                                  variant_probe_coverage.reference_percent_coverage,
                                                  self.klen_to_dna_len(variant_probe_coverage.reference_klen))
                ngaps_likelihood_alt=log_lik_probability_of_N_gaps(expected_depth/2 + (expected_depth/2 * self.error_rate/3),
                                                  variant_probe_coverage.alternate_percent_coverage,
                                                  self.klen_to_dna_len(variant_probe_coverage.alternate_klen))
                ngaps_likelihood=ngaps_likelihood_ref+ngaps_likelihood_alt
                logger.debug("het kmer_count_likelihood: %f %s" % (kmer_count_likelihood,variant_probe_coverage.var_name))
                logger.debug("het ngaps_likelihood_ref: %f %s " % (ngaps_likelihood_ref,variant_probe_coverage.var_name))
                logger.debug("het ngaps_likelihood_alt: %f %s " % (ngaps_likelihood_alt,variant_probe_coverage.var_name))
                logger.debug("het ngaps_likelihood: %f %s " % (ngaps_likelihood,variant_probe_coverage.var_name))


                het_liks.append(kmer_count_likelihood + ngaps_likelihood)
            return max(het_liks)


class DepthCoverageGenotypeModel(GenotypeModel):

    def __init__(self, expected_depths, contamination_depths, error_rate, minor_freq, kmer_size=K):
        super(DepthCoverageGenotypeModel, self).__init__(
            expected_depths, contamination_depths, error_rate, minor_freq, kmer_size)

    def hom_ref_lik(self, variant_probe_coverage):
        if variant_probe_coverage.reference_percent_coverage < 100 * \
                percent_coverage_from_expected_coverage(max(self.expected_depths)):
            return MIN_LLK
        else:
            hom_ref_likes = []
            # Either alt+cov or alt_covg + contam_covg
            for expected_depth in self.expected_depths:
                hom_ref_likes.append(
                    log_lik_R_S_coverage(
                        variant_probe_coverage.reference_median_depth,
                        variant_probe_coverage.alternate_median_depth,
                        expected_depth,
                        expected_depth *
                        self.error_rate /
                        3))
                for contamination in self.contamination_depths:
                    hom_ref_likes.append(
                        log_lik_R_S_coverage(
                            variant_probe_coverage.reference_median_depth,
                            variant_probe_coverage.alternate_median_depth,
                            expected_depth + contamination,
                            (expected_depth + contamination) * self.error_rate / 3))
            return max(hom_ref_likes)

    def hom_alt_lik(self, variant_probe_coverage):
        if variant_probe_coverage.alternate_percent_coverage < 100 * \
                percent_coverage_from_expected_coverage(max(self.expected_depths)):
            return MIN_LLK
        else:
            hom_alt_liks = []
            # Either alt+cov or alt_covg + contam_covg
            for expected_depth in self.expected_depths:
                hom_alt_liks.append(
                    log_lik_R_S_coverage(
                        variant_probe_coverage.alternate_median_depth,
                        variant_probe_coverage.reference_median_depth,
                        expected_depth,
                        expected_depth *
                        self.error_rate /
                        3))
                for contamination in self.contamination_depths:
                    hom_alt_liks.append(
                        log_lik_R_S_coverage(
                            variant_probe_coverage.alternate_median_depth,
                            variant_probe_coverage.reference_median_depth,
                            expected_depth + contamination,
                            (expected_depth + contamination) * self.error_rate / 3))
            return max(hom_alt_liks)

    def het_lik(self, variant_probe_coverage):
        if variant_probe_coverage.alternate_percent_coverage < 100 or variant_probe_coverage.reference_percent_coverage < 100:
            return MIN_LLK
        else:
            het_liks = []
            for expected_depth in self.expected_depths:
                het_liks.append(
                    log_lik_R_S_coverage(
                        variant_probe_coverage.alternate_median_depth,
                        variant_probe_coverage.reference_median_depth,
                        expected_depth * self.minor_freq,
                        expected_depth * (
                            1 - self.minor_freq)))
            return max(het_liks)
