import datetime
import json
import logging
logger = logging.getLogger(__name__)


class VariantProbeCoverage(object):

    def __init__(self, reference_coverages,
                 alternate_coverages,
                 var_name=None,
                 params={}):
        self.reference_coverages = reference_coverages
        self.alternate_coverages = alternate_coverages
        self.var_name = var_name
        self.params = params
        if self.reference_coverages and self.alternate_coverages:
            self.best_alternate_coverage = self._choose_best_alternate_coverage()
            self.best_reference_coverage = self._choose_best_reference_coverage()

    def _choose_best_coverage(self, coverages):
        coverages.sort(
            key=lambda x: x.k_count,
            reverse=True)
        current_best = coverages[0]
        for probe_coverage in coverages[1:]:
            if probe_coverage.k_count < current_best.k_count:
                current_best = current_best
            else:
                if probe_coverage.percent_coverage > current_best.percent_coverage:
                    current_best = probe_coverage
                elif probe_coverage.min_depth > current_best.min_depth:
                    current_best = probe_coverage
                elif probe_coverage.min_depth <= current_best.min_depth:
                    if probe_coverage.median_depth > current_best.median_depth:
                        current_best = probe_coverage
        return current_best

    def _choose_best_alternate_coverage(self):
        return self._choose_best_coverage(self.alternate_coverages)

    def _choose_best_reference_coverage(self):
        best_reference_coverage = self._choose_best_coverage(
            self.reference_coverages)
        return best_reference_coverage

    @property
    def coverage_dict(self):
        return {"reference": self.best_reference_coverage.coverage_dict,
                "alternate": self.best_alternate_coverage.coverage_dict
                }

    def __str__(self):
        d = self.coverage_dict
        d['variant'] = self.var_name
        return json.dumps(d)

    def __repr__(self):
        return self.__str__()

    @property
    def reference_coverage(self):
        return self.best_reference_coverage

    @property
    def reference_percent_coverage(self):
        return self.best_reference_coverage.percent_coverage

    @property
    def reference_kmer_count(self):
        return self.best_reference_coverage.k_count

    @property
    def reference_median_depth(self):
        return self.best_reference_coverage.median_depth

    @property
    def reference_min_depth(self):
        return self.best_reference_coverage.min_depth

    @property
    def reference_klen(self):
        return self.best_reference_coverage.klen

    @property
    def alternate_percent_coverage(self):
        return self.best_alternate_coverage.percent_coverage

    @alternate_percent_coverage.setter
    def alternate_percent_coverage(self, value):
        self.best_alternate_coverage.percent_coverage = value

    @property
    def alternate_median_depth(self):
        return self.best_alternate_coverage.median_depth

    @property
    def alternate_kmer_count(self):
        return self.best_alternate_coverage.k_count

    @property
    def alternate_min_depth(self):
        return self.best_alternate_coverage.min_depth

    @property
    def alternate_klen(self):
        return self.best_alternate_coverage.klen        
