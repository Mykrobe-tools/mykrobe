

class SequenceProbeCoverage(object):

    def __init__(
            self,
            name,
            probe_coverage,
            percent_coverage_threshold=30,
            version=1,
            length=None):
        self.name = name
        self.probe_coverage = probe_coverage
        self.percent_coverage_threshold = percent_coverage_threshold
        self.version = version
        self.length = length

    @property
    def median_depth(self):
        return self.probe_coverage.median_depth

    @property
    def percent_coverage(self):
        return self.probe_coverage.percent_coverage

    @property
    def min_depth(self):
        return self.probe_coverage.min_depth

    @property
    def coverage_dict(self):
        return self.probe_coverage.coverage_dict
