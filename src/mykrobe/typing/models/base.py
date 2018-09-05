import json


class ProbeCoverage(object):

    "Summary of kmer coverage of sequence. e.g output of color covearges"

    def __init__(self, percent_coverage, median_depth, min_depth, k_count, klen):
        self.percent_coverage = percent_coverage
        self.median_depth = median_depth
        self.min_depth = min_depth
        self.k_count = k_count
        self.klen = klen

    @property
    def coverage_dict(self):
        return {"percent_coverage": round(self.percent_coverage, 2),
                "median_depth": round(self.median_depth, 2),
                "min_non_zero_depth": round(self.min_depth, 2),
                "kmer_count": self.k_count,
                "klen": self.klen,
                }

    def __str__(self):
        return json.dumps(self.coverage_dict)

    def __repr__(self):
        return json.dumps(self.coverage_dict)
